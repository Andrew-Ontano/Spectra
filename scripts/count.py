#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral
import pandas as pd
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    maxSize = 30000000 if 30000000 % args.spacing != 0 else (30000000 // args.spacing) * args.spacing

    if args.verbose:
        logger.setLevel(logging.INFO)

    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

    sequences = {}
    if args.memory:
        try:
            sequences = SeqIO.index(args.input_sequence, args.sequence_format)
            if len(sequences.keys()) == 0:
                logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()
    else:
        try:
            sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
            if len(sequences.keys()) == 0:
                logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()
    bases = ["A", "C", "G", "T"]
    queries = [f"{a}{b}{c}" for c in bases for b in bases for a in bases]

    with open(args.output, 'w', newline='') as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter='\t')
        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)
        callableProcess = spectral.windowCount if args.overlap else spectral.windowCountNoOverlap

        for sequence in sequences:
            headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
            sequenceLength = len(sequences[sequence])
            if args.memory or sequenceLength >= maxSize:
                logging.info(f"Sequence {sequence} is large. Breaking into smaller segments")
                indices = range(0, sequenceLength, maxSize)
                for sequenceIndex in indices:
                    toProcess = [[sequences[sequence][i:i + args.width].seq.upper(), queries, i, i + args.width, headers] for i in range(sequenceIndex, sequenceIndex + maxSize, args.spacing) if len(sequences[sequence][i:i + args.width]) > 0]
                    rows = map(callableProcess, toProcess)
                    tsvWriter.writerows(rows)
                logging.info(f"Sequence {sequence} windows written to output file")
            else:
                toProcess = [[sequences[sequence][i:i+args.width].seq.upper(), queries, i, i + args.width, headers] for i in range(0, sequenceLength+args.spacing, args.spacing) if len(sequences[sequence][i:i+args.width]) > 0]
                rows = map(callableProcess, toProcess)
                tsvWriter.writerows(rows)
                logging.info(f"Sequence {sequence} windows written to output file")

    if args.complement:
        logging.info("Simplifying forward and r-c counts")
        spectra = pd.read_csv(args.output, delimiter='\t')
        spectra = spectral.simplify(spectra, dim=len(queries))
        spectra.to_csv(args.output, sep='\t', index=False)

    logging.info(f'Execution time in seconds: {time.time() - startTime}')
