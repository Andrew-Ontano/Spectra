#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral
import pandas as pd
import multiprocessing
from collections import namedtuple

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])

def window_tasks_generator(sequences, queries, width, spacing, libraries, input_sequence):
    for sequence_id in sequences:
        headers = sequence_id.split("_") if libraries else [os.path.basename(input_sequence), sequence_id]
        seq_record = sequences[sequence_id]
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        for i in range(0, seq_len, spacing):
            window = seq_str[i:i+width]
            if not window:
                continue
            yield WindowTask(window, queries, i, min(i+width, seq_len), headers)

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

    if "," in args.query:
        queries = [a.upper() for a in args.query.split(',')]
    else:
        queries = [args.query.upper()]

    try:
        if args.memory:
            sequences = SeqIO.index(args.input_sequence, args.sequence_format)
        else:
            sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))

        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
        exit()

    if args.complement:
        newQueries = []
        for query in queries:
            queryRC = spectral.rc(query)
            if query not in newQueries:
                newQueries.append(query)
            if queryRC not in newQueries:
                newQueries.append(queryRC)
        queries = newQueries

    with open(args.output, 'w', newline='') as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter='\t')
        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

        callableProcess = spectral.windowCount if args.overlap else spectral.windowCountNoOverlap
        tasks = window_tasks_generator(sequences, queries, args.width, args.spacing, args.libraries, args.input_sequence)

        if args.threads > 1:
            pool = multiprocessing.Pool(processes=args.threads)
            for row in pool.imap(callableProcess, tasks):
                tsvWriter.writerow(row)
            pool.close()
            pool.join()
        else:
            for task in tasks:
                tsvWriter.writerow(callableProcess(task))

    if args.complement:
        logging.info("Simplifying forward and r-c counts")
        spectra = pd.read_csv(args.output, delimiter='\t')
        spectra = spectral.simplify(spectra, dim=len(queries))
        spectra.to_csv(args.output, sep='\t', index=False)

    logging.info(f'Execution time in seconds: {time.time() - startTime}')
