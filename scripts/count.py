#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral
import multiprocessing
from collections import namedtuple

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])

def window_tasks_generator(input_sequence, sequence_format, queries_all, width, spacing, libraries, minimum_size, chunk_size):
    for seq_record in SeqIO.parse(input_sequence, sequence_format):
        sequence_name = seq_record.id
        headers = sequence_name.split("_") if libraries else [os.path.basename(input_sequence), sequence_name]
        sequenceLength = len(seq_record)

        if sequenceLength < minimum_size:
            continue

        # Process in chunks to save memory
        for chunk_start in range(0, sequenceLength, chunk_size):
            chunk_end = min(chunk_start + chunk_size, sequenceLength)
            # Ensure chunk includes enough overlap for windows starting near the end of the chunk
            # but here width and spacing are used. Actually, for sliding windows,
            # we should be careful.
            # If we just take the chunk and process windows in it:
            sub_seq = str(seq_record.seq[chunk_start:chunk_end]).upper()

            for i in range(0, len(sub_seq), spacing):
                window = sub_seq[i:i+width]
                if not window or len(window) < width and (chunk_start + i + len(window) < sequenceLength):
                    # Skip incomplete windows unless it's the very end of the sequence
                    continue

                start = chunk_start + i
                end = start + len(window)
                yield WindowTask(window, queries_all, start, end, headers)

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

    # Pre-calculate queries
    queries_all = spectral.setMers(args.mer_size)
    if args.complement:
        queries = spectral.mapCanonicalMers(queries_all)
    else:
        queries = {queries_all[a]: [a] for a in range(len(queries_all))}

    # Prepare for output
    tsvHeaders = ["Library", "Sequence", "Start", "End"] + list(queries.keys())

    # Process windows in parallel
    callableProcess = spectral.windowCount if args.overlap else spectral.windowCountNoOverlap

    with open(args.output, 'w', newline='') as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter='\t')
        tsvWriter.writerow(tsvHeaders)

        tasks = window_tasks_generator(
            args.input_sequence,
            args.sequence_format,
            queries_all,
            args.width,
            args.spacing,
            args.libraries,
            args.minimum_size,
            args.chunk_size
        )

        if args.threads > 1:
            pool = multiprocessing.Pool(processes=args.threads)
            # Use imap to maintain order and be memory efficient
            for row in pool.imap(callableProcess, tasks):
                if args.complement:
                    tsvWriter.writerow(spectral.collapseRC(row, queries, dim=len(queries_all)))
                else:
                    tsvWriter.writerow(row)
            pool.close()
            pool.join()
        else:
            for task in tasks:
                row = callableProcess(task)
                if args.complement:
                    tsvWriter.writerow(spectral.collapseRC(row, queries, dim=len(queries_all)))
                else:
                    tsvWriter.writerow(row)

    logging.info(f'Execution time in seconds: {time.time() - startTime}')
