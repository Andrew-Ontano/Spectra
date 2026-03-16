#!/usr/bin/env python3
import argparse
import math
import tempfile
import heapq
import os
import logging

def read_dump(filename):
    """Generator to read jellyfish dump file and yield (kmer, count) tuples."""
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                yield parts[0], int(parts[1])

def stream_compute_reductions(raw_file, asm_file, out_file, chunk_size=5000000, extreme=None):
    """
    Stream merge raw and asm files, compute log-fold reductions,
    write them to disk in sorted chunks, then external-merge them into ranked output.
    """
    logger.info(f"Reading {raw_file} and {asm_file}...")
    tmp_files = []
    total_matched = 0

    # Pass 1: stream merge + chunked sort
    raw_gen = read_dump(raw_file)
    asm_gen = read_dump(asm_file)

    try:
        r_kmer, r_count = next(raw_gen)
    except StopIteration:
        r_kmer = None
    try:
        a_kmer, a_count = next(asm_gen)
    except StopIteration:
        a_kmer = None

    buffer = []

    while r_kmer is not None and a_kmer is not None:
        if r_kmer == a_kmer:
            reduction = math.log10(a_count + 1) - math.log10(r_count + 1)
            buffer.append((reduction, r_kmer, r_count, a_count))
            total_matched += 1

            if len(buffer) >= chunk_size:
                buffer.sort(key=lambda x: x[0])
                tf = tempfile.NamedTemporaryFile(delete=False, mode="w")
                for red, k, rcount, acount in buffer:
                    tf.write(f"{red}\t{k}\t{rcount}\t{acount}\n")
                tf.close()
                tmp_files.append(tf.name)
                buffer.clear()
                logger.info(f"Processed {total_matched} matches...")

            try:
                r_kmer, r_count = next(raw_gen)
            except StopIteration:
                r_kmer = None
            try:
                a_kmer, a_count = next(asm_gen)
            except StopIteration:
                a_kmer = None

        elif r_kmer < a_kmer:
            try:
                r_kmer, r_count = next(raw_gen)
            except StopIteration:
                r_kmer = None
        else:
            try:
                a_kmer, a_count = next(asm_gen)
            except StopIteration:
                a_kmer = None

    # flush last buffer
    if buffer:
        buffer.sort(key=lambda x: x[0])
        tf = tempfile.NamedTemporaryFile(delete=False, mode="w")
        for red, k, rcount, acount in buffer:
            tf.write(f"{red}\t{k}\t{rcount}\t{acount}\n")
        tf.close()
        tmp_files.append(tf.name)
        buffer.clear()

    # Pass 2: external merge of sorted chunks
    logger.info(f"Merging {len(tmp_files)} sorted chunks...")
    def file_iter(fname):
        with open(fname) as f:
            for line in f:
                red, k, rc, ac = line.strip().split("\t")
                yield float(red), k, int(rc), int(ac)

    iterators = [file_iter(f) for f in tmp_files]
    merged = heapq.merge(*iterators, key=lambda x: x[0])

    # Pass 3: write ranked table
    logger.info(f"Writing ranked table to {out_file}...")
    with open(out_file, "w") as out:
        out.write("kmer\tRaw\tAsm\treduction\treductionRank\n")

        num_to_keep = None
        if extreme is not None:
            num_to_keep = int(total_matched * extreme / 100.0)
            logger.info(f"Filtering for extreme {extreme}%: keeping top and bottom {num_to_keep} k-mers.")

        rank = 1
        for red, k, rc, ac in merged:
            if num_to_keep is None or rank <= num_to_keep or rank > (total_matched - num_to_keep):
                out.write(f"{k}\t{rc}\t{ac}\t{red:.6f}\t{rank}\n")
            rank += 1

    # cleanup
    for f in tmp_files:
        os.remove(f)


# CLI arguments
parser = argparse.ArgumentParser(description="Rank kmers by log-fold change between raw and assembly")
parser.add_argument("-r", "--raw_dump", required=True, help="Raw jellyfish dump file")
parser.add_argument("-a", "--asm_dump", required=True, help="Assembly jellyfish dump file")
parser.add_argument("-o", "--output", required=True, help="Output ranked table (TSV)")
parser.add_argument("-c", "--chunk_size", type=int, default=5000000,
                    help="Number of kmers to hold in memory before spilling to disk [default 5,000,000]")
parser.add_argument("-e", "--extreme", type=float, default=None,
                    help="Keep only the top/bottom PERCENT of reductions (e.g. 5 = keep 5%% lowest and 5%% highest)")
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

# Logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
if args.verbose:
    logger.setLevel(logging.INFO)

stream_compute_reductions(args.raw_dump, args.asm_dump, args.output, args.chunk_size, extreme=args.extreme)

logger.info(f"Ranked table written to {args.output}")
