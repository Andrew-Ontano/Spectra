#!/usr/bin/env python3

# provide info about a given spectra file

import ruptures as rpt
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import spectral

logging.basicConfig(level=logging.INFO)

# plots breakpoints calculated through ruptures
def plotBreakpoints(sequence, data, dataAlgo):
    fig, ax_arr = rpt.display(data, dataAlgo, figsize=(30, 150))
    plt.savefig(f'breakpoints_breakdown_{sequence}.png')
    plt.close()

# function to iterate over the windows and convert them to frequencies
# raw counts misrepresent kmer diversity in gappy alignments
def padWindows(spectra, tally, mer):
    for index, row in spectra.iterrows():
        spectra.loc[index, mer] = [a/tally[index] for a in row[mer]]
    return spectra

def execute(args):
    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    indexLength = 4
    spectraDimensions = len(spectra.columns) - indexLength
    if args.is_binned:
        spectraDimensions -= 1
        results = spectral.getBreakpointFrequencies(spectra, args.frequency, index=indexLength, dim=spectraDimensions)
        if args.output_tsv:
            results.to_csv(args.output_tsv, index=False, sep='\t')
        else:
            print(results)
        exit()
    if args.frequency:
        # penalty is temporarily locked to 1 or lower for frequency-based breakpoints
        if args.penalty > 1:
            args.penalty = 0.5
        spectra = spectral.countToFrequency(spectra, index=indexLength, dim=spectraDimensions)
    breakpoints = spectral.getBreakpoints(spectra, penalty=args.penalty, index=indexLength, dim=spectraDimensions)
    spectra = spectral.applyBreakpoints(spectra, breakpoints)

    if args.output_tsv:
        spectra.to_csv(args.output_tsv, sep='\t', index=False)
    else:
        print(spectra)


