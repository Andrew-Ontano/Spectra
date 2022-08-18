import pandas
import ruptures as rpt
import pandas as pd
import numpy as np
from Bio import Seq

# Shorthand for string reverse-complement using Biopython
def rc(sequence):
    sequence = Seq.Seq(sequence)
    return str(sequence.reverse_complement())

#
def simplify(spectra, index=4, dim=64):
    simpleQueries = {}
    for query in list(spectra.columns)[index:index + dim]:
        if query not in list(simpleQueries.values()) or rc(query) not in list(simpleQueries.keys()):
            simpleQueries[query] = rc(query)
    for mer in simpleQueries:
        spectra[mer] += spectra[simpleQueries[mer]]
    return spectra.drop(columns=list(simpleQueries.values()))

# Counts the spectra of a sequence
def windowCount(seq):
    return seq[4] + [seq[2] + 1, seq[2] + len(seq[0]) if seq[2] + len(seq[0]) < seq[3] else seq[3]] + [
               seq[0].count_overlap(a) for a in seq[1]]

# Calculate breakpoints from spectra
# For 64 literal counts, a penalty of 1,000,000 is ideal, but for frequencies a penalty of 0.5 is ideal
def getBreakpoints(spectra, index=4, dim=64, penalty=1000000, min_size=5):
    spectra = spectra.groupby(['Library', 'Sequence'])
    output = []
    # for each grouping, calculate rupture breakpoints, then convert them to starting window indices
    for group in spectra:
        data = group[1].iloc[0:len(group[1]), index:index + dim].to_numpy()
        if len(data) > min_size * 2:
            dataAlgo = rpt.KernelCPD(min_size=min_size).fit(data).predict(pen=penalty)
            output.append((group[0], group[1].iloc[[a-1 for a in dataAlgo], 2:3]['Start'].tolist()))
        else:
            output.append((group[0], []))
    return output

# Use breakpoints to append as bin column to spectra
def applyBreakpoints(spectra, breakpoints):
    for bkp in breakpoints:
        breakCount = 0
        index = 1
        for bkpStart in bkp[1]:
            breakName = f'{bkp[0][0]}_{bkp[0][1]}_{breakCount}'
            indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                    spectra['Start'] >= index) & (spectra['End'] < bkpStart)].index)
            spectra.loc[indices, 'Bin'] = breakName
            index = bkpStart
            breakCount += 1
        indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                    spectra['Start'] >= index)].index)
        spectra.loc[indices, 'Bin'] = breakName
    return spectra

def getBreakpointFrequencies(spectra, frequency, index=4, dim=64):
    spectra = spectra.groupby(['Bin', 'Library', 'Sequence'])
    outputs = pd.DataFrame()
    for group in spectra:
        frequencies = getGlobalFrequencies(group[1], frequency, index=index, dim=dim)
        outputs.loc[len(outputs.index), ['Library', 'Sequence', 'Bin', 'Length', 'Start', 'End'] + list(frequencies.keys())] = [group[0][1], group[0][2], str(group[0][0]), str(max(group[1]['End']) - min(group[1]['Start']) + 1), str(min(group[1]['Start'])), str(max(group[1]['End']))] + list(frequencies.values())
    return outputs

# Transform spectra counts to spectra frequencies
def countToFrequency(spectra, index=4, dim=64, len=3):
    for row in spectra.iterrows():
        denominator = row[1]['End'] - row[1]['Start'] - (len - 2)
        spectra.iloc[row[0], index:index + dim] = np.array(spectra.iloc[row[0], index:index + dim]) / denominator
    return spectra

# Transform spectra counts to spectra frequencies
def frequencyToCount(spectra, index=4, dim=64, len=3):
    return

# Calculate global frequencies across spectra
def getGlobalFrequencies(spectra, frequency=False, index=4, dim=64):
    if frequency:
        return dict(zip(spectra.columns[index:index + dim], np.array(spectra.iloc[0:len(spectra), index:index + dim].sum()/len(spectra))))
    else:
        counts = np.array(spectra.sum())
        return dict(zip(spectra.columns[index:index + dim], np.divide(counts[index:index + dim], counts[3] - counts[2] - len(spectra.Library))))

# Reduce frequencies from a given set of global frequencies
def reduceFrequencies():
    return
