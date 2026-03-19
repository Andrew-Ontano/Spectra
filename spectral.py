import ruptures as rpt
import pandas as pd
import numpy as np
from Bio import Seq
from collections import Counter
import itertools

# Fast reverse complement using translation table
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")

def setMers(merSize=3):
    bases = ["A", "C", "G", "T"]
    return ["".join(a) for a in list(itertools.product(bases, repeat=merSize))]

# Shorthand for string reverse-complement
def rc(sequence):
    if isinstance(sequence, Seq.Seq):
        sequence = str(sequence)
    return sequence.translate(RC_TRANS)[::-1]

def canonical(kmer):
    return min(kmer, rc(kmer))

def mapCanonicalMers(kmers):
    kmer_index = {k: i for i, k in enumerate(kmers)}
    seen = set()
    result = {}
    for kmer, i in kmer_index.items():
        rcMer = rc(kmer)
        if kmer in seen or rcMer in seen:
            continue
        j = kmer_index[rcMer]
        # canonical = lexicographically smaller of the pair
        result[canonical(kmer)] = [i, j]
        seen.add(kmer)
        seen.add(rcMer)
    return result

def collapseRC(row, mers, index=4, dim=None):
    if dim is None:
        dim = len(row) - index
    meta = row[:index]
    counts = row[index:index+dim]
    collapsedCounts = []
    for canon, (i, j) in mers.items():
        collapsedCounts.append(counts[i] + counts[j])
    return meta + collapsedCounts

# Validate takes an input df, and returns a tuple of (0) approved spectra descriptor fields and (1) triplet or query names
def validate(df):
    columns = df.columns
    validatedDescriptors = []
    for descriptor in ['Library', 'Sequence', 'Block', 'Length', 'Start', 'End']:
        if descriptor in columns:
            validatedDescriptors.append(descriptor)
    validatedNames = [a for a in columns if a[0] in ['A', 'C', 'G', 'T']]
    return validatedDescriptors, validatedNames

# Simplify bi-directional frequencies/counts to unidirectional
def simplify(spectra, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index
    simpleQueries = {}
    for query in list(spectra.columns)[index:index + dim]:
        if query not in list(simpleQueries.values()) or rc(query) not in list(simpleQueries.keys()):
            simpleQueries[query] = rc(query)
    for mer in simpleQueries:
        if simpleQueries[mer] != mer:
            spectra[mer] += spectra[simpleQueries[mer]]
    return spectra.drop(columns=[simpleQueries[a] for a in simpleQueries if a != simpleQueries[a]])

# Counts the spectra of a sequence
def windowCount(seq):
    windowSeq, queries, start, end, headers = seq
    windowSeq = str(windowSeq)
    merLen = len(queries[0]) if queries else 3
    counts = Counter(windowSeq[i:i + merLen] for i in range(len(windowSeq) - merLen + 1))
    return headers + [start + 1, min(start + len(windowSeq), end)] + [counts.get(q, 0) for q in queries]

def windowCountNoOverlap(seq):
    windowSeq, queries, start, end, headers = seq
    windowSeq = str(windowSeq)
    return headers + [start + 1, min(start + len(windowSeq), end)] + [windowSeq.count(q) for q in queries]

# Calculate breakpoints from spectra
# For literal counts, a high penalty is ideal, but for frequencies a lower penalty is ideal
def getBreakpoints(spectra, index=4, dim=None, penalty=1000000, min_size=5):
    if dim is None:
        dim = len(spectra.columns) - index
    spectra_grouped = spectra.groupby(['Library', 'Sequence'])
    output = []
    # for each grouping, calculate rupture breakpoints, then convert them to starting window indices
    for group in spectra_grouped:
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
        breakName = ''
        if bkp[1]:
            for bkpStart in bkp[1]:
                breakName = f'{bkp[0][0]}_{bkp[0][1]}_{breakCount:02d}'
                indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                        spectra['Start'] >= index) & (spectra['End'] < bkpStart)].index)
                spectra.loc[indices, 'Block'] = breakName
                index = bkpStart
                breakCount += 1
            indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                        spectra['Start'] >= index)].index)
            spectra.loc[indices, 'Block'] = breakName
        else:
            breakName = f'{bkp[0][0]}_{bkp[0][1]}_{breakCount:02d}'
            indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1])].index)
            spectra.loc[indices, 'Block'] = breakName
    return spectra

def getBreakpointFrequencies(spectra, frequency, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index
    spectra_grouped = spectra.groupby(['Block', 'Library', 'Sequence'])
    outputs = pd.DataFrame()
    for group in spectra_grouped:
        frequencies = getGlobalFrequencies(group[1], frequency, index=index, dim=dim)
        outputs.loc[len(outputs.index), ['Library', 'Sequence', 'Bin', 'Length', 'Start', 'End'] + list(frequencies.keys())] = [group[0][1], group[0][2], str(group[0][0]), str(max(group[1]['End']) - min(group[1]['Start']) + 1), str(min(group[1]['Start'])), str(max(group[1]['End']))] + list(frequencies.values())
    return outputs

# Transform spectra counts to spectra frequencies
def countToFrequency(spectra, index=4, dim=None, merLen=None):
    if dim is None:
        dim = len(spectra.columns) - index
    if merLen is None:
        merLen = len(spectra.columns[index])
    for column in range(index, index + dim):
        spectra[spectra.columns[column]] = spectra[spectra.columns[column]].astype("float")
    for rowIndex, row in spectra.iterrows():
        denominator = row['End'] - row['Start'] - (merLen - 2)
        if denominator > 0:
            spectra.iloc[rowIndex, index:index + dim] = spectra.iloc[rowIndex, index:index + dim] / denominator
        else:
            spectra.iloc[rowIndex, index:index + dim] = 0
    return spectra

# Transform spectra counts to spectra frequencies (Placeholder)
def frequencyToCount(spectra, index=4, dim=None, merLen=None):
    return spectra

# Calculate global frequencies across spectra
def getGlobalFrequencies(spectra, frequency=False, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index
    if frequency:
        return dict(zip(spectra.columns[index:index + dim], np.array(spectra.iloc[0:len(spectra), index:index + dim].sum()/len(spectra))))
    else:
        counts = np.array(spectra.sum())
        # denominator calculation might need to be generalized as well
        # Original: counts[3] - counts[2] - len(spectra.Library)
        # 3 is End, 2 is Start. Correct.
        # len(spectra.Library) is number of windows. Correct for 3-mers (L - 2).
        # General: L - (merLen - 1).
        # If we assume 3-mer it's L - 2.
        # Let's try to get merLen.
        merLen = len(spectra.columns[index])
        denominator = counts[3] - counts[2] - (merLen - 2) * len(spectra)
        return dict(zip(spectra.columns[index:index + dim], np.divide(counts[index:index + dim], denominator)))

# Reduce frequencies from a given set of global frequencies
def reduceFrequencies():
    return

# needs work. Currently does nothing
def spectraRC(spectra, index=4, dim=None, merLen=None):
    if dim is None:
        dim = len(spectra.columns) - index
    if merLen is None:
        merLen = len(spectra.columns[index])
    newSpectra = spectra.copy()
    mers = list(spectra.columns)[index:index + dim]
    for mer in mers:
        rc_mer = rc(mer)
        if rc_mer in spectra.columns:
            newSpectra[rc_mer] = spectra[mer]
    return newSpectra
