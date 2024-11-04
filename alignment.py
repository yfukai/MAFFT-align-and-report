import pandas as pd
import numpy as np
import itertools
from Bio import Align, Seq
from Bio.Align import substitution_matrices
from os import path
from scipy.spatial import distance

nuc_codes = {
   "A":["A"],
    "T":["T"],
    "G":["G"],
    "C":["C"],
    "R":["A","G"],
    "Y":["C","T"],
    "S":["G","C"],
    "W":["A","T"],
    "K":["G","T"],
    "M":["A","C"],
    "B":["C","G","T"],
    "D":["A","G","T"],
    "H":["A","C","T"],
    "V":["A","C","G"],
    "N":["A","C","G","T"],
}

match_score = 1
shared_score = 0.8
mismatch_score = -1
letters = "".join(nuc_codes.keys())
#print(letters)

score_matrix = np.zeros((len(letters), len(letters)))

for (j1, l1), (j2, l2) in itertools.product(enumerate(letters), repeat=2):
    if l1 == l2:
        if len(l1) == 1:
            score_matrix[j1, j2] = match_score
        else:
            score_matrix[j1, j2] = shared_score
    elif len(set(nuc_codes[l1]) & set(nuc_codes[l2])) > 0:
        score_matrix[j1, j2] = shared_score
    else:
        score_matrix[j1, j2] = mismatch_score

def get_aligner():
    custom_substitution_matrix = substitution_matrices.Array(
        letters, 2, score_matrix
    )
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -2
    aligner.extend_gap_score = - 4
    aligner.substitution_matrix = custom_substitution_matrix
    aligner.query_left_open_gap_score = - 10
    aligner.query_left_extend_gap_score = - 10
    aligner.query_right_open_gap_score = - 10
    aligner.query_right_extend_gap_score = - 10
    aligner.mode = 'local'
    #assert Seq.Seq("M").reverse_complement() == Seq.Seq("K")
    return aligner


def alignment_to_sequence(alignment,start,end):
    """
    indices ... 
    Return the sequence index of each lettter in the alignment.
    This property returns a 2D NumPy array with the sequence index of each letter in the alignment. 
    Gaps are indicated by -1. 
    The array has the same number of rows and columns as the alignment, as given by self.shape.    """
    val_dict = {k:v for v, k in zip(*alignment.indices)} # v is the query index, k is the target index
    target_positions = [val_dict.get(p,-1)
                        for p in range(start,end)]
    target_positions = [tp if tp<len(alignment.target) else -1 for tp in target_positions]
    
    return "".join([alignment.target[tp]
                    if tp >= 0 else "N" 
                    for tp in target_positions])

def find_best_alignment(aligner, sequence, target_seq):
    seq = Seq.Seq(sequence)
    alignments = []
    strands = []
    for strand in ["+", "-"]:
        if strand == "-":
            _seq = seq.reverse_complement()
        else:
            _seq = seq
        _alignments = aligner.align(_seq,target_seq)
        alignments.extend(_alignments)
        strands.extend([strand]*len(_alignments))
    best_alignment = max(alignments, key=lambda x: x.score)
    best_strand = strands[alignments.index(best_alignment)]
    return best_strand, best_alignment

def hamming(s1,s2):
    assert len(s1) == len(s2)
    return sum([a!=b for a,b in zip(s1,s2)])/len(s1)

def find_seq_with_minimum_distance(sequence, barcode_sequences):
    assert len(sequence) == len(barcode_sequences[0])
    distances = [hamming(sequence,b) for b in barcode_sequences]
    min_ind = np.argmin(distances)
    return min_ind, distances[min_ind]

def alignment_indices_to_pos(indices):
    query_indices = indices[0].astype(np.float32).copy()
    reference_linker_indices = indices[1].astype(np.float32).copy()
    query_indices[query_indices<0] = np.nan
    reference_linker_indices[reference_linker_indices<0] = np.nan
    return np.nanmean(query_indices-reference_linker_indices)