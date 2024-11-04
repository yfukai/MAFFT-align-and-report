import numpy as np
import itertools
from Bio import Align
from Bio.Align import substitution_matrices


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

class DNAPairwiseAligner(Align.PairwiseAligner):
    def __init__(self, scoring=None, match_score=1,shared_score=0.8,mismatch_score=-1, 
                 open_gap_score=-2, extend_gap_score=-4, 
                 query_open_gap_score=-10, 
                 query_extend_gap_score=-10, 
                 mode="local",
                 **kwargs):
        letters = "".join(nuc_codes.keys())

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

        custom_substitution_matrix = substitution_matrices.Array(
            letters, 2, score_matrix
        )
        super().__init__(scoring, **kwargs)
        self.open_gap_score = open_gap_score
        self.extend_gap_score = extend_gap_score
        self.substitution_matrix = custom_substitution_matrix
        self.query_left_open_gap_score = query_open_gap_score
        self.query_left_extend_gap_score = query_extend_gap_score
        self.query_right_open_gap_score = query_open_gap_score
        self.query_right_extend_gap_score = query_extend_gap_score
        self.mode = mode
