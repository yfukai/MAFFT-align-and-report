import os
from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np
import itertools
from Bio import Align
import pandas as pd
from Bio.Align import substitution_matrices
from aligner import DNAPairwiseAligner
from utils import extract_good_region, summarize_mismatch

def align_sequence_files(ab1_file, gb_file, rev_complement=False, feature_names=None, feature_type="misc_feature",quality_threshold=40):
    record = SeqIO.read(ab1_file, "abi")
    quality = record.letter_annotations["phred_quality"]
    good_region = extract_good_region(quality, threshold=quality_threshold)
    seq_good = record.seq[good_region[0]:good_region[1]]
    if rev_complement:
        seq_good = seq_good.reverse_complement()
    
    ref_record = SeqIO.read(gb_file, "genbank")
    aligner = DNAPairwiseAligner()
    align_res = aligner.align(ref_record.seq,seq_good)
    best_alignment = max(align_res,key=lambda x: x.score)

    if feature_names:
        target_features = [feature for feature in ref_record.features 
                           if any(feature_name in feature.qualifiers["label"][0] for feature_name in feature_names)
                           and feature.type == feature_type]

        feature_sequences = []
        for feature in target_features:
            start = feature.location.start
            end = feature.location.end
            indices_starts = np.nonzero(best_alignment.indices[0] == start)[0]
            indices_ends = np.nonzero(best_alignment.indices[0] == end)[0]
            if len(indices_starts) == 0 and len(indices_ends) == 0:
                alignment_start_index = None
                alignment_end_index = None
                insertion_count = None
                deletion_count = None
                mismatch_count = None
                query_seq_wo_insertion = ""
                query_seq_w_insertion = ""
                target_seq = ""
                query_seq = ""
                fraction_aligned = 0
            else:
                if len(indices_starts) == 0:
                  indices_starts = [0]
                if len(indices_ends) == 0:
                    indices_ends = [len(best_alignment.indices[0])]
                assert len(indices_starts) == 1
                alignment_start_index = int(indices_starts[0])
                assert len(indices_ends) == 1
                alignment_end_index = int(indices_ends[0])
                insertion_count, deletion_count, mismatch_count = summarize_mismatch(
                    best_alignment,alignment_start_index,alignment_end_index
                )
                target_seq = best_alignment[0][alignment_start_index:alignment_end_index]
                query_seq = best_alignment[1][alignment_start_index:alignment_end_index]
                query_seq_w_insertion = "".join([q.upper() if t != "-" else q.lower() for t,q in zip(target_seq,query_seq)])
                query_seq_wo_insertion = "".join([q for q in query_seq_w_insertion if q=="-" or q.isupper()])
                fraction_aligned = len(query_seq_wo_insertion)/(end-start)

            feature_sequences.append({
                "feature_name":feature.qualifiers["label"][0],
                "feature_start":start,
                "feature_end":end,
                "fraction_aligned":fraction_aligned,
                "insertion":insertion_count,
                "deletion":deletion_count,
                "mismatch": mismatch_count,
                "alignment_start_index":alignment_start_index,
                "alignment_end_index":alignment_end_index,
                "target_seq":target_seq,
                "query_seq":query_seq,
                "aligned_query_seq_wo_insertion":query_seq_wo_insertion,
                "aligned_query_seq_w_insertion":query_seq_w_insertion,
            })
        feature_df = pd.DataFrame(feature_sequences)
        for k in ["feature_start","feature_end","insertion","deletion","mismatch","alignment_start_index","alignment_end_index"]:
            feature_df[k] = feature_df[k].astype(pd.Int64Dtype())
    else:
        feature_df = None
    return best_alignment, feature_df
