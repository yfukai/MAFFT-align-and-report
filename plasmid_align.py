# %%
import os
from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np
import itertools
from Bio import Align
from Bio.Align import substitution_matrices
from aligner import DNAPairwiseAligner
from utils import extract_good_region, summarize_mismatch



# %%
rev_complement = False
ab1_file = "testdata/230825-mixed5#26-SPYF009-PREMIX.ab1"
gb_file = "testdata/pyf042-601x6-ligated-backbone-a2-bsmbi-bbsi.gb"
export_directory = "testdata/output/"
feature_name = "spacer"
feature_type = "misc_feature"

os.makedirs(export_directory, exist_ok=True)

def extract_and_align(ab1_file,gb_file,feature_name=None,feature_type="misc_feature"):
    record = SeqIO.read(ab1_file, "abi")
    quality = record.letter_annotations["phred_quality"]
    good_region = extract_good_region(quality)
    seq_good = record.seq[good_region[0]:good_region[1]]

    ref_record = SeqIO.read(gb_file, "genbank")

    aligner = DNAPairwiseAligner()
    align_res = aligner.align(ref_record.seq,seq_good)
    best_alignment = max(align_res,key=lambda x: x.score)

    target_features = [feature for feature in ref_record.features 
                       if feature_name in feature.qualifiers["label"][0] 
                       and feature.type == feature_type]

    feature_sequences = []
    for feature in target_features:
        start = feature.location.start
        end = feature.location.end
        indices_starts = np.nonzero(best_alignment.indices[0] == start)[0]
        indices_ends = np.nonzero(best_alignment.indices[1] == end)[0]
        if len(indices_starts) == 0 and len(indices_ends) == 0:
            insertion_count = None
            deletion_count = None
            mismatch_count = None
            query_seq_wo_insertion = None
            query_seq_w_insertion = None
        else:
            if len(indices_starts) == 0:
              indices_starts = [0]
            if len(indices_ends) == 0:
                indices_ends = [len(best_alignment.indices[0])]

            assert len(indices_starts) == 1
            indices_start = indices_starts[0]
            assert len(indices_ends) == 1
            indices_end = indices_ends[0]

        feature_sequences.append({
            "feature_name":feature.qualifiers["label"][0],
            "feature_start":start,
            "feature_end":end,
            "insertion":insertion_count,
            "deletion":deletion_count,
            "mismatch": mismatch_count,
            "query_seq_wo_insertion":query_seq_wo_insertion,
            "query_seq_w_insertion":query_seq_w_insertion,
        })

    return quality, good_region, best_alignment, feature_sequences

# %%
plt.plot(record.letter_annotations["phred_quality"])
plt.axvline(good_region[0],color="k")
plt.axvline(good_region[1],color="k")
plt.xlabel("Position (bp)")
plt.ylabel("Phred quality")
plt.savefig(export_directory)

# %%
len(best_align[0])

# %%
len(best_align.indices[0])

# %%
# %%
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
rev_complement = False
ab1_file = "testdata/example2/241023-pYF061-#13-SPYF009-PREMIX_ID1.ab1"
ab1_file = "testdata/example2/241023-pYF061-#2-SPYF009-PREMIX_ID1.ab1"
gb_file = "testdata/example2/pyf061-601x6-30bp-linker-bsmbi.gb"

export_directory = "testdata/output/"
feature_names = ["linker","601"]
feature_type = "misc_feature"
record = SeqIO.read(ab1_file, "abi")
quality = record.letter_annotations["phred_quality"]
good_region = extract_good_region(quality)
seq_good = record.seq[good_region[0]:good_region[1]]

ref_record = SeqIO.read(gb_file, "genbank")

aligner = DNAPairwiseAligner()
align_res = aligner.align(ref_record.seq,seq_good)
best_alignment = max(align_res,key=lambda x: x.score)

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
        indices_start = None
        indices_end = None
        insertion_count = None
        deletion_count = None
        mismatch_count = None
        query_seq_wo_insertion = None
        query_seq_w_insertion = None
        target_seq = None
        query_seq = None
        fraction_aligned = None
    else:
        if len(indices_starts) == 0:
          indices_starts = [0]
        if len(indices_ends) == 0:
            indices_ends = [len(best_alignment.indices[0])]
        assert len(indices_starts) == 1
        indices_start = int(indices_starts[0])
        assert len(indices_ends) == 1
        indices_end = int(indices_ends[0])
        insertion_count, deletion_count, mismatch_count = summarize_mismatch(
            best_alignment,indices_start,indices_end
        )
        target_seq = best_alignment[0][indices_start:indices_end]
        query_seq = best_alignment[1][indices_start:indices_end]
        query_seq_w_insertion = "".join([q.upper() if t != "-" else q.lower() for t,q in zip(target_seq,query_seq)])
        query_seq_wo_insertion = "".join([q for q in query_seq_w_insertion if q=="-" or q.isupper()])
        fraction_aligned = len(query_seq_wo_insertion)/(end-start)
    
    feature_sequences.append({
        "feature_name":feature.qualifiers["label"][0],
        "feature_start":start,
        "feature_end":end,
        "fraction_aligned":fraction_aligned,
        "indices_start":indices_start,
        "indices_end":indices_end,
        "insertion":insertion_count,
        "deletion":deletion_count,
        "mismatch": mismatch_count,
        "target_seq":target_seq,
        "query_seq":query_seq,
        "query_seq_wo_insertion":query_seq_wo_insertion,
        "query_seq_w_insertion":query_seq_w_insertion,
    })
    #export pairwisealigner result to SAM file
#    with open(export_directory+feature.qualifiers["label"][0]+".sam","w") as f:
#        f.write(best_alignment.to_sam())
pd.DataFrame(feature_sequences)
# %%

# %%
