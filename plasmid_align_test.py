# %%
import importlib
import align_sequence_files
from utils import summarize_mismatch
importlib.reload(align_sequence_files)
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
# %%
ab1_file = "testdata/example2/241023-pYF061-#13-SPYF009-PREMIX_ID1.ab1"
ab1_file = "testdata/example2/241023-pYF061-#2-SPYF009-PREMIX_ID1.ab1"
gb_file = "testdata/example2/pyf061-601x6-30bp-linker-bsmbi.gb"
rev_complement = False
feature_names = ["linker","601"]
feature_type = "misc_feature"
best_alignment, features_df = align_sequence_files.align_sequence_files(ab1_file, gb_file, rev_complement, feature_names, feature_type)
features_df



# %%
ab1_file = "testdata/example2/241023-pYF061-#2-SPYF010-PREMIX_ID4.ab1"
gb_file = "testdata/example2/pyf061-601x6-30bp-linker-bsmbi.gb"
rev_complement = True
feature_names = ["linker","601"]
feature_type = "misc_feature"
best_alignment, features_df = align_sequence_files.align_sequence_files(
    ab1_file, gb_file, rev_complement, feature_names, feature_type,quality_threshold=25)
insertion_count, deletion_count, mismatch_count = summarize_mismatch(
    best_alignment
)
display(features_df)
print(insertion_count, deletion_count, mismatch_count)

# %%
record = SeqIO.read(ab1_file, "abi")
# %%
record.annotations["abif_raw"].keys()

# %%
from collections import defaultdict

channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
trace = defaultdict(list)
for c in channels:
    trace[c] = record.annotations["abif_raw"][c]
# %%
plt.plot(trace["DATA9"][:1000], color="blue")
plt.plot(trace["DATA10"][:1000], color="red")
plt.plot(trace["DATA11"][:1000], color="green")
plt.plot(trace["DATA12"][:1000], color="yellow")
plt.show()
# %%
