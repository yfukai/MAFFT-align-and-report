import numpy as np
from aligner import nuc_codes

def extract_good_region(phred_quality,threshold=40):
    phred_quality = np.array(phred_quality)
    good_region = np.where(phred_quality>threshold)[0]
    # get connected regions
    good_regions = []
    start = good_region[0]
    for i in range(1,len(good_region)):
        if good_region[i] != good_region[i-1]+1:
            good_regions.append((start,good_region[i-1]))
            start = good_region[i]
    good_regions.append((start,good_region[-1]))
    return max(good_regions,key=lambda x: x[1]-x[0])
assert extract_good_region([10,50,60,70,80,90,10,20,30,40,50,60,70,80,90,100,90,80,70,60,50,40,30,20,10]) == (10,20)

def summarize_mismatch(alignment,start=None, end=None):
    total_deletion = 0
    total_insertion = 0
    total_mismatch = 0
    for i, (t, q) in enumerate(zip(alignment[0][start:end], alignment[1][start:end])):
        if t == "-":
            total_insertion += 1
        elif q == "-":
            total_deletion += 1
        elif q not in nuc_codes[t]:
            total_mismatch += 1
    return total_insertion, total_deletion, total_mismatch

    
    