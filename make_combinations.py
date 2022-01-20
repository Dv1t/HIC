import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

# %%
chomosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
for chr in chomosomes:

    chrom_number = chr.__str__()
    sv_master_table = pd.read_csv("result2.csv", delimiter="\t")
    sv_master_table = sv_master_table[(sv_master_table.chrom1 == sv_master_table.chrom2) &
                                      (sv_master_table.chrom1 == chrom_number)]

    left_bins = sv_master_table[["start1", "chromo_label1", "unique_id"]]
    right_bins = sv_master_table[["start2", "chromo_label2", "unique_id"]]

    right_bins = right_bins.rename(columns={"start2": "start", "chromo_label2": "chromo_label"})
    left_bins = left_bins.rename(columns={"start1": "start", "chromo_label1": "chromo_label"})


    all_bins = pd.concat([left_bins, right_bins], axis=0)

    combinations_2d = [[start1, start2, chromo_label1, chromo_label2, unique_id1 == unique_id2, chrom_number]
                       for sv_id1, (start1, chromo_label1, unique_id1)
                            in enumerate(zip(all_bins['start'], all_bins['chromo_label'], all_bins['unique_id']))
                       for sv_id2, (start2, chromo_label2, unique_id2)
                            in enumerate(zip(all_bins['start'], all_bins['chromo_label'], all_bins['unique_id']))
                       if sv_id1 > sv_id2]

    mixed_bins = pd.DataFrame(combinations_2d,
                              columns=["start1", "start2", "chromo_label1", "chromo_label2", "same_patient", "chrom"])

    mixed_bins.to_csv("prost_sv_combinations" + chrom_number + ".csv", sep="\t")