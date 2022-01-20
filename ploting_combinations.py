import math
import numpy as np
import  pandas as pd
import matplotlib.pyplot as plt

import combination_classification
from  cooler_extended import  CoolerExtended
import seaborn as sns
from scipy import stats
#%%

filepath = "data/HiC40к80k200к400к800к.mcool::/resolutions/200000"
c = CoolerExtended(filepath)
bases_in_bin = c.binsize
sv_master_table = pd.read_csv("result2.csv",delimiter = "\t")

#%%
chomosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
sv_hic_score = []
chr_hic_score = []
hic_sv_sv_d = []
hic_sv_sv_s = []
hic_sv_chr_d = []
hic_sv_chr_s = []
hic_chr_chr_d = []
hic_chr_chr_s = []
chr_shift_hic_score = []
sv_shift_hic_score = []

for cr_number in chomosomes:
    chr_number = cr_number.__str__()
    #getting normilized hiс matrix by chromosome
    matrix = c.hic_matrices_normalized["chr"+chr_number]
    #getting SVs
    sv_for_chr = sv_master_table[sv_master_table["chrom1"]==chr_number][sv_master_table["chrom2"]==chr_number][sv_master_table["chromo_label1"]!="High confidence"]
    sv_bins_x = sv_for_chr["start1"] // bases_in_bin
    sv_bins_y = sv_for_chr["start2"] // bases_in_bin
    sv_ok = (sv_bins_x - sv_bins_y).abs() > 1
    sv_bins_x = sv_bins_x[sv_ok]
    sv_bins_y = sv_bins_y[sv_ok]

    sv_hic_score.extend(c.get_hic_score(sv_for_chr, chr_number, True))

    chr_for_hic = sv_master_table[sv_master_table["chrom1"]==chr_number][sv_master_table["chrom2"]==chr_number][sv_master_table["chromo_label1"]=="High confidence"]

    if not chr_for_hic.empty:
        chr_hic_score.extend(c.get_hic_score(chr_for_hic, chr_number, True))
        chr_bins_x = []
        chr_bins_y = []
        for x,y in zip(chr_for_hic["start1"].tolist(),chr_for_hic["start2"].tolist()):
                if abs((x//bases_in_bin)-(y//bases_in_bin))>1:
                    chr_bins_x.append(x//bases_in_bin)
                    chr_bins_y.append(y//bases_in_bin)
        chr_shift_bins_x = []
        chr_shift_bins_y = []
        for x,y in zip(chr_bins_x, chr_bins_y):
            chr_shift_bins_x.append((x+10) % matrix.shape[0])
            chr_shift_bins_y.append((y+10) % matrix.shape[0])
        chr_shift_bins = pd.DataFrame(list(zip(chr_shift_bins_x, chr_shift_bins_y)))
        chr_shift_bins.columns = ["start1", "start2"]
        chr_shift_hic_score.extend(c.get_hic_score(chr_shift_bins, chr_number))

    sv_shift_bins_x = []
    sv_shift_bins_y = []
    for x,y in zip(sv_bins_x,sv_bins_y):
        sv_shift_bins_x.append((x+10)%matrix.shape[0])
        sv_shift_bins_y.append((y+10)%matrix.shape[0])
    sv_shift_bins = pd.DataFrame(list(zip(sv_shift_bins_x, sv_shift_bins_y)))
    sv_shift_bins.columns = ["start1", "start2"]
    sv_shift_hic_score.extend(c.get_hic_score(sv_shift_bins, chr_number))

    mixed_bins = combination_classification.get_combinations_data_frame("sv_combinations_csv/sv_combinations", cr_number)

    sv_sv_d = mixed_bins[mixed_bins.combination_class == "sv_sv_d"]
    sv_sv_s = mixed_bins[mixed_bins.combination_class == "sv_sv_s"]
    sv_chr_d = mixed_bins[mixed_bins.combination_class == "sv_chromo_d"]
    sv_chr_s = mixed_bins[mixed_bins.combination_class == "sv_chromo_s"]
    chr_chr_s = mixed_bins[mixed_bins.combination_class == "chromo_chromo_s"]
    chr_chr_d = mixed_bins[mixed_bins.combination_class == "chromo_chromo_d"]

    hic_sv_sv_d.extend(c.get_hic_score(sv_sv_d, chr_number, True))
    hic_sv_sv_s.extend(c.get_hic_score(sv_sv_s, chr_number, True))
    hic_sv_chr_d.extend(c.get_hic_score(sv_chr_d, chr_number, True))
    hic_sv_chr_s.extend(c.get_hic_score(sv_chr_s, chr_number, True))
    hic_chr_chr_d.extend(c.get_hic_score(chr_chr_d, chr_number, True))
    hic_chr_chr_s.extend(c.get_hic_score(chr_chr_s, chr_number, True))
#%%
sns.set_theme(style="whitegrid")
sns.set(rc={'figure.figsize':(15,15)})

type_and_score = [[val, type]
                  for vals, type in zip(
                                        [chr_hic_score, sv_hic_score, chr_shift_hic_score, sv_shift_hic_score,hic_sv_sv_d, hic_sv_sv_s, hic_sv_chr_d, hic_sv_chr_s,hic_chr_chr_d,hic_chr_chr_s],
                                        ['chromo','sv','chromo_shift_5mb','sv_shift_5mb','sv_sv_d','sv_sv_s','sv_chr_d','sv_chr_s','chr_chr_d','chr_chr_s']
                                        )
                  for val in vals]

df = pd.DataFrame(type_and_score, columns=['score', 'type'])

plot2 = sns.boxplot(data=df, x='type', y='score')
plot2.set_title("All chromosomes")
#plot2.set_yscale("log")
plot2.set_ylim([0,3])
plot2.axhline(1)
plot2.figure.savefig("panc_boxplot_combinations500k_high.pdf")

#%%
print("chr_hic_score",np.nanmean(chr_hic_score),np.percentile((chr_hic_score),50))
print("sv_hic_score",np.nanmean(sv_hic_score),np.percentile((sv_hic_score),50))
print("chr_shift_hic_score",np.nanmean(chr_shift_hic_score),np.percentile((chr_shift_hic_score),50))
print("sv_shift_hic_score",np.nanmean(sv_shift_hic_score),np.percentile((sv_shift_hic_score),50))
print("hic_sv_sv_d",np.nanmean(hic_sv_sv_d),np.percentile((hic_sv_sv_d),50))
print("hic_sv_sv_s",np.nanmean(hic_sv_sv_s),np.percentile((hic_sv_sv_s),50))
print("hic_sv_chr_d",np.nanmean(hic_sv_chr_d),np.percentile((hic_sv_chr_d),50))
print("hic_sv_chr_s",np.nanmean(hic_sv_chr_s),np.percentile((hic_sv_chr_s),50))
print("hic_chr_chr_d",np.nanmean(hic_chr_chr_d),np.percentile((hic_chr_chr_d),50))
print("hic_chr_chr_s",np.nanmean(hic_chr_chr_s),np.percentile((hic_chr_chr_s),50))
