import pandas as pd
import combination_classification
from cooler_extended import CoolerExtended
import seaborn as sns
import numpy as np
import math

resolution = 40
tissue_type = "prost"
filepath = f"data/HiC40к80k200к400к800к.mcool::/resolutions/{resolution*1000}"
c = CoolerExtended(filepath)
bases_in_bin = c.binsize
sv_master_table = pd.read_csv(f"result_{tissue_type}.csv", delimiter="\t")

chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
all_hic_df = pd.DataFrame(columns=['score', 'type'])


def getShiftedBins(sv):
    sv_bins_x = sv["start1"]
    sv_bins_y = sv["start2"]
    sv_shift_bins_x = ((sv_bins_x//bases_in_bin+20) % matrix.shape[0])*bases_in_bin
    sv_shift_bins_y = ((sv_bins_y//bases_in_bin+20) % matrix.shape[0])*bases_in_bin
    sv_shift_bins = pd.DataFrame({"start1": sv_shift_bins_x, "start2": sv_shift_bins_y})
    return sv_shift_bins


for cr_number in chromosomes:
    chr_number = cr_number.__str__()
    # getting normilized hiс matrix by chromosome
    matrix = c.hic_matrices_normalized["chr" + chr_number]
    # getting SVs
    sv_for_chr = pd.DataFrame()
    sv_for_chr = sv_master_table[(sv_master_table["chrom1"] == chr_number)&
                                 (sv_master_table["chrom2"] == chr_number)]

    sv_for_chr = sv_for_chr.drop(columns=['unique_id', 'sv_file_id', 'chrom1', 'chrom2',
                                          'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass',
                                          'svmethod', 'chromo_label1', 'chromo_label2'])
    sv_for_chr['combination_class'] = 'sv'

    sv_shifted = getShiftedBins(sv_for_chr)
    sv_shifted['combination_class'] = 'sv_shift_5mb'

    mixed_bins = combination_classification.get_combinations_data_frame(
        f"{tissue_type}_sv_combinations/{tissue_type}_sv_combinations", cr_number, False)

    random_bins = pd.DataFrame(np.random.randint(0, matrix.shape[0]*resolution*1000, size=(10000, 2)), columns=['start1', 'start2'])
    random_bins['combination_class'] = 'random'

    bins_array = sv_for_chr.append(mixed_bins)

    bins_array = bins_array.groupby(['combination_class'], as_index=False).agg({'start1': list, 'start2': list})

    hic_array = [[c.get_hic_score(row, chr_number), row['combination_class']] for index, row
                 in bins_array.iterrows()]
    hic_df = pd.DataFrame(hic_array, columns=['score', 'type'])
    all_hic_df = all_hic_df.append(hic_df)
all_hic_df = all_hic_df.explode('score')
#%%
sns.set_theme(style="whitegrid")
sns.set(rc={'figure.figsize': (20, 20)})
plot2 = sns.boxplot(data=all_hic_df.explode('score'), x='type', y='score',
                    order=["sv", "comb_s", "comb_d"])
plot2.set_title("All chromosomes", fontsize=30)
# plot2.set_yscale("log")
plot2.set_xticklabels(plot2.get_xticklabels(), fontsize=30)
plot2.set_yticklabels(plot2.get_yticklabels(), fontsize=30)
plot2.set_ylabel(plot2.get_ylabel(), fontsize=30)
plot2.set_xlabel(plot2.get_xlabel(), fontsize=30)
plot2.set_ylim([0, 3])
plot2.axhline(1)
plot2.figure.savefig(f"{tissue_type}_boxplot_combinations{resolution}k_no_chromo.pdf")
#%%
all_hic_df = all_hic_df.groupby(['type'], as_index=False).agg({'score': list})
mean_table = pd.read_csv(f"{tissue_type}_mean_table_nosh.csv", delimiter="\t")
#mean_table = pd.DataFrame()
mean_table[f'{str(resolution)}k'] = [np.nanmean(sv_score) for sv_score in all_hic_df["score"]]
mean_table['combination_class'] = [sv_type for sv_type in all_hic_df["type"]]
mean_table.to_csv(f"{tissue_type}_mean_table_nosh.csv", sep="\t")

#median_table = pd.DataFrame()
median_table = pd.read_csv(f"{tissue_type}_median_table_nosh.csv", delimiter="\t")
def median(array):
    a = [0 if math.isnan(x) else x for x in array]
    return np.percentile(a, 50)
median_table[f'{str(resolution)}k'] = [median(sv_score) for sv_score in all_hic_df["score"]]
median_table['combination_class'] = [sv_type for sv_type in all_hic_df["type"]]
median_table.to_csv(f"{tissue_type}_median_table_nosh.csv", sep="\t")

