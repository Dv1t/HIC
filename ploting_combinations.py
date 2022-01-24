import pandas as pd
import combination_classification
from cooler_extended import CoolerExtended
import seaborn as sns
import numpy as np

resolution = 40
tissue_type = "prost"
filepath = f"data/HiC40к80k200к400к800к.mcool::/resolutions/{resolution*1000}"
c = CoolerExtended(filepath)
bases_in_bin = c.binsize
sv_master_table = pd.read_csv(f"result_{tissue_type}.csv", delimiter="\t")

chomosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
all_hic_df = pd.DataFrame(columns=['score', 'type'])


def getShiftedBins(sv):
    sv_bins_x = sv["start1"]
    sv_bins_y = sv["start2"]
    sv_shift_bins_x = (sv_bins_x//bases_in_bin+10 % matrix.shape[0])*bases_in_bin
    sv_shift_bins_y = (sv_bins_y//bases_in_bin+10 % matrix.shape[0])*bases_in_bin
    sv_shift_bins = pd.DataFrame({"start1": sv_shift_bins_x, "start2": sv_shift_bins_y})
    return sv_shift_bins


for cr_number in chomosomes:
    chr_number = cr_number.__str__()
    # getting normilized hiс matrix by chromosome
    matrix = c.hic_matrices_normalized["chr" + chr_number]
    # getting SVs
    sv_for_current_chr = sv_master_table[sv_master_table["chrom1"] == chr_number][
        sv_master_table["chrom2"] == chr_number]

    sv_for_chr = sv_for_current_chr[sv_for_current_chr["chromo_label1"] != "High confidence"]
    chr_for_hic = sv_for_current_chr[sv_for_current_chr["chromo_label1"] == "High confidence"]

    sv_for_chr = sv_for_chr.drop(columns=['unique_id', 'sv_file_id', 'chrom1', 'chrom2',
                                          'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass',
                                          'svmethod', 'chromo_label1', 'chromo_label2'])
    chr_for_hic = chr_for_hic.drop(columns=['unique_id', 'sv_file_id', 'chrom1', 'chrom2',
                                            'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass',
                                            'svmethod', 'chromo_label1', 'chromo_label2'])
    sv_for_chr['combination_class'] = 'sv'
    chr_for_hic['combination_class'] = 'chromo'

    sv_shifted = getShiftedBins(sv_for_chr)
    sv_shifted['combination_class'] = 'sv_shift_5mb'

    chr_shifted = getShiftedBins(chr_for_hic)
    chr_shifted['combination_class'] = 'chromo_shift_5mb'

    mixed_bins = combination_classification.get_combinations_data_frame(
        f"{tissue_type}_sv_combinations/{tissue_type}_sv_combinations", cr_number)

    bins_array = sv_for_chr.append(chr_for_hic)
    bins_array = bins_array.append(sv_shifted)
    bins_array = bins_array.append(chr_shifted)
    bins_array = bins_array.append(mixed_bins)

    bins_array = bins_array.groupby(['combination_class'], as_index=False).agg({'start1': list, 'start2': list})

    hic_array = [[c.get_hic_score(row, chr_number), row['combination_class']] for index, row
                 in bins_array.iterrows()]

    hic_df = pd.DataFrame(hic_array, columns=['score', 'type'])
    if all_hic_df.empty:
        all_hic_df = all_hic_df.append(hic_df)
    else:
        all_hic_df['score'].append(hic_df['score'])

sns.set_theme(style="whitegrid")
sns.set(rc={'figure.figsize': (15, 15)})

plot2 = sns.boxplot(data=all_hic_df.explode('score'), x='type', y='score')
plot2.set_title("All chromosomes")
# plot2.set_yscale("log")
plot2.set_ylim([0, 3])
plot2.axhline(1)
plot2.figure.savefig(f"{tissue_type}_boxplot_combinations{resolution}k.pdf")

for sv_type, sv_score in zip(all_hic_df["type"],all_hic_df["score"]):
    print(resolution, sv_type, np.nanmean(sv_score))
print('............................................')
