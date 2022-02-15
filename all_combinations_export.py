import numpy as np
import pandas as pd
from cooler_extended import CoolerExtended
import combination_classification

resolution = 800
tissue_type = "prost"
filepath = f"data/HiC40к80k200к400к800к.mcool::/resolutions/{resolution*1000}"
c = CoolerExtended(filepath)
bases_in_bin = c.binsize
sv_master_table = pd.read_csv(f"result_{tissue_type}.csv", delimiter="\t")

chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
all_hic_df = pd.DataFrame(columns=['type', f'score_{resolution}k', 'chr', 'start1', 'start2'])
for cr_number in chromosomes:
    chr_number = cr_number.__str__()

    sv_for_current_chr = sv_master_table[(sv_master_table["chrom1"] == chr_number)&
                                 (sv_master_table["chrom2"] == chr_number)]

    sv_for_chr = sv_for_current_chr
    chr_for_hic = sv_for_current_chr[sv_for_current_chr["chromo_label1"] == "High confidence"]

    sv_for_chr = sv_for_chr.drop(columns=['unique_id', 'sv_file_id', 'chrom1', 'chrom2',
                                          'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass',
                                          'svmethod', 'chromo_label1', 'chromo_label2'])
    chr_for_hic = chr_for_hic.drop(columns=['unique_id', 'sv_file_id', 'chrom1', 'chrom2',
                                            'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass',
                                            'svmethod', 'chromo_label1', 'chromo_label2'])
    sv_for_chr['combination_class'] = 'sv'
    chr_for_hic['combination_class'] = 'chromo'

    mixed_bins = combination_classification.get_combinations_data_frame(
        f"{tissue_type}_sv_combinations/{tissue_type}_sv_combinations", cr_number)

    mixed_bins_no_shatter = combination_classification.get_combinations_data_frame(
        f"{tissue_type}_sv_combinations/{tissue_type}_sv_combinations", cr_number, False)

    bins_array = sv_for_chr.append(chr_for_hic)
    bins_array = bins_array.append(mixed_bins)
    bins_array = bins_array.append(mixed_bins_no_shatter)

    bins_array = bins_array.groupby(['combination_class'], as_index=False).agg({'start1': list, 'start2': list})

    hic_array = [[c.get_hic_score(row, chr_number), row['combination_class'], row['start1'], row['start2']]
                 for index, row
                 in bins_array.iterrows()]
    hic_df = pd.DataFrame(hic_array, columns=[f'score_{resolution}k', 'type', 'start1', 'start2'])
    hic_df['chr'] = cr_number
    hic_df = hic_df[(hic_df.type=='sv')|(hic_df.type=='comb_s')|(hic_df.type=='comb_d')|(hic_df.type=='chromo')
                    |(hic_df.type=='chromo_chromo_d')|(hic_df.type=='chromo_chromo_s')]
    all_hic_df = all_hic_df.append(hic_df)

all_hic_df = all_hic_df.explode(f'score_{resolution}k')
all_hic_df.to_csv('lesha_table800.csv', sep=',',  index=False)

#%%

import numpy as np

#df = pd.read_csv('lesha_table.csv', delimiter=',')
#print(df.shape[0])
#проверяем что всё сходится
for ctype, df_ctype in all_hic_df.groupby('type'):
    print(ctype, np.mean(df_ctype.score_800k))