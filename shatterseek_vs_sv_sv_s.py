#%%
import pandas as pd
import combination_classification
from cooler_extended import CoolerExtended


#%%
resolution = 800
tissue_type = "prost"
filepath = f"data/HiC40к80k200к400к800к.mcool::/resolutions/{resolution*1000}"
c = CoolerExtended(filepath)
bases_in_bin = c.binsize
sv_master_table = pd.read_csv(f"result_{tissue_type}.csv", delimiter="\t")

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
               '20', '21', '22', "X"]
patients = sv_master_table.unique_id
patients.drop_duplicates(inplace=True)
main_df = pd.DataFrame(columns=chromosomes, index=patients)
#%%
all_combinations = []
for chromosome in chromosomes:
    combinations = combination_classification.\
            get_combinations_same_patient_data_frame(f"{tissue_type}_sv_combinations/{tissue_type}_sv_combinations", chromosome)
    combinations[chromosome] = [c.get_single_hic_score(start1, start2, chromosome) for start1, start2 in zip(combinations['start1'],
                                                                                            combinations['start2'])]
    combinations.drop(columns=['start1', 'start2', 'Unnamed: 0'], inplace=True)
    combinations = combinations.groupby('unique_id').mean()
    main_df.update(combinations)
#%%
high_chromo = sv_master_table[(sv_master_table["chromo_label1"] == "High confidence")
                              & (sv_master_table["chromo_label2"] == "High confidence")
                              & (sv_master_table["chrom1"] == sv_master_table["chrom2"])]
high_chromo.drop(columns=['sv_file_id', 'start1', 'start2', 'Unnamed: 0', 'chrom2', 'sv_id', 'pe_support', 'strand1',
                          'strand2', 'svclass', 'svmethod', 'chromo_label1', 'chromo_label2'], inplace=True)
high_chromo = high_chromo.groupby('unique_id').agg({'chrom1':lambda x: ', '.join(list(dict.fromkeys(x)))})
#high_chromo = [chrs.drop_duplicates for chrs in high_chromo]
high_chromo = high_chromo.rename(columns={'chrom1' : "high_chromo"})
main_df["high_chromo"] = high_chromo
main_df.update(high_chromo)
main_df.to_csv("shatterseek_vs_sv_sv_s.csv", sep="\t")

