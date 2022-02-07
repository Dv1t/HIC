import pandas as pd

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
               '20', '21', '22', "X"]
min_bin_dist = 60000

for chrom_number in chromosomes:

    sv_master_table = pd.read_csv("result_prost.csv", delimiter="\t")
    sv_master_table = sv_master_table[(sv_master_table.chrom1 == sv_master_table.chrom2) &
                                      (sv_master_table.chrom1 == chrom_number)]
    sv_master_table = sv_master_table.loc[(sv_master_table["start1"]-sv_master_table["start2"]).abs() > min_bin_dist]
    left_bins = sv_master_table[["start1", "chromo_label1", "unique_id"]]
    right_bins = sv_master_table[["start2", "chromo_label2", "unique_id"]]

    right_bins = right_bins.rename(columns={"start2": "start", "chromo_label2": "chromo_label"})
    left_bins = left_bins.rename(columns={"start1": "start", "chromo_label1": "chromo_label"})


    all_bins = pd.concat([left_bins, right_bins], axis=0)
    all_bins["sv_id"] = all_bins.index #index from left_bins and right_bins dataframes
    all_bins = all_bins.reset_index()

    combinations_2d = [[start1, start2, chromo_label1, chromo_label2, unique_id1 == unique_id2, chrom_number, unique_id1, sv_id1, sv_id2]
                       for (sv_id1, start1, chromo_label1, unique_id1)
                            in (zip(all_bins['sv_id'], all_bins['start'], all_bins['chromo_label'], all_bins['unique_id']))
                       for (sv_id2, start2, chromo_label2, unique_id2)
                            in (zip(all_bins['sv_id'], all_bins['start'], all_bins['chromo_label'], all_bins['unique_id']))
                       if sv_id1 > sv_id2]

    mixed_bins = pd.DataFrame(combinations_2d,
                              columns=["start1", "start2", "chromo_label1", "chromo_label2", "same_patient", "chrom",
                                       "unique_id", 'sv_id1', 'sv_id2'])

    mixed_bins.to_csv(f"prost_sv_combinations/prost_sv_combinations{chrom_number}.csv", sep="\t")