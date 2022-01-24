import pandas as pd


def __get_combination_class(label1, label2, same_patient):
    result = "sv_sv"
    if (label1 == "High confidence") | (label2 == "High confidence"):
        result = "sv_chromo"
    if (label1 == "High confidence") & (label2 == "High confidence"):
        result = "chromo_chromo"
    if same_patient:
        return f"{result}_s"
    else:
        return f"{result}_d"


def get_combinations_data_frame(sv_path, cr_number):
    mixed_bins = pd.read_csv(f"{sv_path}{cr_number.__str__()}.csv", delimiter="\t")
    mixed_bins['combination_class'] = [__get_combination_class(l1, l2, s) for l1, l2, s in
                                       zip(mixed_bins.chromo_label1, mixed_bins.chromo_label2, mixed_bins.same_patient)]
    mixed_bins = mixed_bins.drop(columns=['chromo_label1', 'chromo_label2', 'same_patient', 'chrom'])
    return mixed_bins
