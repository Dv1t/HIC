import pandas as pd

def __get_combination_class(label1, label2, same_patient):
    if (label1 == "High confidence") & (label2 == "High confidence"):
        if same_patient:
            return "chromo_chromo_s"
        else:
            return "chromo_chromo_d"
    if (label1 == "High confidence") | (label2 == "High confidence"):
        if same_patient:
            return "sv_chromo_s"
        else:
            return "sv_chromo_d"
    if same_patient:
        return "sv_sv_s"
    else:
        return "sv_sv_d"

def get_combinations_data_frame(sv_path, cr_number):
    mixed_bins = pd.read_csv(f"{sv_path}{cr_number.__str__()}.csv", delimiter="\t")
    mixed_bins['combination_class'] = [__get_combination_class(l1, l2, s) for l1, l2, s in zip(mixed_bins.chromo_label1, mixed_bins.chromo_label2, mixed_bins.same_patient)]
    return mixed_bins
