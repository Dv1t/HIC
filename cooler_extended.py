import cooler
import numpy as np
import math

import pandas
from scipy import stats
import savitzky_golay

def distribution_at_dist(arr, d):
    n = arr.shape[0]
    return np.array([arr[i, j] for i, j in zip(range(0, n - d), range(d, n))])


def normalize_intra(arr):
    n = arr.shape[0]
    averages_at_dist = [np.nanmean(distribution_at_dist(arr, d)) for d in range(0, n)]
    ans = np.zeros_like(arr, dtype='float64')
    for i in range(n):
        for j in range(n):
            ans[i, j] = arr[i, j] / averages_at_dist[abs(i - j)]
    return ans


class CoolerExtended(cooler.Cooler):
    __threshold = 0.8

    def __init__(self, filepath):
        super().__init__(filepath)
        self.hic_matrices_normalized = {}
        self.bases_in_bin = self.binsize
        for current_chr in self.chromnames:
            mat = self.matrix(balance=False).fetch(current_chr)
            if "chr" not in current_chr:
                current_chr = "chr" + current_chr
            mat_nan = self.__zeros_to_nan(mat)
            mat_norm = normalize_intra(mat_nan)
            self.hic_matrices_normalized[current_chr] = mat_norm

    def __zeros_to_nan(self, arr):
        arr = arr.astype(float)
        n = arr.shape[0]
        for i in range(len(arr)):
            if ((arr[i] == 0).sum(0) / n) >= self.__threshold:
                arr[i] = np.nan
                arr[:, i] = np.nan
        return arr

    def get_hic_score(self, table, chr_number, min_bin_dist=60000, max_bin_dist=math.inf):
        if "chr" not in chr_number:
            chr_number = "chr" + chr_number
        bins_x = np.array(table["start1"])
        bins_y = np.array(table["start2"])

        bins_ok = [(abs(bin_x - bin_y) > min_bin_dist) for bin_x, bin_y
                   in zip(bins_x, bins_y)]
        bins_x = bins_x[bins_ok] // self.bases_in_bin
        bins_y = bins_y[bins_ok] // self.bases_in_bin

        bins_ok = [(abs(bin_x - bin_y) > 1) for bin_x, bin_y
                   in zip(bins_x, bins_y)]
        bins_x = bins_x[bins_ok]
        bins_y = bins_y[bins_ok]

        hic_score = []
        for x, y in zip(bins_x, bins_y):
            try:
                if not math.isnan(self.hic_matrices_normalized[chr_number][x][y]):
                    hic_score.append(float(self.hic_matrices_normalized[chr_number][x][y]))
            except IndexError:
                pass
        return hic_score

    def get_single_hic_score(self, start1, start2, chr_number, min_bin_dist=60000, max_bin_dist=math.inf):
        if "chr" not in str(chr_number):
            chr_number = "chr" + str(chr_number)
        if abs(start1 - start2) <= min_bin_dist:
            return np.nan
        start1 = (start1 // self.bases_in_bin).astype(int)
        start2 = (start2 // self.bases_in_bin).astype(int)
        if abs(start1 - start2) <= 1:
            return np.nan
        try:
            return self.hic_matrices_normalized[chr_number][start1][start2]
        except IndexError:
            return ''

    def calculate_quality(self, chr_number):
        quality_distribution = []
        for dist in range(self.hic_matrices_normalized[chr_number].shape[0]):
            quality_distribution.append(float(stats.skew(np.diagonal(self.hic_matrices_normalized[chr_number], dist), nan_policy='omit')))
        quality_distribution = savitzky_golay.savitzky_golay(np.asarray(quality_distribution), 31, 4)
        return quality_distribution

    def calculate_mean_minus_median(self, chr_number):
        result = []
        for dist in range(self.hic_matrices_normalized[chr_number].shape[0]):
            diagonal = np.diagonal(self.hic_matrices_normalized[chr_number], dist)
            diagonal = diagonal[~np.isnan(diagonal)]
            if len(diagonal) == 0:
                continue
            result.append(np.nanmean(diagonal) - np.percentile(diagonal, 50))
        result = savitzky_golay.savitzky_golay(np.asarray(result), 31, 4)
        return result


