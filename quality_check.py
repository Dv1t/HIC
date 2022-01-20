from  cooler_extended import  CoolerExtended
import matplotlib.pyplot as plt
import numpy as np

filepath = "data/panc100k200k400k500k800k.mcool::/resolutions/100000"
c = CoolerExtended(filepath)

#%%
chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]

chr_skew = []
for current_chr in chromosomes:
    chr_skew.append(c.calculate_quality("chr"+current_chr.__str__()))
    plt.plot(chr_skew[len(chr_skew)-1], linewidth=0.7)
plt.grid(True)
plt.suptitle("Chromosomes skewness by distance to diagonal")
plt.ylabel("Skewness")
plt.xlabel("Distance")
plt.yticks(np.arange(-2, 31, 2.0))
plt.xticks(np.arange(0, len(chr_skew[0])+1, 400.0))
plt.savefig("quality_check_skew_100k_panc.pdf")
plt.show()


#%%

chr_mean_all = []
for current_chr in chromosomes:
    chr_mean = c.calculate_mean_minus_median(f'chr{current_chr}')
    chr_mean_all.append(chr_mean)
    plt.plot(chr_mean)
plt.grid(True)
plt.suptitle(f"Chromosomes mean subtract median by distance to diagonal")
plt.ylabel("Mean subtract median")
plt.xlabel("Distance")
plt.yticks(np.arange(-0.15,1.15, 0.1))
plt.xticks(np.arange(0, len(chr_mean_all[0]) + 1, 20.0))
plt.savefig(f"графики/среднее_минус_медиана_простата_800k/mean_median_all_chr.pdf")
plt.show()

