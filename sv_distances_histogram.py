import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import seaborn as sns

sv_table = pd.read_csv(f"result_prost.csv", delimiter="\t")

sv_table['distance'] = [(abs(start1-start2)) for start1, start2 in zip(sv_table.start1, sv_table.start2)]
#%%
# An "interface" to matplotlib.axes.Axes.hist() method
plot = sns.histplot(x=sv_table.distance, bins='auto', color='#0504aa',
                            alpha=0.7, log_scale=[10, False])
plt.grid(axis='y', alpha=0.75)
plt.xlabel('SV length in bases')
plt.ylabel('Frequency')
#plt.xticks(bins[::2], np.arange(0, max(bins), 10000), rotation=90)
#plt.subplots_adjust(bottom=0.3, top=0.99, right=0.99, left=0.1)
plt.title('SV distances histogram')
plt.axvline(x=60000, color="red")
maxfreq = 1400
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.savefig('prost_sv_dist_histogram_log.pdf')
plt.show()
