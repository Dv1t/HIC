from hic2cool import hic2cool_convert
chomosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]

#if __name__ == '__main__':
 #   for chr in chomosomes:
 #       hic2cool_convert("hic_files/HPNE.chr"+chr.__str__()+".hic", "cool_files_pancreatic_50k/cool"+chr.__str__()+".cool", 100000)

#%%
import cooler
cooler.merge_coolers(output_uri="cool_files_pancreatic_50k/cool_merged.cool",
                        input_uris=['cool_files_pancreatic_50k/cool{}.cool'.format(i) for i in chomosomes],
                        mergebuf = int(20e6))

#%%
from  cooler_extended import  CoolerExtended
import matplotlib.pyplot as plt
import numpy as np
import cooler


c = cooler.Cooler("cool_files_pancreatic_100k/cool_merged.cool")
arr2 = c.matrix(balance=False, sparse=True)[1000:1200, 1000:1200].toarray()

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(np.log10(arr2), cmap='YlOrRd')
fig.colorbar(im)
plt.savefig("temp.pdf")
#%%
import cooler
c = cooler.Cooler("cool_files_pancreatic_50k/cool_merged.cool")
print(c.binsize)
