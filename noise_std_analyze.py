from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run023/run023*.fits')
print(len(flist))
#flist = glob('/gpfs01/astro/workarea/lmezini/code/test.fits')
# read each file and combine into one big array
all_data = np.zeros((8,8))
for i in flist:
    data = eu.io.read(i)
    all_data +=data
avg = all_data/len(flist)
print("avg standard devs: ",avg)
plt.imshow(avg,interpolation='nearest', cmap='gray',vmin = np.min(avg),vmax =np.max(avg)) 
plt.colorbar()
plt.title("Avg Noise STD Bgrms=10")
plt.savefig("avg_noise_std_neigh_bgrms10_2.png")
