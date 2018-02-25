from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run026/run026*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/code/test*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/test.fits')
# read each file and combine into one big array
neigh_data = np.zeros((8,8))
cen_data = np.zeros((8,8))
data = eu.io.read(flist)
for i in range(len(flist)):
    print(data['neigh_pix_noise_std'][i])
    neigh_data +=data['neigh_pix_noise_std'][i]
    cen_data +=data['cen_pix_noise_std'][i]

neigh_avg = neigh_data/len(flist)
cen_avg= cen_data/len(flist)
print("neigh avg standard devs: ",neigh_avg)
print("cen avg standard devs: ",cen_avg)
plt.imshow(neigh_avg,interpolation='nearest', cmap='gray',vmin = np.min(avg),vmax =np.max(avg)) 
plt.colorbar()
plt.title("Avg Noise STD Bgrms=10")
plt.savefig("test.png")
#fitsio.write("avg_noise_std_neigh_bgrms10.fits", avg, clobber=True) 
