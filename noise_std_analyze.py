from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run028/run028*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/code/test*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/test.fits')
# read each file and combine into one big array
neigh_data = np.zeros((8,8))
cen_data = np.zeros((8,8))
data = eu.io.read(flist)
for i in range(len(flist)):
    neigh_data +=data['neigh_pix_noise_std'][i]
    cen_data +=data['cen_pix_noise_std'][i]

neigh_avg = neigh_data/len(flist)
cen_avg= cen_data/len(flist)
print("neigh avg standard devs: ",neigh_avg)
print("cen avg standard devs: ",cen_avg)
plt.imshow(neigh_avg,interpolation='nearest', cmap='gray',vmin = np.min(neigh_avg),vmax =np.max(neigh_avg)) 
plt.colorbar()
plt.title("Neigh AvgNoise STD Bgrms=.001")
plt.savefig("neigh_noise_std_bg001_2.png")
plt.close()

plt.imshow(cen_avg,interpolation='nearest', cmap='gray',vmin = np.min(cen_avg),vmax=np.max(cen_avg))
plt.colorbar()
plt.title("Cen AvgNoise STD Bgrms=.001")
plt.savefig("cen_noise_std_bg001_2.png")
#fitsio.write("avg_noise_std_neigh_bgrms10.fits", avg, clobber=True) 
