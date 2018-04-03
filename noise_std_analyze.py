from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run119/run*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/test.fits')
# read each file and combine into one big array
neigh_data = np.zeros((15,15))
cen_data = np.zeros((15,15))
data = eu.io.read(flist)
for i in range(len(flist)):
    neigh_data +=data['std'][i]
    cen_data +=data['mean'][i]
neigh_avg = neigh_data/len(flist)
cen_avg= cen_data/len(flist)
print("neigh avg standard devs: ",neigh_avg)
print("cen avg standard devs: ",cen_avg)

plt.imshow(cen_avg,interpolation='nearest', cmap='gray',vmin = np.min(cen_avg),vmax =np.max(cen_avg)) 
plt.colorbar()
plt.title("Single Object True-Mod Mean")
plt.savefig("sing_true-mod_mean.png")
plt.close()

plt.imshow(neigh_avg,interpolation='nearest', cmap='gray',vmin = np.min(neigh_avg),vmax=np.max(neigh_avg))
plt.colorbar()
plt.title("Single Object True-Mod Std")
plt.savefig("sing_true-mod_std.png")
plt.close()

plt.hist((neigh_avg.flatten()),histtype='step')
plt.title("Single Object True-Mod Stds")
plt.savefig("sing_true-mod_hist.png")
plt.close()

mean_of_stds = neigh_avg.mean()
err_of_mean_stds = neigh_avg.std()/np.sqrt(neigh_avg.size)

print("mean: %g +/- %g" % (mean_of_stds,err_of_mean_stds))
