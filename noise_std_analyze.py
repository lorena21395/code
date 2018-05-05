from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run166/run166_1*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/code/test.fits')
# read each file and combine into one big array
mean_data = np.zeros((15,15))
std_data = np.zeros((15,15))
data = eu.io.read(flist)
for i in range(len(flist)):
    std_data +=data['std'][i]
    mean_data +=data['mean'][i]
std_avg = std_data/len(flist)
mean_avg= mean_data/len(flist)
print("avg standard devs: ",std_avg)
print("avg means: ",mean_avg)

plt.imshow(std_avg,interpolation='nearest', cmap='gray',vmin = np.min(std_avg),vmax =np.max(std_avg)) 
plt.colorbar()
plt.title("Multi Object True-Mod STD r=9")
plt.savefig("multi_true-mod_std_r9_PSF_match.png")
plt.close()

plt.imshow(mean_avg,interpolation='nearest', cmap='gray',vmin = np.min(mean_avg),vmax=np.max(mean_avg))
plt.colorbar()
plt.title("Multi Object True-Mod Mean r=9")
plt.savefig("multi_true-mod_mean_r9_PSF_match.png")
plt.close()

plt.hist((std_avg.flatten()),100,histtype='step')
plt.title("Multi Object True-Mod Stds r=9")
plt.savefig("multi_true-mod_hist_r9_PSF_match.png")
plt.close()

mean_of_stds = std_avg.mean()
err_of_mean_stds = std_avg.std()/np.sqrt(std_avg.size)

print("mean: %g +/- %g" % (mean_of_stds,err_of_mean_stds))
