from glob import glob
import argparse
from astropy.io import fits
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run195/run195_4-*')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/test.fits')
# read each file and combine into one big array
mean_data = np.zeros((25,25))
std_data = np.zeros((25,25))
data = eu.io.read(flist)
for i in range(len(flist)):
    std_data +=data['std'][i]
    mean_data +=data['mean'][i]
std_avg = std_data/len(flist)
mean_avg= mean_data/len(flist)
print("avg standard devs: ",std_avg)
print("avg means: ",mean_avg)


f, ax = plt.subplots(1,3,figsize=(12,4))
f1=ax[0].imshow(std_avg,interpolation='nearest', cmap='gray',vmin = np.min(std_avg),vmax=np.max(std_avg)) 
plt.colorbar(f1,ax=ax[0])
ax[0].set_title("True-Mod STD")
#plt.savefig("multi_true-mod_std_r9_PSF_match_s02_3.png")
#plt.close()

f2 = ax[1].imshow(mean_avg,interpolation='nearest', cmap='gray',vmin =np.min(mean_avg),vmax=np.max(mean_avg))
plt.colorbar(f2,ax=ax[1])
ax[1].set_title("True-Mod Mean")
#plt.savefig("multi_true-mod_mean_r9_PSF_match_s02_3.png")
#plt.close()

ax[2].hist((std_avg.flatten()),100,histtype='step')
ax[2].set_title("True-Mod Stds")
plt.tight_layout()
#f.savefig("multi_true-mod_sers_r9_n10_size51_noSC.png")
f.savefig('PSF_true-mod_noPSF.png')
#f.savefig('test.png')
plt.close()

mean_of_stds = std_avg.mean()
err_of_mean_stds = std_avg.std()/np.sqrt(std_avg.size)

print("mean: %g +/- %g" % (mean_of_stds,err_of_mean_stds))
"""
hdu = (fits.PrimaryHDU(mean_avg))
hdu1 = fits.HDUList([hdu])
hdu1.writeto('true-mod-r9-N10-mean_avg.fits')

hdu = (fits.PrimaryHDU(std_avg))
hdu1 = fits.HDUList([hdu])
hdu1.writeto('true-mod-r9-N10-std_avg.fits')
"""
