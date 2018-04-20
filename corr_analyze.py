from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu
import scarlet.display
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import mlab, cm
import scipy.ndimage

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run150/run*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test*.fits')

#corr_big = np.zeros((49,49))
corr = np.zeros((29,29))

data = eu.io.read(flist)

for i in range(len(flist)):
    #corr_big +=data['corr_big'][i]
    corr +=data['corr'][i]
#corr_big_avg = corr_big/len(flist)
corr_avg= corr/len(flist)
"""
x = np.arange(0., 49., 1)
y = np.arange(0., 49., 1)
X, Y = np.meshgrid(x, y)
norm = cm.colors.Normalize(vmax=abs(corr_big_avg).max()-50000, vmin=corr_big_avg.min())
cmap = cm.PRGn
cset1 = plt.contourf(X, Y,corr_big_avg,cmap=cm.get_cmap(cmap),norm=norm)
cset2 = plt.contour(X, Y, corr_big_avg, cset1.levels, colors='k')
asinh = scarlet.display.Asinh(img=corr_big_avg, Q=100)

plt.imshow(corr_big_avg,cmap=cmap,norm=norm)

plt.imshow(corr_big_avg,interpolation='nearest', cmap='gray',vmin = np.min(corr_big_avg),vmax =np.max(15))#corr_big_avg))#,norm=asinh)
#plt.contour(corr_big_avg, colors='k', origin='upper')
plt.colorbar()
plt.title("Averaged Noise Correlation")
plt.savefig('test.png')
plt.close()
"""
plt.plot(corr_avg[14,:])
#plt.ylim(-15,15)
plt.tight_layout()
plt.savefig('corr_cross_section.png')


plt.imshow(corr_avg,interpolation='nearest', cmap='gray',vmin = np.min(corr_avg),vmax =np.max(corr_avg))#corr_big_avg))#,norm=asinh)
plt.colorbar()
plt.title("Averaged Noise Correlation")
plt.savefig('avg_noise_corr.png')

