from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu

#parser = argparse.ArgumentParser()
#parser.add_argument("filename",help="filename")
#args = parser.parse_args()

# this gets a list of all files that match the pattern
flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run160/run160_7-*.fits')
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run155/run155_1*.fits'))
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run098/run098-output-00002*.fits'))
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run05*/run*.fits'))
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run06*/run*.fits'))
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run07*/run*.fits'))
#flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run029/run029*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test.fits')
# read each file and combine into one big array
data = eu.io.read(flist)
#data = fitsio.read(args.filename)

w, = np.where(data['flags']==0)
#w, = np.where(data['mod_size_flag']==0)
print("kept %d/%d" % (w.size, data.size))
data = data[w]

g = data['pars'][:, 2:2+2]
g_1p = data['pars_1p'][:, 2:2+2]
g_1m = data['pars_1m'][:, 2:2+2]
g_2p = data['pars_2p'][:, 2:2+2]
g_2m = data['pars_2m'][:, 2:2+2]

e1 = data['pars'][:,2]
e2 = data['pars'][:,3]

dims = data['dims']
#asym = np.where(dims[0]!=dims[1])
#print(dims[asym])
#print(np.min(dims))
#print(np.mean(dims))
#print(np.max(dims))
g_mean = g.mean(axis=0)

# this is used to calibrate the shear, I will explain this
R11vals = (g_1p[:,0] - g_1m[:,0])/0.02
R22vals = (g_2p[:,1] - g_2m[:,1])/0.02

R11 = R11vals.mean()
R22 = R22vals.mean()

shear = g_mean.copy()
shear[0] /= R11
shear[1] /= R22

g_err = g.std(axis=0)/np.sqrt(w.size)
shear_err = g_err.copy()
shear_err[0] /= R11
shear_err[1] /= R22

frac = shear[0]/0.1-1-0.01
frac_err = shear_err[0]/0.1

print("bias: %g +/- %g" % (frac, frac_err))
print("additive: %g +/- %g" % (shear[1], shear_err[1]))
