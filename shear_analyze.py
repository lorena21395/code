from glob import glob
import argparse
import fitsio
import numpy as np
from glob import glob
import esutil as eu

#parser = argparse.ArgumentParser()
#parser.add_argument("filename",help="filename")
#args = parser.parse_args()

# this gets a list of all files that match the pattern
flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run003/run003*.fits')

# read each file and combine into one big array
data = eu.io.read(flist)

#data = fitsio.read(args.filename)
w, = np.where(data['flags']==0)
print("kept %d/%d" % (w.size, data.size))
data = data[w]

g = data['pars'][:, 2:2+2]
g_1p = data['pars_1p'][:, 2:2+2]
g_1m = data['pars_1m'][:, 2:2+2]
g_2p = data['pars_2p'][:, 2:2+2]
g_2m = data['pars_2m'][:, 2:2+2]

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

frac = shear[0]/0.02-1
frac_err = shear_err[0]/0.02

print("bias: %g +/- %g" % (frac, frac_err))
print("additive: %g +/- %g" % (shear[1], shear_err[1]))

