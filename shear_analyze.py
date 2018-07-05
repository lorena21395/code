from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu

#parser = argparse.ArgumentParser()
#parser.add_argument("filename",help="filename")
#args = parser.parse_args()
# this gets a list of all files that match the pattern
flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run192/run192_8-output-*.fits')
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run192/run192_12-output-*.fits'))
#flist2 = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run188/run188_2-*.fits')
#flist3 = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run188/run188_3-*.fits')
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test.fits')

lists = [flist]#,flist2,flist3]
# read each file and combine into one big array

shear1 = []
shear2 = []

shear1_err = []
shear2_err = []

for fl in lists:
    data = eu.io.read(fl)
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

    p1 = data['pars'][:,0]
    p2 = data['pars'][:,1]

    p1_err = p1.std(axis=0)/np.sqrt(w.size)
    p2_err = p2.std(axis=0)/np.sqrt(w.size)
    print("position: %g +/- %g" % (p1.mean(), p1_err))
    print("position: %g +/- %g" % (p2.mean(), p2_err))

    g_mean = g.mean(axis=0)


    # this is used to calibrate the shear, I will explain this
    R11vals = (g_1p[:,0] - g_1m[:,0])/0.02
    R22vals = (g_2p[:,1] - g_2m[:,1])/0.02

    R11 = R11vals.mean()
    R22 = R22vals.mean()

    shear = g_mean.copy()
    shear[0] /= R11
    shear[1] /= R22

    shear1.append(shear[0])
    shear2.append(shear[1])

    g_err = g.std(axis=0)/np.sqrt(w.size)
    shear_err = g_err.copy()
    shear_err[0] /= R11
    shear_err[1] /= R22

    shear1_err.append(shear_err[0])
    shear2_err.append(shear_err[1])

    frac = shear[0]/0.02-1#-0.01
    frac_err = shear_err[0]/0.02

    print("bias: %g +/- %g" % (frac, frac_err))
    print("additive e1: %g +/- %g" % (shear[0], shear_err[0]))
    print("additive e2: %g +/- %g" % (shear[1], shear_err[1]))


###if measuring response from deblending to shear ####
if len(lists) == 3:
    resp = (shear1[0]-shear1[1])/0.02
    S = shear1[2]/resp
    S_err = shear1_err[2]/resp
    frac = S/0.02-1#-0.01
    frac_err = S_err/0.02

    print("bias: %g +/- %g" % (frac, frac_err))
    #print("additive: %g +/- %g" % (shear[1], shear_err[1]))
    print(resp)
    print(shear1,shear2)
