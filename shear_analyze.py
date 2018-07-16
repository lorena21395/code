from glob import glob
import argparse
import fitsio
import numpy as np
import esutil as eu

#import matplotlib.pyplot as plt
#plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument("f1",help="filename path")
parser.add_argument("f2",help="filename path")
parser.add_argument("f3",help="filename path")
parser.add_argument("outfile",help="outfile name")
args = parser.parse_args()

f1 = args.f1
f2 = args.f2
f3 = args.f3
outfile = args.outfile

# this gets a list of all files that match the pattern
flist = glob(f1+'*.fits')
#flist = (glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run198/run198_1-output-*.fits'))
#flist.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run186/run186_21-output-007*.fits'))
flist2 = glob(f2+'*.fits')
#flist2.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run186/run186_22-output-006*.fits'))
#flist2.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run186/run186_22-output-007*.fits'))
flist3 = glob(f3+'*.fits')
#flist3.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run186/run186_16-output-002*.fits'))
#flist3.append(glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run186/run186_16-output-003*.fits'))
#flist = glob('/gpfs01/astro/workarea/lmezini/deblender_tests/code/test.fits')

lists = [flist,flist2,flist3]
# read each file and combine into one big array

shear1 = []
shear2 = []

noresp_shear1 = []
noresp_shear2 = []

shear1_err = []
shear2_err = []

noresp_shear1_err = []
noresp_shear2_err = []

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

    print(R11vals)
    R11 = R11vals.mean()
    R22 = R22vals.mean()

    R11_err = R11vals.std(axis=0)/np.sqrt(w.size)
    R22_err = R22vals.std(axis=0)/np.sqrt(w.size)

    shear = g_mean.copy()
    
    noresp_shear1.append(shear[0])
    noresp_shear2.append(shear[1])
    
    shear[0] /= R11
    shear[1] /= R22

    shear1.append(shear[0])
    shear2.append(shear[1])

    g_err = g.std(axis=0)/np.sqrt(w.size)
    shear_err = g_err.copy()

    noresp_shear1_err.append(shear_err[0])
    noresp_shear2_err.append(shear_err[1])
    
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
    S = noresp_shear1[2]/resp
    frac = S/0.02-1#-0.01
    frac_err = S_err/0.02

    print("bias: %g +/- %g" % (frac, frac_err))
    #print("additive: %g +/- %g" % (shear[1], shear_err[1]))
    print(resp)
    print(shear1,shear2)

dt = [('shear','f8'),('shear_err','f8')]
output = np.zeros(2,dtype = dt)
output['shear'],output['shear_err'] = frac,frac_err
fitsio.write(args.outfile,output,clobber=True)
"""
binwidth = 0.01
bins1 = np.arange(min(p1), max(p1) + binwidth, binwidth)
bins2 = np.arange(min(p2), max(p2) + binwidth, binwidth)

range_p1_min = p1.mean()-5.*p1.std(axis=0)
range_p1_max = p1.mean()+5.*p1.std(axis=0)
range_p2_min = p2.mean()-5.*p2.std(axis=0)
range_p2_max = p2.mean()+5.*p2.std(axis=0)
plt.hist(p1,bins1,histtype='step',label=("p1 mean,std: "+ str(round(p1.mean(),5))+"+/-" +str(round(p1_err,5))+", "+str(round(p1.std(axis=0),5))),range = (range_p1_min,range_p1_max))
plt.hist(p2,bins2,histtype='step',label=("p2 mean,std: "+ str(round(p2.mean(),5))+"+/-" +str(round(p2_err,5))+", "+str(round(p2.std(axis=0),5))),range = (range_p2_min,range_p2_min))
plt.legend(loc=0)
plt.title("Positions for Center Gal (R=11)")
plt.xlabel("position (pixel distance to center)")
plt.savefig('pos_r11.png')
plt.close()

binwidth = 0.01
bins11 = np.arange(min(R11vals), max(R11vals) + binwidth, binwidth)
bins22 = np.arange(min(R22vals), max(R22vals) + binwidth, binwidth)

range_R11_min = R11-5.*R11vals.std(axis=0)
range_R11_max = R11+5.*R11vals.std(axis=0)
range_R22_min = R22-5.*R22vals.std(axis=0)
range_R22_max = R22+5.*R22vals.std(axis=0)

plt.hist(R11vals,bins11,histtype='step',label=("R11 mean,std: "+ str(round(R11,3))+"+/-" +str(round(R11_err,3))+", "+str(round(R11vals.std(axis=0),3))),range=(range_R11_min,range_R11_max))
plt.hist(R22vals,bins22,histtype='step',label=("R22 mean,std: " + str(round(R22,3)) +"+/-" +str(round(R22_err,3))+", "+str(round(R22vals.std(axis=0),3))),range=(range_R22_min,range_R22_max))
plt.legend(loc = 0)
plt.title("Responses for R=11")
plt.xlabel("R")
plt.xlim(range_R11_min-0.02,range_R11_max+0.02)
plt.savefig('response_r11.png')
plt.close()

plt.hist(R11vals,bins11,histtype='step',label=("R11 mean,std: "+ str(round(R11,3))+"+/-" +str(round(R11_err,3))+", "+str(round(R11vals.std(axis=0),3))),range=(range_R11_min,range_R11_max))
plt.hist(R22vals,bins22,histtype='step',label=("R22 mean,std: " + str(round(R22,3)) +"+/-" +str(round(R22_err,3))+", "+str(round(R22vals.std(axis=0),3))),range=(range_R22_min,range_R22_max))
plt.legend(loc = 0)
plt.title("Responses for R=11")
plt.xlabel("R")
plt.xlim(range_R11_min-0.02,range_R11_max+0.02)
plt.yscale('log')
plt.savefig('response_r11_log.png')
plt.close()
"""
