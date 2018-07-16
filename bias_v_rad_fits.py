import numpy as np
import matplotlib.pyplot as plt
import fitsio
import esutil as eu

plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
"""
parser = argparse.ArgumentParser()
parser.add_argument("f1",help="filename path")
parser.add_argument("f2",help="filename path")
parser.add_argument("f3",help="filename path")
parser.add_argument("f4",help="filename path")
parser.add_argument("f5",help="filename path")
parser.add_argument("outfile",help="outfile name")
args = parser.parse_args()

f1 = args.f1
f2 = args.f2
f3 = args.f3
f4 = args.f4
f5 = args.f5

fl1 = glob(f1+'*.fits')
data1 = eu.io.read(fl1)
fl2 = glob(f2+'*.fits')
data2 = eu.io.read(fl2)
fl3 = glob(f3+'*.fits')
data3 = eu.io.read(fl3)
fl4= glob(f4+'*.fits')
data4 = eu.io.read(fl4)
fl5= glob(f5+'*.fits')
data5 = eu.io.read(fl5)

bias = [data1['Shear'],data2['Shear'],data3['Shear'],data4['Shear'],data5['Shear']]
err = [data['Shear_err'],data['Shear_err'],data['Shear_err'],data['Shear_err'],data['Shear_err']]
"""
majorLocator = MultipleLocator(0.001)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.0005)

r = [7,9,11,13,15]

bias = [-0.00178431,0.00419539,-0.00282502,-0.000611995,0.000631424]
err = [0.00230481,0.0021461,0.00206198,0.00290171,0.00290287]

bias2 = [-0.00268918,0.00408031,-0.00366705,-0.00062781,7.7428e-05]
err2 = [0.00230272,0.00214585,0.00206024,0.00290167,0.00290126]

fig,ax = plt.subplots()
plt.axhline(y=0,color='black',linewidth=2)
plt.errorbar(r,bias,yerr=err,xerr=None,label="no meta-resp",color='r')
plt.errorbar(r,bias2,yerr=err2,xerr=None,label="w/ meta-resp",color='r',linestyle=':')
plt.xlabel("Distance to Neighbor (pixels)")
plt.ylabel("Fractional Shear Bias")
plt.title("Bias v radius")

plt.grid(which='minor')
plt.grid(which='major')

ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
ax.set_xlim(6,16)
plt.legend(loc=4)

#plt.savefig('bias_v_rad_MOF_scar_knots_N0.1.png')                              
plt.savefig('bias-v-rad-mof-sers-N10.png')
