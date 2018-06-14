import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.002)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.001)

r = [7,9,11,13,15]
####Calibrate with extra shear response:
## ALL SERSIC GALS
#2 step, bg =0.1
bias = [0.00266184,-0.0248007,-0.00544136,-0.00624419,0.00708624]
err = [0.00127191,0.00113552,0.00107044,0.00104095,0.00138829]
#2 step, bg=0.1, without division by metacal responses:
bias2 = [0.00289931,-0.0227334,-0.00450081,-0.00345485,0.00987126]
err2 = [0.00127221,0.00113792,0.00107146,0.00104387,0.00139212]

#1step, bg=0.1
bias3 = [0.00183218,-0.0135421,-0.0131184,-0.00579225,0.00621719]
err3 = [0.00170727,0.0019614,0.00137348,0.00134146,0.00137743]
#1step, bg = 0.1, w/o division by metacal response
bias4 = [0.00303388,-0.0110562,-0.0123844,-0.00260066,0.0092246]
err4 = [0.00170932,0.00196634,0.0013745,0.00134577,0.00138155]

#MOF, bg = 0.1: 
bias5 = [-0.00283521,-0.0035636,0.00294756,0.00142434,-0.00063632]
err5 = [0.0010371,0.00102218,0.00101716,0.00101186,0.00100993]
#MOF, bg = 0.1, w/o division by metacal response
bias6 = [ -0.00249482,-0.000580047,0.00116197,0.00160032,-0.000313858]
err6 = [0.00073413,0.000724463,0.000718528,0.000715669,0.000778343]

#MOF, bg = 10 w/o division by metacal response:
bias7 = [0.00130908,0.00398234,-0.0133561,-0.0004764,-0.00127988] 
err7 = [0.00325758,0.00304018,0.00314472,0.00290211,0.00289733]

bias8 = [0.000719619,0.00322232,-0.00822339,-0.000963111,-0.00186173]
err8 = [0.00325566,0.00303788,0.00316108,0.00290069,0.00289564]

###regular mof bg= 0.1
bias = [-0.000237279,-0.000476259,-0.000621423,-0.00105133,-6.75232e-05]
err = [0.000735791,0.000724539,0.000717248,0.000713775,0.0010105]
###regular mof bg = 10
#bias = [-0.0245964,-0.0129146,-0.00877397,0.000550745,0.00151485]
#err = [0.0031733,0.00298901,0.00315932,0.00290509,0.00290543]
fig,ax = plt.subplots()
plt.axhline(y=0,color='black',linewidth=2)
#plt.errorbar(r,bias,yerr=err,xerr=None,label="2-step w/ meta",color='r')
#plt.errorbar(r,bias2,yerr=err2,xerr=None,label=" 2-Step",color='r',linestyle=':')
#plt.errorbar(r,bias3,yerr=err3,xerr=None,label="1-Step w/ meta",color='blue')
#plt.errorbar(r,bias4,yerr=err4,xerr=None,label="1-Step",color='blue',linestyle=':')
plt.errorbar(r,bias,yerr=err,xerr=None,label="MOF",color='m',linestyle=':')
plt.errorbar(r,bias6,yerr=err6,xerr=None, label= "MOF w/ cali",color='m')
#plt.errorbar(r,bias5,yerr=err5,xerr=None, label= "MOF w/ cali, no resp",color='m')
#plt.errorbar(r,bias8,yerr=err8,xerr=None, label= "MOF w/ cali",color='m',linestyle = "-")
plt.errorbar(13,-0.00125321,yerr=0.00102325,xerr=None,label = "MOF w/ cali + sep",marker="o",color="c")
plt.errorbar(13,0.000438874,yerr=0.00102498,xerr=None,label="MOF + sep",marker="o",color="r")
plt.xlabel("Distance to Neighbor (pixels)")
plt.ylabel("Fractional Shear Bias")
plt.title("Bias v radius (N 0.1)")

plt.grid(which='minor')
plt.grid(which='major')

ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)

plt.legend(loc=4)

#plt.savefig('bias_v_rad_N0.1_scar_mof_cali.png')
plt.savefig('test.png')
