import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.005)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.0025)

s = [0,0.005,0.01,0.015,0.02]

#e1 component
#1step
addit = [-4.00162e-05,0.00435731,0.00872995,0.0131275,0.0175651]
err = [1.95307e-05,1.96184e-05,2.48956e-05,2.0034e-05,2.0357e-05]
#2step
addit2 = [2.1547e-06,0.00432599,0.00874227,0.013188,0.0175884]
err2 = [1.97703e-05,2.04646e-05,1.99056e-05,2.0173e-05,2.0449e-05]

#e2 componenet
addit3 = [-2.30634e-06,0.004325,0.00860673,0.0128571,0.0171364]
err3 = [1.9362e-05,1.94542e-05,1.9595e-05,1.99103e-05,2.02676e-05]
#x = np.array(addit)/np.array(s)-1
#err = np.array(err)/np.array(s)

#x2 = np.array(addit2)/np.array(s)-1
#err2 = np.array(err2)/np.array(s)

fig,ax = plt.subplots()

plt.errorbar(s,addit,yerr=err,xerr=None,label='e1 shear',color='r')
plt.errorbar(s,addit3,yerr=err3,xerr=None,label='e2 shear',color='r',linestyle=':')
plt.plot(s,s,label='True',color='b')

#plt.axhline(y=0,color='black',linewidth=2)
#plt.errorbar(s,x,yerr=err,xerr=None,label="1step",color='r')
#plt.errorbar(s,x2,yerr=err2,xerr=None,label="2-step",color='r',linestyle=':')
plt.xlabel("e1,e2 shear")
#plt.ylabel("Fractional Shear Bias")
plt.ylabel("additive")
plt.title("True shear v measured shear (r=9)")

plt.grid(which='minor')
plt.grid(which='major')

ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
#ax.set_xlim(6,16)
plt.legend(loc=4)

plt.savefig('test.png')
