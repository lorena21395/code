import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.01)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.0025)


r = [7,9,11,13,15,17,19]
#bias = [-0.0798784,0.00250587,0.00246968,-0.0117825,-0.00859576,0.00115493,0.00448714]
#err = [0.000451956,0.0004099,0.000381296,0.000368168,0.000365477,0.000363522,0.000362772]

#bias = [-0.0777718,0.00361787,0.00249013,-0.0119868,-0.0079754,0.000218821,0.00288052]
#err = [0.00100841,0.000915445,0.000849632,0.000823646,0.00082008,0.000810255,0.000810833]

#bias = [-0.0656239,-0.0582427,-0.0293778,-0.00732776,0.00196363,0.00614466,0.00710499]
#err = [0.00109913,0.000929064,0.000858764,0.000825089,0.000816381,0.000814174,0.000813136]

#bias = [-0.0988779,-0.00859686,-0.00648767,-0.0253978,-0.017486,-0.00946224,-0.00390392]
#err = [0.00171358,0.00166918,0.00162093,0.00160796,0.00161678,0.00161748,0.00161722]

#bias = [-0.0985608,-0.00870379,-0.00439757,-0.0262299,-0.0204152,-0.00927668,-0.00583736]
#err = [0.00115685,0.00105166,0.000976681,0.000940484,0.000932757,0.000933138,0.000927936]

bias = [-0.0981212,-0.00811703,-0.00595669,-0.0257535,-0.0208248,-0.00696838,-0.00558925]
err = [0.00115576,0.00105386,0.000976345,0.000940587,0.000934544,0.000931826,0.000928647]

fig,ax = plt.subplots()

plt.errorbar(r,bias,yerr=err,xerr=None)
plt.xlabel("Neighbor distance from center (pixels)")
plt.ylabel("Bias")
plt.title("Bias vs Radius (Exp,w/ |0)")
#plt.yticks(np.arange(min(bias), max(bias)+0.01, .005))
plt.grid(which='minor')
plt.grid(which='major')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)

# for the minor ticks, use no labels; default NullFormatter
ax.yaxis.set_minor_locator(minorLocator)

plt.show()
plt.savefig("bias_v_rad_exp_WC.png")
