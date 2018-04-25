import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.01)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.005)


r = [7,9,11,13,15,17,19]
#bias = [-0.0798784,0.00250587,0.00246968,-0.0117825,-0.00859576,0.00115493,0.00448714]
#err = [0.000451956,0.0004099,0.000381296,0.000368168,0.000365477,0.000363522,0.000362772]

#bias = [-0.0777718,0.00361787,0.00249013,-0.0119868,-0.0079754,0.000218821,0.00288052]
#err = [0.00100841,0.000915445,0.000849632,0.000823646,0.00082008,0.000810255,0.000810833]

#bias_thresh005 = [-0.0656239,-0.0582427,-0.0293778,-0.00732776,0.00196363,0.00614466,0.00710499]
#err_thresh005 = [0.00109913,0.000929064,0.000858764,0.000825089,0.000816381,0.000814174,0.000813136]

#bias = [-0.0988779,-0.00859686,-0.00648767,-0.0253978,-0.017486,-0.00946224,-0.00390392]
#err = [0.00171358,0.00166918,0.00162093,0.00160796,0.00161678,0.00161748,0.00161722]

#bias = [-0.0985608,-0.00870379,-0.00439757,-0.0262299,-0.0204152,-0.00927668,-0.00583736]
#err = [0.00115685,0.00105166,0.000976681,0.000940484,0.000932757,0.000933138,0.000927936]

#bias = [-0.0981212,-0.00811703,-0.00595669,-0.0257535,-0.0208248,-0.00696838,-0.00558925]
#err = [0.00115576,0.00105386,0.000976345,0.000940587,0.000934544,0.000931826,0.000928647]

#bias_mini = [0.019005,-0.00361754,-0.00383634,-0.00105035,0.000678016,0.00128536,-0.00101486]
#err_mini = [0.00273789,0.00103818,0.000826182,0.00080896,0.00028416,0.000804708,0.000566695]


### low noise results
bias = [-0.0196882,-0.00473535,-0.00195364,-0.00180363,-0.00582074,-0.00701391,-0.00510679]
err = [0.000787389,0.000408626,0.000296113,0.000142127,0.000164697,0.000189147,0.000236306]

bias_thresh005 = [-0.0486585,-0.0383526,-0.0123844,-0.00813976,-0.00756243,-0.00738492,-0.0054329]
err_thresh005 = [0.000892989,0.000433803,0.000232912,0.000246704,0.000267562,0.000163117,0.000192444]

bias_mini = [0.0157917,-0.602205,0.00166221,0.00054305,0.0114081,-0.000130317,0.000355017]
err_mini = [0.0344728,0.0197704,0.000297506,0.00032547,0.00260028,0.000166784,9.60466e-05]


fig,ax = plt.subplots()

plt.axhline(y=0,color='black',linewidth=2)
plt.errorbar(r,bias,yerr=err,xerr=None,label="Scarlet",color='r')
plt.errorbar(r,bias_mini,yerr=err_mini,xerr=None,label="Minimof",color='b')
plt.errorbar(r,bias_thresh005,yerr=err_thresh005,xerr=None,label="Scarlet (EFT=0.05)",color='m')
plt.xlabel("Distance to Neighbor (pixels)")
plt.ylabel("Fractional Shear Bias")
plt.title("Bias vs Radius (low noise)")
#plt.yticks(np.arange(min(bias), max(bias)+0.01, .005))
plt.grid(which='minor')
plt.grid(which='major')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
plt.legend()
# for the minor ticks, use no labels; default NullFormatter
ax.yaxis.set_minor_locator(minorLocator)
plt.ylim(-0.075,0.075)
plt.show()
plt.savefig("bias_v_rad_mm_scar_ln.png")
