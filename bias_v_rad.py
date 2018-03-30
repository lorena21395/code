import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.01)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.0025)


r = [7,9,11,13,15,17,19]
bias = [-0.0798784,0.00250587,0.00246968,-0.0117825,-0.00859576,0.00115493,0.00448714]
err = [0.000451956,0.0004099,0.000381296,0.000368168,0.000365477,0.000363522,0.000362772]

fig,ax = plt.subplots()

plt.errorbar(r,bias,yerr=err,xerr=None)
plt.xlabel("Neighbor distance from center (pixels)")
plt.ylabel("Bias")
plt.title("Bias vs Radius")
#plt.yticks(np.arange(min(bias), max(bias)+0.01, .005))
plt.grid(which='minor')
plt.grid(which='major')y
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)

# for the minor ticks, use no labels; default NullFormatter
ax.yaxis.set_minor_locator(minorLocator)

plt.show()
plt.savefig("bias_v_rad.png")
