import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.1)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.025)


r = [7,9,11,13,15,17,19]
"""
bias = [-0.0812011,0.00226535,0.0023871,-0.0123014,-0.00806869,0.00111734,0.004131]
err = [0.00063597,0.000587096,0.000554516,0.000540192,0.000538317,0.000536182,0.000535371]

bias_thresh005 = [-0.0640488,-0.0607346,-0.0263503,-0.00830626,0.00142396,0.006777,0.00594998]
err_thresh005 = [0.000683949,0.000593845,0.000558472,0.000543158,0.000538711,0.000536819,0.000535745]

bias_mini = [0.019005,-0.00361754,-0.00383634,-0.00105035,0.000678016,0.00128536,-0.00101486]
err_mini = [0.00273789,0.00103818,0.000826182,0.00080896,0.00028416,0.000804708,0.000566695]
"""

### low noise results
#bias = [-0.0196882,-0.00473535,-0.00195364,-0.00180363,-0.00582074,-0.00701391,-0.00510679]
#err = [0.000787389,0.000408626,0.000296113,0.000142127,0.000164697,0.000189147,0.000236306]

#bias_thresh005 = [-0.0486585,-0.0383526,-0.0123844,-0.00813976,-0.00756243,-0.00738492,-0.0054329]
#err_thresh005 = [0.000892989,0.000433803,0.000232912,0.000246704,0.000267562,0.000163117,0.000192444]

#bias_mini = [0.0157917,-0.602205,0.00166221,0.00054305,0.0114081,-0.000130317,0.000355017]
#err_mini = [0.0344728,0.0197704,0.000297506,0.00032547,0.00260028,0.000166784,9.60466e-05]

##bgrms = 0.1 Results

#bias = [0.0124162,0.00696045,-0.00701927,-0.0119736,-0.00758878,-0.00711046,-0.00307]
#err = [0.00879894,0.00131622,0.0013015,0.00129649,0.00129442,0.00129586,0.00141067]

#bias_thresh005 = [-0.0204691,-0.0265281,-0.0236849,-0.0156729,-0.0135327,-0.0110527,-0.00689438]
#err_thresh005 = [0.00141067,0.00132056,0.00130825,0.00129899,0.00129403,0.00129247,0.0012936]

#bias_mini = [0.015217,0.000734902,0.00168429,0.000579046,0.000622958,-0.0006141,-0.0015698]
#err_mini = [0.0107746,0.00187754,0.00130606,0.00130084,0.00129568,0.00129241,0.00129223]

##no noise re-added bgrms=10

#bias2 = [-0.0823561,0.00203271,0.00273231,-0.0118959,-0.00938221,0.00135266,0.00484167]
#err2 = [0.000609999,0.000572152,0.000547146,0.000538184,0.000537953,0.000535862,0.000534611]

#bias_thresh0052 = [-0.0647605,-0.0598775,-0.0268613,-0.00755885,0.00148567,0.00528927,0.00729856]
#err_thresh0052 = [0.000657842,0.000576765,0.00054923,0.00053998,0.000536613,0.000536545,0.000535982]

##no noise re-added bgrms = 0.1

#bias2 = [0.00450584,0.00832763,-0.00700911,-0.0122488,-0.00803526,-0.00730444,-0.00356617]
#err2 = [0.00136115,0.00131194,0.00130362,0.00130211,0.00129494,0.00129254,0.00129553]

#bias_thresh0052 = [-0.0177946,-0.0260808,-0.021333,-0.0173084,-0.010885,-0.0138278,-0.00560525]
#err_thresh0052 = [0.00140998,0.00132141,0.00130829,0.00129973,0.00129998,0.00129502,0.00129432]

##no second step, bgrms =10
#bias = [-0.198448,0.3075,0.256529,0.064536,0.0286707,0.0960161,0.0414348]
#err = [0.00196558,0.00280918,0.00237365,0.00180178,0.0015737,0.00122952,0.000689568]

#bias_thresh005 = [0.0187861,0.0445087,-0.0345543,0.00924453,-0.00162122,-0.0290263,0.0227536]
#err_thresh005 = [0.00173017,0.0017998,0.001966,0.0016946,0.00177524,0.00158451,0.00123914]

"""
#use old code
bias2 = [-0.0456805,0.0932809,0.0220271,-0.00466206,-0.00913659,0.0023212,-0.00034029]
err2 = [0.0051878,0.00507743,0.00489331,0.00462416,0.0040779,0.00421858,0.00396608]

bias_thresh0052 = [-0.0733619,-0.0694567,-0.0563895,-0.0243661,-0.00124443,0.00191828,-0.00440019]
err_thresh0052 = [0.00467769,0.00443167,0.00452556,0.00415059,0.00421953,0.00468915,0.0041015]
"""
"""
##Sersic gals
"""
bias = [-0.119586,-0.041375,-0.0157337,-0.0191576,-0.0236748,-0.0171168,-0.0142049]
err = [0.00152048,0.00145598,0.0013534,0.00131648,0.00132543,0.00132846,0.00114863]

#bias_mini = [-0.106728,-0.0188556,-0.00881271,-0.00431783,-0.00182994,0.00150557,-0.000325426]
#err_mini = [0.00692365,0.00187336,0.00141832,0.00133865,0.0013372,0.00133648,0.0013343]

#used cm model for mini

bias_mini = [-0.0167178,0.29282,0.0556112,0.00428428,-0.00447817,-0.00324964,-0.00286702]
err_mini = [0.00511843,0.00418356,0.00219845,0.00108241,0.000943684,0.000944699,0.000942955]

## Sersic with bgrms = 0.1
#bias = [-0.0652111,-0.0304299,-0.0165877,-0.00791151,-0.00485946,-0.00242685,-0.00582617]
#err = [0.00230512,0.00227304,0.00226535,0.00225719,0.00225828,0.00227093,0.00226765]

#bias_mini = [-0.159687,-0.0164582,-0.00409156,0.00378743,4.3599e-05,0.000101924,0.000573448]
#err_mini = [0.035359,0.00488523,0.00230586,0.00229654,0.00230798,0.0022896,0.00227844]

#used cm model for mini
#bias_mini = [-0.169015,0.267239,-0.00207964,-0.000695256,-0.000599733,-0.000912097,-0.000706023]
#err_mini = [0.0262119,0.0192209,0.00177944,0.00160845,0.00160396,0.00160149,0.00160665]

fig,ax = plt.subplots()

plt.axhline(y=0,color='black',linewidth=2)
plt.errorbar(r,bias,yerr=err,xerr=None,label="Scarlet",color='r')
#plt.errorbar(r,bias2,yerr=err2,xerr=None,label="Scarlet",linestyle = ":",color='r')
plt.errorbar(r,bias_mini,yerr=err_mini,xerr=None,label="Minimof",color='b')
#plt.errorbar(r,bias_thresh005,yerr=err_thresh005,xerr=None,label="Scarlet (EFT=0.05)",color='m')
#plt.errorbar(r,bias_thresh0052,yerr=err_thresh0052,xerr=None,label="Scarlet (EFT=0.05)",color='m',linestyle=":")
plt.xlabel("Distance to Neighbor (pixels)")
plt.ylabel("Fractional Shear Bias")
plt.title("Bias vs Radius (N=10, Sersic)")
#plt.yticks(np.arange(min(bias), max(bias)+0.01, .005))
plt.grid(which='minor')
plt.grid(which='major')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
plt.legend()
# for the minor ticks, use no labels; default NullFormatter
ax.yaxis.set_minor_locator(minorLocator)
#plt.ylim(-0.075,0.075)
plt.show()
plt.savefig("bias_v_rad_n10_sers_mini_scar2.png")
