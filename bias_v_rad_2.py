import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator = MultipleLocator(0.01)
majorFormatter = FormatStrFormatter('%s')
minorLocator = MultipleLocator(0.005)

r = [7,9,11,13,15]

"""
####Calibrate with extra shear response:
## ALL SERSIC GALS
#2 step, bg =0.1
bias = [0.00266184,-0.0248007,-0.00544136,-0.00624419,0.00708624]
err = [0.00127191,0.00113552,0.00107044,0.00104095,0.00138829]
#2 step, bg=0.1, without division by metacal responses:
bias2 = [0.00289931,-0.0227334,-0.00450081,-0.00345485,0.00987126]
err2 = [0.00127221,0.00113792,0.00107146,0.00104387,0.00139212]

#1step, bg=0.1
bias9 = [0.00183218,-0.0135421,-0.0131184,-0.00579225,0.00621719]
err9 = [0.00170727,0.0019614,0.00137348,0.00134146,0.00137743]

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
"""
###regular mof bg= 0.1
bias2 = [-0.000237279,-0.000476259,-0.000621423,-0.00105133,-6.75232e-05]
err2 = [0.000735791,0.000724539,0.000717248,0.000713775,0.0010105]

###regular mof bg = 10
bias4 = [-0.0245964,-0.0129146,-0.00877397,0.000550745,0.00151485]
err4 = [0.0031733,0.00298901,0.00315932,0.00290509,0.00290543]


"""
############
#knots only (MOF)
###########
#N = 0.1
bias = [-0.00906016,0.00320047,0.000231944,0.00544578,-0.00609954]
err = [0.00364751,0.00333022,0.0030994,0.00296978,0.00291995]
#N = 10
#bias2 = [-0.00691184,-0.0016747,-0.00334637,3.37054e-05,0.00701863]
#err2 = [0.00416599,0.00391162,0.00362538,0.00347558,0.00342445]

#KNOTS ONLY (SCARLET)
#N = 0.1
#1step
bias3 = [-0.0181289,-0.0993413,-0.0641069,-0.0312261,-0.00394522]
err3= [0.00376232,0.00327911,0.0030641,0.00296848,0.00293451]

#2step
bias4 = [-0.0150497,-0.0860171,-0.0518112,-0.0213856,-0.00203856]
err4 = [0.00358457,0.00320502,0.0030408,0.00296391,0.00293567]
"""
"""
###Trying different constraints###
###Scarlet, sersic type galaxies, bgrms = 0.1###

#default psf constraints, no source constraints
bias = [-0.245301,-0.299176,-0.235179,-0.200164,-0.129979] 
err = [0.00195663,0.0013471,0.00111105,0.00113133,0.00154654]

#no psf constraints, no source constraints
bias2 = [-0.246145,-0.299217, -0.235316,-0.200361,-0.130444]
err2 = [0.00196074,0.00134868,0.00111116,0.0011295,0.00154619]

#default psf constraints, simple source constraint
bias3 = [-0.235998,-0.213525,-0.122551,-0.0578107,-0.0319187]
err3 = [0.00214287,0.00174524,0.00204687,0.0013753,0.00113997]

#no psf constraints, simple source constraint
bias4 = [-0.239269,-0.213929,-0.122214,-0.0575771,-0.0317016]
err4 = [0.00218168,0.0017496,0.00207052,0.00137847,0.00114131]

#default psf constraints, simple + monotonicty for sources
bias5 = [-0.269635,-0.205867,-0.118536,-0.0365999, 0.00252565]
err5 = [0.00157017,0.0010794,0.00110085,0.00104087,0.00109344]

#none for psf, simple + monotonicity for sources
bias6 = [-0.269641,-0.205904,-0.11858,-0.0366356,0.00262338]
err6 = [0.0015712,0.00108116,0.00110036,0.00104084,0.00109373]

#default for psf, simple and symmetry for sources
bias7 = [-0.12889,-0.104401,-0.0714764,-0.0809106,-0.0231729]
err7 = [0.00110826,0.00103729,0.00104761,0.00102487,0.00105384]

#none for psf, simple and symmetry for sources
bias8 = [-0.132263,-0.105286,-0.0725884,-0.0803108,-0.0224682]
err8 = [0.00110814,0.00103429,0.00104042,0.00102166,0.00105028]

#simple for psf, default for sources

bias9 = [-0.184006,-0.116358,-0.0654654,-0.0229924,0.011997]
err9 = [0.00139057,0.00175697,0.00130063,0.00131826,0.00138534]
"""

#####
#MOF 2band polarizations (6/28/2018)
####
bias = [0.00272439,0.00140518,-0.000214843,-1.30061e-05,-0.000252071]
err = [0.00104082,0.00102796,0.00101565,0.00101317,0.00101174]

#RGB bands (6/29/2018)
bias3 = [0.000707012,0.00112226,0.000425393,0.000560498,0.000993042]
err3 = [0.00104322,0.00102752,0.00101707,0.00101164,0.00101319]

#RGB bands with bgrms = 10
bias5 = [-0.0873927,-0.0698742,-0.0547051,-0.052327,-0.0563293]
err5 = [0.00675345,0.00645022,0.00640499,0.00647275,0.00647933]
fig,ax = plt.subplots()
plt.axhline(y=0,color='black',linewidth=2)
plt.errorbar(r,bias,yerr=err,xerr=None,label="MOF (2 band N0.1)",color='b')
plt.errorbar(r,bias2,yerr=err2,xerr=None,label="MOF (1 band N0.1)",color='b',linestyle=':')
plt.errorbar(r,bias3,yerr=err3,exerr=None,label="MOF (RGB N0.1)",color='b',linestyle = '--')
plt.errorbar(r,bias4,yerr=err4,xerr=None,label="MOF (1 band N10)",color='m',linestyle=':')
#plt.errorbar(r,bias6,yerr=err6,xerr=None,label="None/Simp+mono",color='m')
plt.errorbar(r,bias5,yerr=err5,xerr=None,label="MOF (RGB N10)",color='m',linestyle='--')
#plt.errorbar(r,bias7,yerr=err7,xerr=None, label= "Default/Simp+sym",color='c',linestyle=':')
#plt.errorbar(r,bias8,yerr=err8,xerr=None, label= "None/Simp+sym",color='c')
#plt.errorbar(r,bias9,yerr=err9,xerr=None, label= "Simple/Default",color='g')
#plt.errorbar(13,-0.00125321,yerr=0.00102325,xerr=None,label = "MOF w/ cali + sep",marker="o",color="c")
#plt.errorbar(13,0.000438874,yerr=0.00102498,xerr=None,label="MOF + sep",marker="o",color="r")

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
plt.savefig('bias-v-rad-MOF-multiband.png')
