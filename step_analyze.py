from glob import glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import fitsio
import numpy as np
import esutil as eu

#flist = glob('test*.fits')
flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run078/run*.fits')
data = eu.io.read(flist)
print(data)
init_step = data['init_step']
fin_step = data['final_step']
#mid_step = data['mid_step']
print(len(init_step),len(fin_step))
w = np.where(data['init_step']== 10000)
print(len(w[0])/1000.)
#w = np.where(data['mid_step']== 10000)
#print(len(w[0])/1000.)
w = np.where(data['final_step'] == 10000)
print(len(w[0])/1000.)


#print(data['init_step'])
#print(data['final_step'])
plt.hist(init_step,1000,range=[min(init_step), max(init_step)],histtype="step")
plt.title("Init Steps w/ erel=1e-5,1e4 max")
plt.savefig("steps_init_erel_e-5_maxe4_bg10_2.png")
plt.close()
plt.hist(fin_step,1000,range=[min(fin_step), max(fin_step)],histtype="step")
plt.title("Fin Steps w/ erel=1e-5,1e4 max")
plt.savefig("steps_fin_erel_e-5_maxe4_bg10_2.png")
plt.close()
#plt.hist(mid_step,1000,range=[min(mid_step), max(mid_step)],histtype="step")
#plt.title("Init Steps w/ erel=1e-5,1e4 max")
#plt.savefig("steps_mid_erel_e-5_maxe4_bg10.png")
#plt.close()
