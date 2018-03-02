from glob import glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import fitsio
import numpy as np
import esutil as eu

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run075/run*.fits')
data = eu.io.read(flist)
plt.hist(data,1000,range=[min(data), max(data)],histtype="step")
plt.title("Steps w/ erel=1e-5,bg10")
plt.savefig("steps_erel_00001_bg10.png")
