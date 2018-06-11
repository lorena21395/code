from glob import glob
import fitsio
import numpy as np
import esutil as eu
import argparse
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#parser = argparse.ArgumentParser()                                             
#parser.add_argument("filename",help="filename")                                
#args = parser.parse_args()

# this gets a list of all files that match the pattern                          
flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run171/run171_21-*.fits')

data = eu.io.read(flist)
#data = fitsio.read(args.filename)                                              

w, = np.where(data['flags']==0)
#w, = np.where(data['mod_size_flag']==0)                                        
print("kept %d/%d" % (w.size, data.size))
data = data[w]
#y, = np.where(data['flags']==1)                                                
#print("metacal flag",y.size)                                                   

e1 = data['pars'][:,2]
e2 = data['pars'][:,3]

plt.hist(e1,50,histtype='step',label = 'pars 2')
plt.hist(e2,50,histtype='step',label = 'pars 3')
plt.legend()
plt.title("Ellipticities for minimof r=19 (N=0.1)")
plt.xlabel("Ellipticity")

plt.savefig("ellipticities_mini_r19_n01.png")
