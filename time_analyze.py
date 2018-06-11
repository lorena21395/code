from glob import glob
import fitsio
import numpy as np
import esutil as eu
import matplotlib.pyplot as plt                                                
plt.switch_backend('agg')

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run177/run177_2-*.fits')

data = eu.io.read(flist)
w = np.where(data['flags']==0)
data = data[w]
times = []
j = 0
for i in range(100):
    times.append(data['time'][j])
    j += 200

flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run177/run177_3-*.fits')

data = eu.io.read(flist)
w = np.where(data['flags']==0)
data = data[w]
times2 = []
j = 0
for i in range(100):
    times2.append(data['time'][j])
    j += 200


flist = glob('/gpfs01/astro/workarea/lmezini/scarlet-tests/run177/run177_5-*.fits')

data = eu.io.read(flist)
w = np.where(data['flags']==0)
data = data[w]
times3 = []
j = 0
for i in range(98):
    times3.append(data['time'][j])
    j += 200
print("scarlet",np.sum(np.array(times))/len(flist),np.array(times).std()/10.)
print("mof",np.sum(np.array(times2))/len(flist),np.array(times2).std()/10.)
print("minimof",np.sum(np.array(times3))/98,np.array(times3).std()/np.sqrt(100))

plt.grid(axis='x',linestyle=':')
plt.hist(times,50,histtype='step',color = 'r',label="scarlet")
plt.hist(times2,50,histtype='step',color = 'b',label="mof")
plt.hist(times3,50,histtype='step',color = 'm',label="minimof")
plt.title("Runtimes(s) N=0.1")
plt.xticks(range(0,26,2))
plt.legend()
plt.savefig('runtimes_N0.1.png')
