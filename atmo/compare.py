import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys

def load_file(i):
    latlons = np.load('../proc/latlons.npy')
    output = np.fromfile('../ignored/output.%s.bin' % i, dtype=np.float32)
    output.shape = (len(output)//3, 3)
    return output[:,0], output[:,1], output[:,2]

files = list(map(load_file, sys.argv[1:]))
for f in files:
    plt.plot(f[2])
plt.figure()

m = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=60,\
            llcrnrlon=-130,urcrnrlon=-50,resolution='l')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

for f in files:
    xpred,ypred = m(f[1]-360, f[0])
    plt.plot(xpred,ypred)

plt.show()
