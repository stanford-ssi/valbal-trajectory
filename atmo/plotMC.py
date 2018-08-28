import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys

#For plotting Monte Carlo's outs 


def load_file(i):
    output = np.fromfile('../ignored/output.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//3, 3)
    return output

files = list(map(load_file, range(0,10)))
plt.subplot(2,1,1)
for f in files:
    plt.plot(f[:,2],c="blue",alpha=0.3)
plt.grid()
plt.xlabel("time")
plt.ylabel("altitude")
plt.subplot(2,1,2)

m = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=60,\
            llcrnrlon=-130,urcrnrlon=-50,resolution='i')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

for f in files:
    xpred,ypred = m(f[:,1]-360, f[:,0])
    plt.plot(xpred,ypred,c="blue",alpha=0.3)

plt.show()
