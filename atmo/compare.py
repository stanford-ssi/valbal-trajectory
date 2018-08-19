import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys

latlons = np.load('../proc/latlons.npy')
output = np.fromfile('../ignored/output.%d.bin' % int(sys.argv[1]), dtype=np.float32)
output.shape = (len(output)//2, 2)

m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=40,\
            llcrnrlon=-130,urcrnrlon=-90,resolution='l')
xpath,ypath = m(latlons[:,1], latlons[:,0])
xpred,ypred = m(output[:,1]-360, output[:,0])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])


plt.plot(xpath,ypath)
plt.plot(xpred,ypred)
#plt.plot(output[:,1])
plt.show()
