import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys
import pandas as pd 
import datetime

def load_file(i):
    output = np.fromfile('../ignored/sim/output.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//3, 3)
    return output
files = list(map(load_file, range(1,1000)))
times = np.arange(files[0].shape[0])*10/60
plt.subplot(2,1,1)
for file in files:
	plt.plot(times,file[:,2],color='blue',alpha=0.1)

plt.plot(times,files[0][:,2],color='red')
plt.plot(times,files[-1][:,2],color='green')

plt.subplot(2,1,2)
m = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=70,
            llcrnrlon=-180,urcrnrlon=80,resolution='c')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

for i,f in enumerate(files[1:]):
	xpred,ypred = m(f[:,1]-360, f[:,0])
	plt.plot(xpred,ypred,c="blue",alpha=0.1)
xpred,ypred = m(files[0][:,1]-360, files[0][:,0])
plt.plot(xpred,ypred,c="red")
xpred,ypred = m(files[-1][:,1]-360, files[-1][:,0])
plt.plot(xpred,ypred,c="green")

xt,yt = m( 13.405491, 52.516655)
plt.plot(xt,yt,"r*")
plt.show()


