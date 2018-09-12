import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys
import pandas as pd 
import datetime
#For plotting Monte Carlo's outs 

df = pd.read_hdf('../ignored/flights/ssi71_alt.h5') 
df2 = pd.read_hdf('../ignored/flights/ssi71_location.h5')
def load_file(i):
    output = np.fromfile('../ignored/sim/output.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//3, 3)
    return output

files = list(map(load_file, range(1,2000)))

t0 = df.index[0]
times = t0+(200+np.arange(files[0].shape[0]))*datetime.timedelta(seconds=60*10)
#df.reindex(times)
#df2.reindex(times)
plt.subplot(2,1,1)
for f in files[1:]:
    plt.plot(times,f[:,2],c="blue",alpha=0.1)
plt.plot(df.altitude_barometer.index, df.altitude_barometer.values, color='red')# df.altitude_barometer.plot(c="red") 
plt.grid()
plt.xlabel("time")
plt.ylabel("altitude")
plt.subplot(2,1,2)

m = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=60,\
            llcrnrlon=-180,urcrnrlon=10,resolution='i')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

for i,f in enumerate(files[1:]):
	xpred,ypred = m(f[:,1]-360, f[:,0])
	plt.plot(xpred,ypred,c="blue",alpha=0.3)
xpred,ypred = m(files[0][:,1]-360, files[0][:,0])
plt.plot(xpred,ypred,c="red",alpha=1)
xreal,yreal = m(df2.long_gps.values,df2.lat_gps.values)
plt.plot(xreal,yreal,c="red")
plt.show()
