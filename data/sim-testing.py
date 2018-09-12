import numpy as np 
import atmotools as at 
import pandas as pd 
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 

if 0:
	df = pd.read_hdf('../ignored/ssi63.h5')
	[print(x) for x in df.keys()]
	dfs = df[['raw_pressure_1','raw_pressure_2','raw_pressure_3','raw_pressure_4','lat_gps','long_gps']]
	dfs.to_hdf('ssi63_position.h5','df',complib='zlib', mode='w', complevel=5)
	exit()
if 1:
	df = pd.read_hdf('ssi63_position.h5')
	at.procWindData(df.index[0],df.index[-1],db="gfs_anl_0deg5")
	#windkey,winddata = at.makeWindArray(df.index[0],df.index[-1],overwrite=True)
	exit()

windobj = "../ignored/GFS_anl_0deg5_objs/2017-12-09_18_to_2017-12-15_05.pickle"
winddata = "../ignored/GFS_anl_0deg5_objs/2017-12-09_18_to_2017-12-15_05.npy"
#windobj["data"] = winddata

df = pd.read_hdf('ssi63_position.h5')
data = np.load(winddata)
key = pickle.load(open(windobj,'rb'))

fhrs = (df.index[::20*60*10] - key["start_time"]).days*24 + (df.index[::20*60*10] - key["start_time"]).seconds/3600 
fhrs = fhrs.values
fl = ((df.raw_pressure_1.values + df.raw_pressure_2.values + df.raw_pressure_3.values + df.raw_pressure_4.values)/4/100)[::20*60*10]
flons = df.long_gps.values[::20*60*10]
flats = df.lat_gps.values[::20*60*10]

print(key["alts"])

lats = key["lats"]
lons = key["lons"]
hrs = key["times"]
levels = key["levels"]
glati = lambda lat : int(round((lats[0] - lat)/(lats[0] - lats[-1])*(lats.size-1)))
gloni = lambda lon : int(round((lons[0] - lon)/(lons[0] - lons[-1])*(lons.size-1)))
gti = lambda t : (int(round((hrs[0] - t)/(hrs[0] - hrs[-1])*(hrs.size-1))),int(round((hrs[0] - t)/(hrs[0] - hrs[-1])*(hrs.size-1)*2)) % 2)
gli = lambda l : np.argmin(np.abs(levels - l))
#yolo, 111111 meter/deg ftw

starthr = 40
startind = np.argmin(np.abs(fhrs - starthr))

idx = startind
lon = flons[idx]
lat = flats[idx]
N = fhrs.size - startind
ploc = np.ones((2,N))

for i in range(N):
	vels = data[gloni(lon),glati(lat),gli(fl[idx]),gti(fhrs[idx])[0],gti(fhrs[idx])[1],:]
	print(np.cos(lat/180*np.pi),lat)
	lon += vels[0]*10*60/111111/np.cos(lat/180*np.pi)
	lat += vels[1]*10*60/111111
	ploc[:,i] = [lon,lat]
	idx+=1


lat = key["lats"][90-50:90-25] 
lon = key["lons"][360-125:360-70]
latgr,longr = np.meshgrid(lat,lon)
time = 4*1
m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=40,\
            llcrnrlon=-130,urcrnrlon=-90,resolution='l')
x,y = m(longr-360,latgr)
xpath,ypath = m(df.long_gps.values,df.lat_gps.values)
xpred,ypred = m(ploc[0,:],ploc[1,:])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,81,10.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

plt.plot(xpath,ypath)
plt.plot(xpred,ypred)
plt.show()


