import numpy as np 
import atmotools as at 
import pandas as pd 
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 
import matplotlib


def plotWindOverMap(region,alts,times,windobj,traj=None):
	""" Region is a list of format [lon min, lon max, lat min, lat max]
		alts is a list of format [alt min, alt max]
		Windobj is an dict of wind data of the format output by makeWindArray in atmotools
	"""
	N = len(times)
	fig, axes = plt.subplots(nrows=N, ncols=1)
	hrs = windobj["times"]
	gti = lambda t : (int(round((hrs[0] - t)/(hrs[0] - hrs[-1])*(hrs.size-1))),int(round((hrs[0] - t)/(hrs[0] - hrs[-1])*(hrs.size-1)*2)) % 2)
	for t in range(N):
		m = Basemap(projection='merc',llcrnrlat=region[2],urcrnrlat=region[3],
	            llcrnrlon=region[0],urcrnrlon=region[1],resolution='l',ax=axes[t])
		loninds = np.where(np.logical_and(region[0]<=(windobj["lons"]+180)%360-180,(windobj["lons"]+180)%360-180<=region[1]))[0][:,np.newaxis]
		latinds = np.where(np.logical_and(region[2]<=windobj["lats"],windobj["lats"]<=region[3]))[0][np.newaxis,:]
		altinds = np.where(np.logical_and(alts[0]*1000<=windobj["alts"],windobj["alts"]<=1000*alts[1]))[0]
		latgr,longr = np.meshgrid(windobj["lats"][latinds],(windobj["lons"][loninds]+180)%360-180)
		x,y = m(longr,latgr)
		m.drawcoastlines()
		m.drawcountries()
		m.drawstates()

		alts_ = windobj['alts'][altinds]
		cmap = matplotlib.cm.get_cmap('jet')
		norm = matplotlib.colors.Normalize()
		norm.autoscale(alts_)
		sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
		sm.set_array([])
		ti1,ti2 = gti(times[t])
		for i in range(altinds.size):
			U = windobj["data"][loninds,latinds,altinds[i],ti1,ti2,0]*1e7
			V = windobj["data"][loninds,latinds,altinds[i],ti1,ti2,1]*1e7
			axes[t].quiver(x,y,U,V,units="xy",color=cmap((i / (altinds.size-1))*0.8 + 0.2),headlength=1, headwidth=1,width=10000,pivot='mid',alpha=0.7)
		if traj:
			x,y = m(traj[0,:],traj[1,:])
			axes[t].plot(x,y)

		axes[t].set_title("%04d-%02d-%02d %02d hr + %02d hrs"%(windobj["start_time"].year,windobj["start_time"].month,windobj["start_time"].day,windobj["start_time"].hour,times[t]))
	fig.subplots_adjust(right=0.85)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
	fig.colorbar(sm, ticks=alts_, cax=cbar_ax)
	plt.show()




windobj = "../ignored/GFS_anl_0deg5_objs/2017-12-09_18_to_2017-12-15_05.pickle"
winddata = "../ignored/GFS_anl_0deg5_objs/2017-12-09_18_to_2017-12-15_05.npy"
data = np.load(winddata)
windobj = pickle.load(open(windobj,'rb'))
windobj["data"] = data

plotWindOverMap([-130,-100,30,45],[10,16],[0,24],windobj)