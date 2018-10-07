import numpy as np 
import atmotools as at 
import pandas as pd 
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 
import matplotlib


def plotWindOverMap(region,alts,times,db='gfs_anl_0deg5',traj=None):
	""" Region is a list of format [lon min, lon max, lat min, lat max]
		alts is a list of format [alt min, alt max]
		times a list where of time strings of format "YYYY-MM-DD_HH" or 
		db is the database used
	"""
	N = len(times)
	fig, axes = plt.subplots(nrows=N, ncols=1)
	if N ==1:
		axes = [axes]
	for t in range(N):
		ret = at.procWindData(times[t],times[t])
		time = ret[1][0]
		file = ret[0][0]
		windobj = at.getKeysForBin(file)
		data = at.getArrayFromBin(file,windobj)

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
		for i in range(altinds.size):
			U = data[latinds,loninds,altinds[i],0]
			V = data[latinds,loninds,altinds[i],1]
			axes[t].quiver(x,y,U,V,units="xy",color=cmap((i / (altinds.size-1))*0.8 + 0.2),headlength=1, headwidth=1,width=10000,pivot='mid',alpha=0.7)
		if traj:
			x,y = m(traj[0,:],traj[1,:])
			axes[t].plot(x,y)

		axes[t].set_title(time)
	fig.subplots_adjust(right=0.85)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
	fig.colorbar(sm, ticks=alts_, cax=cbar_ax)
	return axes,fig,m
	#plt.show()
'''
np.set_printoptions(edgeitems=1000,linewidth=1000)
ret = at.procWindData("2018-09-09_00","2018-09-09_00",overwrite=True)
time = ret[1][0]
file = ret[0][0]
windobj = at.getKeysForBin(file)
data = at.getArrayFromBin(file,windobj)
#plotWindOverMap([-100,-40,0,50],[1,20],["2018-09-10_00"])
'''