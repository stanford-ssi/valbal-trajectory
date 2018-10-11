import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import pandas as pd 
import datetime
from windvistools import *
import itertools as it
def load_file(i):
    output = np.fromfile('../ignored/sim/output.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//3, 3)
    return output

def plot1():
	# jank thing for plotting flights on top of a wind vectory field
	files = list(map(load_file, range(1,1000)))
	times = np.arange(files[0].shape[0])*10/60
	for file in files:
		plt.plot(times,file[:,2],color='blue',alpha=0.1)

	plt.plot(times,files[0][:,2],color='red')
	plt.plot(times,files[-1][:,2],color='green')

	axes,fig,m = plotWindOverMap([-100,-40,10,50],[1,20],["2018-09-13_00"])

	for i,f in enumerate(files[1:]):
		xpred,ypred = m(f[:,1]-360, f[:,0])
		axes[0].plot(xpred,ypred,c="blue",alpha=0.1)
	xpred,ypred = m(files[0][:,1]-360, files[0][:,0])
	axes[0].plot(xpred,ypred,c="red")
	xpred,ypred = m(files[-1][:,1]-360, files[-1][:,0])
	axes[0].plot(xpred,ypred,c="green")
	axes[0].plot(xpred[::6*24],ypred[::6*24],'g*')

	ll = [59.916193, 30.325234]
	xt,yt = m(ll[1],ll[0])
	axes[0].plot(xt,yt,"r*")
	plt.show()

def plotruns():
	files = list(map(load_file, range(50)))
	files2 = list(map(load_file, range(50,100)))
	fig, ax = plt.subplots(2,2,gridspec_kw = {'height_ratios':[1,.63], "hspace":0,},figsize=(8,5))
	gs = gridspec.GridSpec(2, 2, width_ratios=[0.1, 1],height_ratios=[1, 0.56])
	gs.update(left=0.05,wspace=0,hspace=0.05)
	ax0 = plt.subplot(gs[0,0:2])
	ax1 = plt.subplot(gs[1,1])
	ax1.plot(0,14,c="#be1e2d")
	ax1.plot(0,14,c="blue")
	for i,f in enumerate(files):
		f2 = files2[i]
		N = f.shape[0]
		N2 = f2.shape[0]
		ax1.plot(np.arange(N2)/6,f2[:,2]/1000,c="#be1e2d",alpha=0.05)
		ax1.plot(np.arange(N)/6,f[:,2]/1000,c="blue",alpha=0.05)
	ax1.legend(["initial","optimized"],loc=(.18,.05))
	ax1.plot(np.arange(files[0].shape[0])/6,files[0][:,2]/1000,c="blue")	
	ax1.set_xlabel("flight time (hr)")
	ax1.set_ylabel("altitude (km)")
	ax1.grid()
	all_vals = np.concatenate(files)

	m = Basemap(projection='merc',llcrnrlat=np.min(all_vals[:,0])-5,urcrnrlat=np.max(all_vals[:,0])+5,
            llcrnrlon=np.min(all_vals[:,1])-5,urcrnrlon=np.max(all_vals[:,1])+5,resolution='l',ax=ax0,)
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")
	for i,f in enumerate(files):
		if (i+1)%3 != 0:
			continue
		f2 = files2[i]
		xpred,ypred = m(f[:,1], f[:,0])
		x2,y2 = m(f2[:,1], f2[:,0])
		ax0.plot(x2,y2,color="#be1e2d",alpha=0.3)
		ax0.plot(x2[-1],y2[-1],"*",c="#be1e2d")
		ax0.plot(xpred,ypred,color="blue",alpha=0.1)
		ax0.plot(xpred[-1],ypred[-1],"*",c="blue")
	plt.savefig("../ignored/figs/plot1.png")

def plotflights():
	#used for plotting and comparing with real flights
	files = list(map(load_file, range(1,200)))
	df = pd.read_hdf('../ignored/flights/ssi67_location.h5')
	print(df.keys())
	m = Basemap(projection='merc',llcrnrlat=np.min(df.lat_gps.values) -5,urcrnrlat=np.max(df.lat_gps.values)+5,
            llcrnrlon=np.min(df.long_gps.values)-5,urcrnrlon=np.max(df.long_gps.values)+5,resolution='l')
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	for i,f in enumerate(files[1:]):
		xpred,ypred = m(f[:,1]-360, f[:,0])
		plt.plot(xpred,ypred,c="blue",alpha=0.5)
	x,y = m(df.long_gps.values,df.lat_gps.values)
	plt.plot(x,y,c="red")
	plt.show()



plotruns()