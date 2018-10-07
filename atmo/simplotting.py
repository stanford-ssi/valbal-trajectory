import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys
import pandas as pd 
import datetime
from windvistools import *
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
	files = list(map(load_file, range(500)))
	plt.subplot(2,1,1)
	for i,f in enumerate(files):
		plt.plot(f[:,2],c="blue",alpha=0.01)
	plt.plot(files[0][:,2],c="red")	

	plt.subplot(2,1,2)
	print(files[0])
	m = Basemap(projection='merc',llcrnrlat=np.min(files[0][:,0])-5,urcrnrlat=np.max(files[0][:,0])+5,
            llcrnrlon=np.min(files[0][:,1])-5,urcrnrlon=np.max(files[0][:,1])+5,resolution='l')
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	for i,f in enumerate(files):
		xpred,ypred = m(f[:,1], f[:,0])
		plt.plot(xpred,ypred,color="blue",alpha=0.05)
		plt.plot(xpred[-1],ypred[-1],"r*")
	plt.show()

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