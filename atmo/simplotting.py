import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import pandas as pd 
import datetime
from windvistools import *
import itertools as it
#import gmplot
import os
import time
from datetime import datetime  
from datetime import timedelta
def load_file(i,name="output"):
    output = np.fromfile('../ignored/sim/'+name+'.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//5, 5)
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
	#plt.savefig("fig.png")

def plot2():
	# This function generates a nice plot of optimized montecarlo trajectories.
	# It was used to make the fig that's in the vb data sheet on Oct. 12th 2018.
	fnum = 3
	files = list(map(load_file, range(100*fnum,50+100*fnum)))
	files2 = list(map(load_file, range(50+100*fnum,100+100*fnum)))
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
	ax1.legend(["unoptimized","optimized"])#,loc=(.18,.05))

	ax1.plot(np.arange(files[0].shape[0])/6,files[0][:,2]/1000,c="blue")	
	#ax1.set_xlabel("hours from " + datetime.fromtimestamp(1541045606).strftime("%Y-%m-%d %H:%M:%S"))
	ax1.set_xlabel("hours from launch")
	ax1.set_ylabel("altitude (km)")
	ax1.grid()
	all_vals = np.concatenate(files + files2)

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
		#ax0.plot(x2[::6*10],y2[::6*10],"*",c="#be1e2d",alpha=0.3)
		ax0.plot(xpred,ypred,color="blue",alpha=0.1)
		ax0.plot(xpred[-1],ypred[-1],"*",c="blue")
		#ax0.plot(xpred[::6*10],ypred[::6*10],"*",c="blue",alpha=0.3)
	#plt.savefig("../ignored/figs/plot1.png")
	plt.show()
	#plt.savefig("fig.png")


def plot3():
	# this plot generates trajectories from CE gradient descent
	N_iter = 700
	N_runs = 100
	fig, ax = plt.subplots(3,1,figsize=(8,5))
	gs = gridspec.GridSpec(3,1, width_ratios=[1],height_ratios=[1, 0.56,.56])
	gs.update(left=0.1,wspace=0.1,hspace=0.2)
	ax0 = plt.subplot(gs[0,0])
	ax1 = plt.subplot(gs[1,0])
	ax2 = plt.subplot(gs[2,0])
	files = list(map(load_file, range(0,N_iter*N_runs)))
	all_vals = np.concatenate(files)
	m = Basemap(projection='merc',llcrnrlat=np.min(all_vals[:,0])-5,urcrnrlat=np.max(all_vals[:,0])+5,
	            llcrnrlon=np.min(all_vals[:,1])-5,urcrnrlon=np.max(all_vals[:,1])+5,resolution='l',ax=ax0,)
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")
	ax1.set_xlabel("flight time (hr)")
	ax1.set_ylabel("altitude (km)")
	ax1.grid()
	ax2.set_ylabel("objective value")
	ax2.set_ylabel("iteration")
	for j in range(N_runs):
		print(j)
		files = list(map(load_file, range(N_iter*j+699,N_iter*(j+1))))
		obj = np.fromfile('../ignored/sim/opt.%03d.bin' % j, dtype=np.float32)
		for i,f in enumerate(files):
			N = f.shape[0]
			ax1.plot(np.arange(N)/6,f[:,2]/1000,c="blue",alpha=0.05)
		f = files[-1]	
		ax1.plot(np.arange(N)/6,f[:,2]/1000,c="blue",alpha=0.1)

		for i,f in enumerate(files):
			xpred,ypred = m(np.mod(f[:,1]+180,360)-180, f[:,0])
			ax0.plot(xpred,ypred,color="blue",alpha=0.1)
			ax0.plot(xpred[-1],ypred[-1],"*",c="green")
		f = files[-1]
		xpred,ypred = m(f[:,1], f[:,0])
		ax0.plot(xpred,ypred,color="blue",alpha=0.1)
		ax0.plot(xpred[-1],ypred[-1],"*",c="red")
		ax2.plot(obj)
		#plt.savefig("../ignored/figs/plot1.png")
	plt.show()

def plot4():
	# this function actually has nothing to do with simulation, but rather 
	# gerates a nice figure of ssi63 for the valbal datasheet
	data = np.load('../ignored/flights/ssi63_latlons.npy')
	lat=36.845679
	lon=-121.402538
	m = Basemap(projection='merc',llcrnrlat=lat-0.3,urcrnrlat=lat+.7,
            llcrnrlon=lon-1,urcrnrlon=lon+1,resolution='i')
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")
	m.shadedrelief()
	x,y =m(data[:,1],data[:,0])
	plt.plot(x,y)
	plt.plot(x[0],y[0],"*")
	plt.show()


def plot5():
	# this plot generates trajectories from CE gradient descent
	N_iter = 700
	N_runs = 100
	runnums = np.arange(16,18)
	fnums = runnums*N_iter - 1
	if 0:
		cmd = 'scp john@$D:~/SSI/valbal/valbal-trajectory/ignored/sim/\{'
		for j in runnums:
			cmd += "opt.%03d.bin," % j
		for i in fnums:
			cmd += "output.%03d.bin," % i
		cmd = cmd[:-1]
		cmd +=  "\} ../ignored/sim/"
		print(cmd)
		os.system(cmd)
	fig, ax = plt.subplots(3,1,figsize=(8,5))
	gs = gridspec.GridSpec(3,1, width_ratios=[1],height_ratios=[1, 0.56,.56])
	gs.update(left=0.1,wspace=0.1,hspace=0.2)
	ax0 = plt.subplot(gs[0,0])
	ax1 = plt.subplot(gs[1,0])
	ax2 = plt.subplot(gs[2,0])
	files = list(map(load_file, fnums))
	all_vals = np.concatenate(files)
	m = Basemap(projection='merc',llcrnrlat=np.min(all_vals[:,0])-5,urcrnrlat=np.max(all_vals[:,0])+5,
	            llcrnrlon=np.min(all_vals[:,1])-5,urcrnrlon=np.max(all_vals[:,1])+5,resolution='l',ax=ax0,)
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")
	ax1.set_xlabel("flight time (hr)")
	ax1.set_ylabel("altitude (km)")
	ax1.grid()
	ax2.set_ylabel("objective value")
	ax2.set_ylabel("iteration")
	for j in runnums:
		print(j)
		files = list(map(load_file,[j*N_iter - 1]))
		obj = np.fromfile('../ignored/sim/opt.%03d.bin' % j, dtype=np.float32)
		for i,f in enumerate(files):
			N = f.shape[0]
			ax1.plot(np.arange(N)/6,f[:,2]/1000,alpha=1)
		f = files[-1]	
		ax1.plot(np.arange(N)/6,f[:,2]/1000,alpha=1)

		for i,f in enumerate(files):
			xpred,ypred = m(np.mod(f[:,1]+180,360)-180, f[:,0])
			ax0.plot(xpred,ypred,alpha=1)
			ax0.plot(xpred[-1],ypred[-1],"*",c="green")
		f = files[-1]
		xpred,ypred = m(f[:,1], f[:,0])
		ax0.plot(xpred,ypred,alpha=1)
		ax0.plot(xpred[-1],ypred[-1],"*",c="red")
		ax2.plot(obj)
		#plt.savefig("../ignored/figs/plot1.png")
	plt.show()

def plotssi73():
	## for messing with SSI73
	files = list(map(load_file, range(1000)))
	f=files[0]
	t = (np.arange(f.shape[0])*30 + 1612)/60/60
	locs = np.where(np.logical_and(t>(1+2.5/6),t<1+3/6))[0]
	print(locs)
	gmap3 = gmplot.GoogleMapPlotter(f[0,0],f[0,1], 10) 
	for f in files:
		#gmap3.scatter((f[-1,0],f[-1,0]),(f[-1,1],f[-1,0]),'#FF0000',size = 500, marker = False ) 
		gmap3.plot(f[locs,0],f[locs,1],'# FF0000',size = 500, marker = False ) 
	gmap3.draw("map.html")
	exit()
	plt.subplot(211)
	for f in files:
		plt.plot(t,f[:,2],c="blue",alpha=0.1)
	plt.title("t=0 is at unix time 1540060749")
	plt.xlabel("time (hr)")
	plt.ylabel("alt (m)")
	plt.grid()
	plt.subplot(212)
	all_vals = np.concatenate(files)
	s = 2
	m = Basemap(projection='merc',llcrnrlat=np.min(all_vals[:,0])-s,urcrnrlat=np.max(all_vals[:,0])+s,
            llcrnrlon=np.min(all_vals[:,1])-s,urcrnrlon=np.max(all_vals[:,1])+s,resolution='l')
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	for f in files:
		xpred,ypred = m(f[:,1], f[:,0])
		plt.plot(xpred,ypred,c="blue",alpha=0.05)
		plt.plot(xpred[-1],ypred[-1],"*",c="green",alpha=0.1)
		#plt.plot(f[:,1],f[:,0],c="blue",alpha=0.03)
		#plt.plot(f[-1,1],f[-1,0],"g*")
		#plt.plot(f[0,1],f[0,0],"r*")
	#plt.title("red dot is start")
	plt.grid()
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

def plotmc():
	#used for plotting and comparing with real flights
	files = list(map(load_file, range(1,2000)))
	minlat = files[0][:,0].min()
	maxlat = files[0][:,0].max()
	minlon = files[0][:,1].min()
	maxlon = files[0][:,1].max()
	print(minlat,maxlat, minlon, maxlon)
	m = Basemap(projection='merc',llcrnrlat=minlat-5,urcrnrlat=maxlat+5,
            llcrnrlon=minlon-5,urcrnrlon=maxlon+5,resolution='l')
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	for i,f in enumerate(files[1:]):
		xpred,ypred = m(f[:,1], f[:,0])
		plt.plot(xpred,ypred,c="blue",alpha=0.5)
	plt.show()

def plot6():
	#was used for plotting stuff from the evaluator
	clist = ["green","red","purple","orange"]
	ptimes = [0,24,48,72]
	target = np.array([[19.805777, 13.702091+360]]);
	files = list(map(load_file, range(1,87)))
	const = list(map(load_file, [0]))
	files2 = []
	for i in ptimes:
		files2.append(list(map(load_file, range(100+i*100,100+(i+1)*100))))

	fig, ax = plt.subplots(2,2,gridspec_kw = {'height_ratios':[1,.63], "hspace":0,},figsize=(8,5))
	gs = gridspec.GridSpec(2, 2, width_ratios=[0.1, 1],height_ratios=[1, 0.56])
	gs.update(left=0.05,wspace=0,hspace=0.05)
	ax0 = plt.subplot(gs[0,0:2])
	ax1 = plt.subplot(gs[1,1])
	ax1.plot(0,14,c="#be1e2d")
	ax1.plot(0,14,c="blue")

	all_vals = np.concatenate(files + [i1 for i2 in files2 for i1 in i2])
	real_vals = np.concatenate(files)

	for j,fs in enumerate(files2):
		c = clist[j]
		vals = fs[0]
		N = vals.shape[0]
		ax1.plot(ptimes[j] + np.arange(N)/6,vals[:,3]/1000 - vals[:,4]/1000,"--",c=c,alpha=0.5)
		ax1.plot(ptimes[j] + np.arange(N)/6,vals[:,3]/1000 + vals[:,4]/1000,"--",c=c,alpha=0.5)
		for i in range(0,50):
			if (i+1)%8 != 0:
				continue
			vals = fs[i];
			N = vals.shape[0]
			#ax1.plot(ptimes[j] + np.arange(N)/6,vals[:,2]/1000,c=c,alpha=0.1)

	N = real_vals.shape[0]
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 + real_vals[:,4]/1000,"--",c="blue",alpha=0.5)
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 - real_vals[:,4]/1000,"--",c="blue",alpha=0.5)
	ax1.plot(np.arange(N)/6,real_vals[:,2]/1000,c="blue",alpha=1)
	#ax1.legend(["initial","optimized"],loc=(.18,.05))
	ax1.set_xlabel("hours from " + datetime.fromtimestamp(1543492801).strftime("%Y-%m-%d %H:%M:%S"))
	ax1.set_ylabel("altitude (km)")
	ax1.grid()

	all_lls = np.vstack((all_vals[:,:2],target));
	m = Basemap(projection='merc',llcrnrlat=np.min(all_lls[:,0])-5,urcrnrlat=np.max(all_lls[:,0])+5,
            llcrnrlon=np.min(all_lls[:,1])-5,urcrnrlon=np.max(all_lls[:,1])+5,resolution='l',ax=ax0,)
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")

	for j,fs in enumerate(files2):
		for i in range(0,50):
			if (i+1)%8 != 0:
				continue
			vals = fs[i];
			xpred,ypred = m(vals[:,1], vals[:,0])
			c = clist[j]
			ax0.plot(xpred,ypred,color=c,alpha=0.2)
			#ax0.plot(xpred[-1],ypred[-1],"*",c=c,alpha=0.2)
		ax0.plot(xpred[0],ypred[0],"*",c=c,alpha=1)	

	cv = const[0]
	xpred,ypred = m(cv[:,1], cv[:,0])
	ax0.plot(xpred,ypred,color="black",alpha=1)
	ax0.plot(xpred[-1],ypred[-1],"*",c="black")	

	xpred,ypred = m(real_vals[:,1], real_vals[:,0])
	ax0.plot(xpred,ypred,color="blue",alpha=1)
	ax0.plot(xpred[-1],ypred[-1],"*",c="blue")	

	xt,yt = m(target[0,1],target[0,0])
	ax0.plot(xt,yt,"o",c="red")
	plt.show()

def plot7():
	# used for plotting multiple random starts
	fig, ax = plt.subplots(2,2,gridspec_kw = {'height_ratios':[1,.63], "hspace":0,},figsize=(8,5))
	gs = gridspec.GridSpec(2, 2, width_ratios=[0.1, 1],height_ratios=[1, 0.56])
	gs.update(left=0.05,wspace=0,hspace=0.05)
	ax0 = plt.subplot(gs[0,0:2])
	ax1 = plt.subplot(gs[1,1])
	nstart = 5;
	fnums = it.chain(*[[100*j + i for i in range(50)] for j in range(nstart)])
	fnums2 = it.chain(*[[100*j + 50 + i for i in range(50)] for j in range(nstart)])
	endfiles = list(map(load_file, fnums))
	startfiles = list(map(load_file, fnums2))
	all_vals = np.concatenate(endfiles + startfiles)
	m = Basemap(projection='merc',llcrnrlat=np.min(all_vals[:,0])-5,urcrnrlat=np.max(all_vals[:,0])+5,
            llcrnrlon=np.min(all_vals[:,1])-5 ,urcrnrlon=np.max(all_vals[:,1])+5,resolution='l',ax=ax0,)
	m.drawcoastlines(color="grey")
	m.drawcountries(color="grey")
	m.drawstates(color="grey")

	for j in range(nstart):
		endf =  endfiles[50*j:50*(j+1)]
		startf =  startfiles[50*j:50*(j+1)]
		avgs = []
		Ns = []
		for f in [startf,endf]:
			sizes = list(map(lambda x: x.shape[0],f))
			avg = np.zeros((max(sizes),f[0].shape[1]))
			for a in f:
				a = np.pad(a,((0,avg.shape[0]-a.shape[0]),(0,0)),mode='edge')
				print("yooo",a.shape,avg.shape)
				avg += a
			avg /= len(endf)
			avgs.append(avg)
			Ns.append(avg.shape[0])		
		np.save("../ignored/misc/start.%03d"%j,avgs[0])
		np.save("../ignored/misc/end.%03d"%j,avgs[0])

		for i in [0,1]:
			c=["black","blue"][i]
			x,y = m(avgs[i][:,1], avgs[i][:,0])
			ax0.plot(x,y,alpha=0.3,color=c)
			ax1.plot(np.arange(Ns[i])/6,avgs[i][:,3]/1000,c=c,alpha=0.5)
			#ax1.plot(np.arange(Ns[i])/6,avgs[i][:,3]/1000 + avgs[i][:,4]/1000,"--",c=c,alpha=0.5)
	ax1.grid()
	ax1.set_xlabel("hours from launch")
	ax1.set_ylabel("altitude (km)")
	ax1.legend(["unoptimized","optimized"])#,loc=(.18,.05))
	plt.show()
	#plt.savefig("fig.png")

def plot_opt():
	nvars = 3
	for j in range(5):
		output = np.fromfile('../ignored/sim/opt.%03d.bin' % int(j), dtype=np.float32).reshape(-1,nvars)
		for i in range(nvars):
			plt.subplot(nvars,1,i+1)
			if i==2:
				plt.plot(np.diff(output[:,i]))
			else:
				plt.plot(output[:,i])
			plt.grid()
	plt.show()


#plot_opt()

#plotmc()
plot2()

