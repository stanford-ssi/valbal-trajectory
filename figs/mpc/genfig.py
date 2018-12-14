import numpy as np 
import pandas as pd 
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib import gridspec
from datetime import datetime  
from datetime import timedelta
import matplotlib.lines as li
def load_file2(i,name="output"):
    output = np.fromfile('data2/'+name+'.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//5, 5)
    return output

def f2():
	#was used for plotting stuff from the evaluator
	clist = ["green","red","purple","orange"]
	ptimes = [0,24,48,72]
	target = np.array([[42.733156, 13.449431+360]]);
	files = list(map(load_file2, range(1,79)))
	const = list(map(load_file2, [0]))
	files2 = []
	for i in ptimes:
		files2.append(list(map(load_file2, range(100+i*100,100+(i+1)*100))))

	fig, ax = plt.subplots(2,3,gridspec_kw = {'height_ratios':[1,.63], "hspace":0,},figsize=(10,5))
	gs = gridspec.GridSpec(2,3, width_ratios=[0.06, .14, 1],height_ratios=[1, 0.56])
	gs.update(left=0.05,wspace=0,hspace=0.05)
	ax0 = plt.subplot(gs[0,2:3])
	ax1 = plt.subplot(gs[1,1:3])
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
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 + real_vals[:,4]/1000,"--",c="blue",alpha=1)
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 - real_vals[:,4]/1000,"--",c="blue",alpha=1)
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
	ax0.plot(xpred,ypred,color="black",alpha=1,label="pre-optimized \ntrajectory")
	ax0.plot(xpred[-1],ypred[-1],"*",c="black")	

	xpred,ypred = m(real_vals[:,1], real_vals[:,0])
	ax0.plot(xpred,ypred,color="blue",alpha=1,label="evalaton \nsimulation \ntrajectory")
	ax0.plot(xpred[-1],ypred[-1],"*",c="blue")	

	xt,yt = m(target[0,1],target[0,0])
	p = ax0.plot(xt,yt,".",c="red")


	legend_elements = [
	li.Line2D([0], [0], color='b', lw=1.5, label="evalaton \nsimulation \ntrajectory"),
	li.Line2D([0], [0], color='black', lw=1.5, label="pre-optimized \ntrajectory"),
	li.Line2D([0], [0], linestyle="--",color='b', alpha=0, lw=2, label="------------------"),
	li.Line2D([0], [0], color='b', lw=1.5, label="evalaton \nsimulation \naltitude"),
	li.Line2D([0], [0], linestyle="--",color='b', lw=1.5, label="controller \ncounds"),
	li.Line2D([0], [0], marker='o',markerfacecolor='r',color='w', label='goal')
	]
	print(p)
	ax0.legend(handles=legend_elements,loc='uppper right', bbox_to_anchor=(0, 1.03),ncol=1)
	plt.savefig("mpc2.png")

def load_file1(i,name="output"):
    output = np.fromfile('data1/'+name+'.%03d.bin' % int(i), dtype=np.float32)
    output.shape = (len(output)//5, 5)
    return output

def f1():
	#was used for plotting stuff from the evaluator
	clist = ["green","red","purple","orange"]
	ptimes = [0,24,48,72]
	target = np.array([[19.805777, 13.702091+360]]);
	files = list(map(load_file1, range(1,89)))
	const = list(map(load_file1, [0]))
	files2 = []
	for i in ptimes:
		files2.append(list(map(load_file1, range(100+i*100,100+(i+1)*100))))

	fig, ax = plt.subplots(2,3,gridspec_kw = {'height_ratios':[1,.63], "hspace":0,},figsize=(10,5))
	gs = gridspec.GridSpec(2,3, width_ratios=[0.06, .14, 1],height_ratios=[1, 0.56])
	gs.update(left=0.05,wspace=0,hspace=0.05)
	ax0 = plt.subplot(gs[0,2:3])
	ax1 = plt.subplot(gs[1,1:3])
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
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 + real_vals[:,4]/1000,"--",c="blue",alpha=1)
	ax1.plot(np.arange(N)/6,real_vals[:,3]/1000 - real_vals[:,4]/1000,"--",c="blue",alpha=1)
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
	ax0.plot(xpred,ypred,color="black",alpha=1,label="pre-optimized \ntrajectory")
	ax0.plot(xpred[-1],ypred[-1],"*",c="black")	

	xpred,ypred = m(real_vals[:,1], real_vals[:,0])
	ax0.plot(xpred,ypred,color="blue",alpha=1,label="evalaton \nsimulation \ntrajectory")
	ax0.plot(xpred[-1],ypred[-1],"*",c="blue")	

	xt,yt = m(target[0,1],target[0,0])
	p = ax0.plot(xt,yt,".",c="red")


	legend_elements = [
	li.Line2D([0], [0], color='b', lw=1.5, label="evalaton \nsimulation \ntrajectory"),
	li.Line2D([0], [0], color='black', lw=1.5, label="pre-optimized \ntrajectory"),
	li.Line2D([0], [0], linestyle="--",color='b', alpha=0, lw=2, label="------------------"),
	li.Line2D([0], [0], color='b', lw=1.5, label="evalaton \nsimulation \naltitude"),
	li.Line2D([0], [0], linestyle="--",color='b', lw=1.5, label="controller \ncounds"),
	li.Line2D([0], [0], marker='o',markerfacecolor='r',color='w', label='goal')
	]
	print(p)
	ax0.legend(handles=legend_elements,loc='uppper right', bbox_to_anchor=(0, 1.03),ncol=1)
	plt.savefig("mpc1.png")

f2()