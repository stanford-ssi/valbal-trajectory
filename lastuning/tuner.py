import subprocess
import numpy as np 
import matplotlib.pyplot as plt 
from multiprocessing import Pool
import functools
from datetime import datetime

def run(val,arg,seed):
	args = ("./lasTune",arg,str(val),"--seed",str(seed));
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)	
	popen.wait()
	output = popen.stdout.read()
	return float(output.decode("utf-8")[:-1])

def frun(arg,seed=0):
	copier = functools.partial(run,arg=arg,seed=seed)
	return copier

def tol_gen():
	tol_ = np.linspace(0,5000,30)
	balavg_ = np.zeros(tol_.size)
	N = 20
	for j in range(N):
		p = Pool(8)
		bal_ = np.array(p.map(frun("--tol"),tol_))
		#plt.plot(tol_,bal_,"g*")
		balavg_ += bal_/N
	np.savetxt("../ignored/misc/ballast_use.txt",np.hstack((balavg_,tol_)))
	plt.plot(tol_,balavg_)
	plt.grid()
	plt.ylabel("ballast use (g/hr)")
	plt.xlabel("tollerance (m)")	
	plt.show()

def tol_fit():
	data = np.loadtxt("../ignored/misc/ballast_use.txt")
	N = int(data.size/2)
	bal_ = data[:N]
	tol_ = data[N:]/1000
	x = np.linspace(tol_[0],tol_[-1],1000)
	c = np.polyfit(tol_,bal_,7)
	y = np.polyval(c,x)
	print(c)
	np.save("../ignored/misc/ballast_coeffs",c)
	b = 0.03*750/tol_ + 0.04;
	plt.plot(tol_,bal_)
	#plt.plot(x,y)
	plt.plot(tol_,b)
	plt.show()

def test_param(param,val_,N):
	balavg_ = np.zeros(val_.size)
	for j in range(N):
		p = Pool(8)
		bal_ = np.array(p.map(frun(param),val_))
		#plt.plot(tol_,bal_,"g*")
		balavg_ += bal_/N
		print(bal_[0])
	#np.savetxt("../ignored/misc/ballast_use.txt",np.hstack((balavg_,val_)))
	plt.plot(val_,balavg_)
	plt.grid()
	plt.ylabel("ballast use (g/hr)")
	plt.xlabel(param)	
	plt.show()

tol_fit()
