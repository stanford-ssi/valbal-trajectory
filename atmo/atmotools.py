from datetime import datetime as dt 
from datetime import timedelta
import pygrib as gb
import numpy as np 
import matplotlib.pyplot as plt
import math
from urllib.request import urlretrieve
import os
import pickle

s = "2017-05-21_4"

time = dt.strptime(s,"%Y-%m-%d_%H")
time=time.replace(hour=math.floor(time.hour/3)*3,minute=0,second=0)
print(time+timedelta(1/8))

def fetchGFSAnalysisData(start,end):
	""" Fetch data from gfs database for times inbetween start and end
	"""
	remote="ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/"
	local="../ignored/GFS-anl-0deg5/"
	t = dt.strptime(start,"%Y-%m-%d_%H");
	t = t.replace(hour=math.floor(t.hour/3)*3,minute=0,second=0)
	end_t = dt.strptime(end,"%Y-%m-%d_%H");
	print("Downloading files")
	filelist = []
	times = []
	while t < end_t:
		datestr = "%04d%02d%02d" % (t.year, t.month, t.day) 
		dpath = "%04d%02d/"% (t.year, t.month) + datestr + "/"
		fpath = "gfsanl_3_" + datestr + "_%02d00"%t.hour
		for i in [0,3,6]:
			fulfpath = fpath + "_00%d"%i + ".grb2"
			path = dpath + fulfpath;
			if os.path.exists(local+fulfpath):
				print("Local file "+local+fulfpath+" found, skipping")
			else:	
				print("Fetching from",remote+path+"...",end='',flush=True)
				urlretrieve(remote+path,local+fulfpath)
				print("   done")
			filelist.append(fulfpath)
			times.append(t)
		t +=  timedelta(1/4)
	return filelist,times

def makeWindArray(start,end):
	""" gets data and returns in a mulitdimentional array of format:
		[lon,lat,alt,time,pred,vaules]
	"""
	local="../ignored/GFS-anl-0deg5/"
	filepath = "../ignored/GFS_anl_0deg5_objs/" + start + "_to_" +end + ".pickle"

	if os.path.exists(filepath): print("File " +filepath+ " already exists, skipping") ; return
	windobj = {}
	files,times = fetchGFSAnalysisData(start,end)
	grb = gb.open(local+files[0])
	lat,lon = grb.select(shortName="u",typeOfLevel='isobaricInhPa',level=250)[0].latlons() #arbitrary data, doing this for latlons
	windobj["lats"] = lat[:,0] 
	nlats = lat[:,0].size
	windobj["lons"] = lon[0,:]
	nlons = lon[0,:].size
	windobj["values"] = ["u","v"]
	nvalues = 2;
	levels = getGRIBlevels(grb,windobj["values"][0])
	windobj["levels"] = levels
	alts = p2a(levels)
	windobj["alts"] = alts
	nalts = alts.size
	timeshr = np.unique([(t - times[0]).seconds/3600 for t in times])
	ntimes = timeshr.size
	windobj["times"] = timeshr
	windobj["start_time"] = times[0]
	tdelatshr = np.array([0.0,3.0,6.0]);
	windobj["tdeltas"] = tdelatshr;
	ntdeltas = 3;
	data = np.zeros((nlons,nlats,nalts,ntimes,ntdeltas,nvalues))
	dtind = 0
	tind = 0
	print("Starting to build array. This may take a while :(")
	total = nvalues*ntdeltas*ntimes
	ctr = 0
	for i in range(len(files)):
		path = local+files[i]
		grb = gb.open(path)
		for valind,val in enumerate(windobj["values"]):
			dat = grb.select(shortName=val,typeOfLevel='isobaricInhPa',level=levels)
			data[:,:,:,tind,dtind,valind] = np.array(list(map(lambda x : x.values,dat))).T
			ctr+=1
			print("%d/%d"%(ctr,total))
		dtind += 1
		if dtind == 3 : tind += 1; dtind=0
	windobj["data"] = data
	print("done")
	print("Saving to ",filepath)
	pickle.dump(windobj,open(filepath,'wb'))

def getGRIBlevels(grib,shortname='v'):
	""" retruns array of all levels in a GRIB file
	"""
	levels = []
	for message in grib:
		if message.shortName == shortname:
			if message.typeOfLevel == "isobaricInhPa":
				levels.append(message.level)
	levels = np.unique(levels)
	return levels

def p2a(x):
	""" Presure to altitude hectopascals to meters 
	"""
	return (1-(x/1013.25)**0.190284)*145366.45


wind =  pickle.load(open("../ignored/GFS_anl_0deg5_objs/2018-05-1_0_to_2018-05-2_0.pickle","rb"))
print(wind)

exit()
makeWindArray("2018-05-1_0","2018-05-2_0")

grb = gb.open("../ignored/gfsanl_3_20060407_0000_000.grb")


Udat = grb.select(shortName="u",typeOfLevel='isobaricInhPa',level=250)[0]
Vdat = grb.select(shortName="v",typeOfLevel='isobaricInhPa',level=250)[0]
U = Udat.values
V = Vdat.values
print(V)
lat, lon = Udat.latlons()
lat = lat[90-50:90-25,:45] 
lon = lon[:25,360-115:360-70]
plt.quiver(lon,lat,U[90-50:90-25,360-115:360-70],V[90-50:90-25,360-115:360-70])
plt.show()