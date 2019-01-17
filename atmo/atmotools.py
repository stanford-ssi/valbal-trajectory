
from datetime import datetime  
from datetime import timedelta
import pandas as pd
import pygrib as gb
import numpy as np 
import matplotlib.pyplot as plt
import math
from urllib.request import urlretrieve
import urllib.request as url
import os
import pickle
import sys
import argparse
import re
import json
import hashlib

def getRecentGFS():
	index = str(url.urlopen("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/").read()).split("\\n")
	times = []
	for row in index:
		timestr = re.findall('(?<=href\=\"gfs\.).*(?=/\")',row)
		if timestr:
			times.append(timestr[0])
	list.sort(times,key=int)
	recent=times[-1]
	pred_time = datetime.strptime(recent,"%Y%m%d%H")
	return pred_time

def getFile(src,dst,dry_run=False):
	if os.path.exists(dst):
		size = os.path.getsize(dst)
		print("Local file "+dst+" found and is %.2fMB, skipping"%(size/1000000))
		if size < 50000000:
			raise Warning("file size smaller than expected: %.2fMB "%(size/1000000))
	else:	
		print("Fetching from",src+"...",end='',flush=True)
		if not dry_run:
			urlretrieve(src,dst)
			size = os.path.getsize(dst)
			print("   done, %.2fMB"%(size/1000000))
		else:
			print("   jk dry run.")

def gfsFileToTime(f):
	if "pred" in f:
		base = f.split("/")[-2]
		return datetime.strptime(base,"%Y%m%d_%H") + timedelta(hours=int(f[-3:]))
	else:
		f = f.split("/")[-1]
		date = f.split("_")[2]
		hr = f.split("_")[3][:2]
		plus = f.split("_")[4][:3]
		return datetime.strptime(date,"%Y%m%d") + timedelta(hours=int(hr)) + timedelta(hours=int(plus))


def setupAtmoData(start,end,db='gfs_anl_0deg5',pred_time=None):
	if not type(start) == type("boop"):
		if start == None:
			start = datetime.utcnow();
		if end == None:
			end = datetime.utcnow() + timedelta(hours=100);
		start -= timedelta(hours=1)
		start = "%04d-%02d-%02d_%02d"%(start.year,start.month,start.day,start.hour)
		end = "%04d-%02d-%02d_%02d"%(end.year,end.month,end.day,end.hour)
	if "hr" in start:
		if end == None:
			end = datetime.utcnow() + timedelta(hours=100);
		start = datetime.utcnow() + timedelta(hours=int(start.split("hr")[0]))
		start = "%04d-%02d-%02d_%02d"%(start.year,start.month,start.day,start.hour)
		end = "%04d-%02d-%02d_%02d"%(end.year,end.month,end.day,end.hour)
		print(start)
	t = datetime.strptime(start,"%Y-%m-%d_%H");
	if int(np.round(t.hour/6)*6) > 23:
		t += timedelta(hours=24)
	t = t.replace(hour=int(np.round(t.hour/6)*6)%24,minute=0,second=0)
	end_t = datetime.strptime(end,"%Y-%m-%d_%H") + timedelta(1/4);
	start = "%04d-%02d-%02d_%02d"%(t.year,t.month,t.day,t.hour)
	end = "%04d-%02d-%02d_%02d"%(end_t.year,end_t.month,end_t.day,end_t.hour)
	if db == 'gfs_anl_1deg':
		remote="https://nomads.ncdc.noaa.gov/data/gfsanl/"
		local="../ignored/raw/gfs_anl_1deg/"
		fstartname = "gfsanl_3_"
	if db == 'gfs_anl_0deg5':
		remote="https://nomads.ncdc.noaa.gov/data/gfsanl/"
		local="../ignored/raw/gfs_anl_0deg5/"
		fstartname = "gfsanl_4_"
	if db == "euro_fc":
		raise Warning("lol good luck son, this data is not easy to get")
		return 
	if db == "gfs_pred_0deg5":
		# If no pred time is passed, we use the most recent data. If the start time, t, is less than the 
		# the predition time, we use that instead
		if pred_time == None:
			pred_time = getRecentGFS()
		else:
			pred_time = datetime.strptime(pred_time,"%Y-%m-%d_%H")
		pred_time = pred_time.replace(microsecond=0,second=0,minute=0,hour=pred_time.hour-(pred_time.hour % 6))
		if t<pred_time:
			pred_time = t
		remote = "http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.%04d%02d%02d%02d/"%(pred_time.year,pred_time.month,pred_time.day,pred_time.hour)
		db = "gfs_pred_0deg5/%04d%02d%02d_%02d"%(pred_time.year,pred_time.month,pred_time.day,pred_time.hour)
		local = "../ignored/raw/"+db+"/"
		fstartname = "gfs.t%02dz.pgrb2.0p50.f"%pred_time.hour
	return 	start,end,t,end_t,db,local,remote,fstartname,pred_time



def fetchWindData(start,end,db='gfs_anl_0deg5',pred_time=None,dry_run=False):
	""" Fetch data from gfs database for times inbetween start and end
	"""
	start,end,t,end_t,db,local,remote,fstartname,pred_time = setupAtmoData(start,end,db=db,pred_time=pred_time)
	if not os.path.exists(local):
		os.makedirs(local)
	print("Downloading files from " + remote)
	print("Saving to "+local)
	filelist = []
	times = []
	if 'gfs_anl_' in db:
		while t < end_t:
			datestr = "%04d%02d%02d" % (t.year, t.month, t.day) 
			dpath = "%04d%02d/"% (t.year, t.month) + datestr + "/"
			fpath = fstartname + datestr + "_%02d00"%t.hour
			for i in [0,3]:
				fulfpath = fpath + "_00%d"%i + ".grb2"
				path = dpath + fulfpath;
				getFile(remote+path,local+fulfpath,dry_run=dry_run)
				filelist.append(fulfpath)
				times.append(t + i*timedelta(1/24))
			t +=  timedelta(1/4)

	if 'gfs_pred_' in db:
		while t < end_t:
			dhours = int(((t - pred_time).total_seconds()/60/60))
			if dhours > 372:
				break
			fpath = fstartname + "%03d"%dhours
			getFile(remote+fpath,local+fpath,dry_run=dry_run)
			times.append(t)
			filelist.append(fpath)
			if dhours < 240:
				t += timedelta(1/8)
			else:
				t += timedelta(1/2)


	return filelist,times,start,end,db


def procWindData(start,end,db='gfs_anl_0deg5',overwrite=False,pred_time=None,altitude_range = [0,30000],aux_data=False,dry_run=False,files=[]):
	if files:
		db = "/".join(files[0].split("/")[:-1]).split("ignored/raw/")[-1]
		times = list(map(gfsFileToTime,files))
		for i in range(len(files)):
			files[i] = files[i].split("/")[-1]
		#print(files[0],times[0])
		#print(db)
		ret = [files,times,0,0,db]	
		if not db:
			raise Warning("bad file paths")
	else:
		ret = fetchWindData(start,end,db,pred_time=pred_time,dry_run=dry_run)
	if dry_run:
		print("quitting, dry run")
		return
	db = ret[4]
	dstpath = "../ignored/proc/" + db + "/"
	auxpath = "../ignored/proc/" + db + "_aux/"
	srcpath = "../ignored/raw/" + db + "/"
	times = ret[1]
	files = ret[0]
	if not os.path.exists(dstpath):
		os.makedirs(dstpath)
	if aux_data and not os.path.exists(auxpath):
		os.makedirs(auxpath)
	grb = gb.open(srcpath+files[0])
	lats,lons = grb.select(shortName="u",typeOfLevel='isobaricInhPa',level=250)[0].latlons() #arbitrary data, doing this for latlons
	lats = lats[:,0]
	lons = lons[0,:]
	levels = getGRIBlevels(grb,altitude_range=altitude_range)
	headertext,jsontext,check = genWindHeader(db.replace("/","_"),lons,lats,levels)
	headerfile = dstpath + db.replace("/","_") + ".h"
	jsonfile = dstpath + "config.json"
	procfiles=[]
	with open(headerfile,"w") as f:
		f.write(headertext)
	with open(jsonfile,"w") as f:
		f.write(jsontext)
	keys = {"lons":lons, "lats":lats, "levels": levels, "alts": p2a(levels)}
	pickle.dump(keys,open(dstpath + "keys.pickle",'wb'))
	print("Altitudes: ",keys["alts"])
	for k,file in enumerate(files):
		if db.split("/")[0]=="gfs_pred_0deg5":
			# check file timestamps
			timestr = db.split("/")[1]
			ftime = datetime.strptime(timestr,"%Y%m%d_%H")
			assert(ftime+timedelta(hours=int(file[-3:]))==times[k])
		outpath = dstpath + '%d.bin' % (times[k]-datetime(1970,1,1)).total_seconds()
		auxoutpath = auxpath + '%d.bin' % (times[k]-datetime(1970,1,1)).total_seconds()
		#print(times[k]) 
		#exit()
		procfiles.append(outpath)
		if os.path.exists(outpath) and not overwrite:
			size = os.path.getsize(outpath)
			if np.zeros((lats.size,lons.size,levels.size,2),dtype=np.int16).size*2 + 4==size: #+4 is for the size of the hash
				print("Local file "+outpath+" found and is %.2fMB"%(size/10**6)+", skipping (%d / %d)" % (k+1, len(files)))
				continue
		print("Saving to",outpath+"...",end='',flush=True)
		data = np.zeros((lats.size,lons.size,levels.size,2),dtype=np.int16)
		if aux_data:
			aux = np.zeros((lats.size,lons.size,levels.size,1),dtype=np.int16)
		path = srcpath+file
		grb = gb.open(path)
		i = 0.
		for row in grb:
			last_level = 0
			if row.shortName == 'u' or row.shortName == 'v':
				if row['typeOfLevel'] == 'isobaricInhPa' and row.level >= levels[0] and row.level <= levels[-1]:
					j = 0 if row.name.startswith('U') else 1
					data[:,:,int(i),j] = row.values*100
					i += 0.5
					assert row.level > last_level , "Levels in grib not increasing"
					assert (j) or (i % 1), "U, V order mismatch in grib"
					last_level = row.level
			if aux_data:
				if row.shortName == "t" and row.level >= levels[0] and row.level <= levels[-1]:
					data[:,:,levels.searchsorted(row.level),0] = row.values
		f = open(outpath,"wb")
		f.write(check)
		f.write(data.flatten().tobytes())
		if aux_data:
			aux.flatten().tofile(auxoutpath)
		print("   done (%d / %d)" % (k+1, len(files)))
	return procfiles,times

def getDataValue(data,time,db):
	srcpath = "../ignored/proc/" + db + "/"	
	keys = pickle.load(open(srcpath+'keys.pickle','rb'))
	ftimes = np.array(list(map(lambda x : int(x.split('.bin')[0]) if len(x.split('.bin'))==2 else 0 ,os.listdir(srcpath))))
	nearest = ftimes[np.argmin(np.abs(ftimes-time))]
	print(time,nearest)



def genWindHeader(dataset,lons,lats,levels):

	filetext = "#ifndef __%s__ \n\r" % dataset
	filetext += "#define __%s__ \n\r \n\r" % dataset
	filetext += "/* \n\r"
	filetext += " * This header is auto-generated by atmotools.py. \n\r"
	filetext += " * For the dataset '%s'. \n\r" % dataset
	filetext += " */ \n\r \n\r"
	filetext += "typedef short wind_t; // Type used to store wind data on disk, in cm/s. \n\r \n\r"
	filetext += "/* wind grid paramters */ \n\r"
	filetext += "const float LON_MIN = %f;\n\r" % lons[0]
	filetext += "const float LON_MAX = %f;\n\r" % lons[-1]
	filetext += "const float LON_D = %f;\n\r" % (lons[1] - lons[0])
	filetext += "const int NUM_LONS = %d;\n\r" % lons.size
	filetext += "const float LAT_MIN = %f;\n\r" % lats[0]
	filetext += "const float LAT_MAX = %f;\n\r" % lats[-1]
	filetext += "const float LAT_D = %f;\n\r" % (lats[1] - lats[0])
	filetext += "const int NUM_LATS = %d;\n\r" % lats.size
	filetext += "const float LEVELS[] = {" + ",".join(map(str,levels*100)) + "}; \n\r"
	filetext += "const int NUM_LEVELS = sizeof(LEVELS)/sizeof(LEVELS[0]); \n\r"
	filetext += "const int NUM_VARIABLES = 2; \n\r \n\r"
	filetext += "#endif \n\r"

	conf = {"LON_MIN":lons[0],"LON_MAX":lons[-1],
	"LON_D":(lons[1] - lons[0]),
	"NUM_LONS":lons.size,
	"LAT_MIN":lats[0],
	"LAT_MAX":lats[-1],
	"LAT_D":(lats[1] - lats[0]),
	"NUM_LATS":lats.size,
	"LEVELS":"[" + ",".join(map(str,levels*100)) + "]",
	"NUM_LEVELS":levels.size,
	"NUM_VARIABLES":2,
	}
	jtxt = str(json.dumps(conf)).encode('utf-8')
	check = hashlib.sha1(jtxt)
	conf["hash"] = check.hexdigest()[:8]
	jtxt = json.dumps(conf)
	check = check.digest()[:4]
	return filetext,jtxt,check;



def makeWindArray(start,end,overwrite=False,altitude_range=[10000,20000],name_modifier=''):
	""" gets data and returns in a mulitdimentional array of format:
		[lon,lat,alt,time,pred,vaules]
	"""
	files,times,start,end = fetchGFSAnalysisData(start,end)
	local="../ignored/GFS-anl-0deg5/"
	filepath = "../ignored/GFS_anl_0deg5_objs/" + start + "_to_" +end + name_modifier+ ".pickle"
	filepath2 = "../ignored/GFS_anl_0deg5_objs/" + start + "_to_" +end + name_modifier+".npy"
	if os.path.exists(filepath) and not overwrite: print("File " +filepath+ " already exists, skipping") ; return filepath,filepath2
	windobj = {}
	grb = gb.open(local+files[0])
	lat,lon = grb.select(shortName="u",typeOfLevel='isobaricInhPa',level=250)[0].latlons() #arbitrary data, doing this for latlons
	windobj["lats"] = lat[:,0] 
	nlats = lat[:,0].size
	windobj["lons"] = lon[0,:]
	nlons = lon[0,:].size
	windobj["values"] = ["u","v"]
	nvalues = 2;
	levels = getGRIBlevels(grb,windobj["values"][0])
	alts = p2a(levels)
	idxs = np.where(np.logical_and(altitude_range[0]<=alts,alts<=altitude_range[1]))
	levels = levels[idxs]
	alts = alts[idxs]
	windobj["levels"] = levels
	windobj["alts"] = alts
	nalts = alts.size
	timeshr = np.unique([(t - times[0]).seconds/3600 + (t - times[0]).days*24 for t in times])
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
		for valind,val in enumerate(windobj["values"]): #this could be sped up by a factor of 2 if u selected both vals are once. But then have to be sure they are returned in the correct order
			dat = grb.select(shortName=val,typeOfLevel='isobaricInhPa',level=levels)
			data[:,:,:,tind,dtind,valind] = np.array(list(map(lambda x : x.values,dat))).T
			ctr+=1
			print("%d/%d"%(ctr,total))
		dtind += 1
		if dtind == 3 : tind += 1; dtind=0
	print("done")
	print("Saving to ",filepath)
	pickle.dump(windobj,open(filepath,'wb'))
	print("Saving to ",filepath2)
	np.save(filepath2,data)
	return filepath,filepath2

def getGRIBlevels(grib,shortname='v',altitude_range = [0,25000]):
	""" retruns array of all levels in a GRIB file
	"""
	levels = []
	for message in grib:
		if message.shortName == shortname:
			if message.typeOfLevel == "isobaricInhPa":
				levels.append(message.level)
	levels = np.unique(levels)
	alts = p2a(levels)
	idxs = np.where(np.logical_and(altitude_range[0]<=alts,alts<=altitude_range[1]))
	levels = levels[idxs]
	return levels

def p2a(x):
	""" Presure to altitude hectopascals to meters 
	"""
	return (1-(x/1013.25)**0.190284)*145366.45*0.3048

def getKeysForBin(filepath):
	keypath = "/".join(filepath.split("/")[:-1])+"/keys.pickle"
	return pickle.load(open(keypath,'rb'))

def getArrayFromBin(filepath,keys):
	return np.fromfile(filepath,dtype=np.int16).reshape(keys["lats"].size,keys["lons"].size,keys["levels"].size,2)

'''
df = pd.read_hdf('../../valbal-controller/ssi67-analysis/ssi67.h5')
print(df.long_gps.values[1000])
print(df.lat_gps.values[1000])
print(df.index[0])
procWindData(df.index[0],df.index[-1] + timedelta(2),db="gfs_anl_0deg5",overwrite=False)
'''

'''
#reprocess old files
for di in os.listdir("../ignored/raw/gfs_pred_0deg5/"):
	#di = "../ignored/raw/gfs_anl_0deg5/"
	if datetime.strptime(di,"%Y%m%d_%H") < datetime.strptime("20181120_00","%Y%m%d_%H"):
		continue
	print(di)
	files = os.listdir("../ignored/raw/gfs_pred_0deg5/"+di+"/")
	for i in range(len(files)):
		files[i] = "../ignored/raw/gfs_pred_0deg5/"+ di+"/"+files[i] 
	procWindData("0","0",files=files,overwrite=False)
'''

#fetchWindData("2018-10-20_00","2018-10-22_00",db="gfs_pred_0deg5")
#procWindData("2018-10-28_12","2018-11-02_00",db="gfs_pred_0deg5",overwrite=True)
#procWindData("2018-10-20_18","2018-10-23_18",db="gfs_pred_0deg5",overwrite=True,altitude_range = [0,30000])
#procWindData("2018-10-28_12","2020-11-01_00",db="gfs_pred_0deg5",overwrite=False)


if __name__== "__main__":
	parser = argparse.ArgumentParser(description="Tools for handling atmosphere data")
	parser.add_argument("activity",type=str,help="activity to perform")
	parser.add_argument("--start",type=str,help="start time for wind files")
	parser.add_argument("--end",type=str,help="start time for wind files")
	parser.add_argument("--db",type=str,help="wind database name",default="gfs_anl_0deg5")
	parser.add_argument("--pred_time",type=str,help="time predictions were made on",default=None)
	parser.add_argument("-o","--overwrite",help="if processed files should be over written",action='store_true')
	args = parser.parse_args()
	if args.activity=="proc":
		procWindData(args.start,args.end,db=args.db,pred_time=args.pred_time,overwrite=args.overwrite)
