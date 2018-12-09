from datetime import datetime  
from datetime import timedelta
import urllib.request as url
import re
from atmotools import *

index = str(url.urlopen("http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/").read()).split("\\n")
times = []
for row in index:
	timestr = re.findall('(?<=href\=\"gfs\.).*(?=/\")',row)
	if timestr:
		times.append(timestr[0])
list.sort(times,key=int)
for t in times[1:]:
	time = datetime.strptime(t,"%Y%m%d%H")
	procWindData(time, time + timedelta(hours=200),db="gfs_pred_0deg5",dry_run=False)