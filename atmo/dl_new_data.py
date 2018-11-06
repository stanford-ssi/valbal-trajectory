import time
from atmotools import *

previous = getRecentGFS()
path = setupAtmoData(previous, previous+ timedelta(hours=500),db="gfs_pred_0deg5")[5]
if not os.path.exists(path):
	print(path)
	procWindData(previous, previous+timedelta(hours=500),db="gfs_pred_0deg5")
while(1):
	print("[%s] Checking for new data... " % datetime.now(),end='',flush=True)
	newest = getRecentGFS()
	if newest > previous:
		print("New data folder found")
		try:
			procWindData(newest, newest + timedelta(hours=500),db="gfs_pred_0deg5")
		except:
			print("Files not ready :(")
	else:
		print("No new data :(")
	time.sleep(60*1)