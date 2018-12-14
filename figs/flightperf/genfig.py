import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# This is for creating the a figure for a real flight test of lasanga performance.  

dpath = "../../../flightdata/"
plt.figure(figsize=(15,3))
df = pd.read_hdf(dpath + "smoldf/ssi71_smol.h5")
[print(k) for k in df.keys()]
v = lambda s : df[s][st:-en]
st = int(60*60*8.1)
en = int(60*60*35.1)
t_ = df.index[st:-en]
h = (df["altitude_barometer"]/1000)[st:-en].fillna(method="backfill").values
hg = (df["altitude_gps"]/1000)[st:-en].fillna(method="backfill").values
hf = gaussian_filter1d(h,400)
good = np.where(np.abs(hf-h) < 0.15)[0]
hr = np.interp(np.arange(h.size),good,h[good])
hrf = gaussian_filter1d(hr,60)
 
for i,t in enumerate(df.index.values[st:-en][np.nonzero(np.diff(v('ballast_time_total')))]):
	if i==0: plt.plot(t,0,c='b',alpha=1,label="ballast event", marker='|', linestyle='None', markersize=10, markeredgewidth=1.5) 
	plt.axvline(t,c='b',alpha=0.1) 
for i,t in enumerate(df.index.values[st:-en][np.nonzero(np.diff(v('valve_time_total')))]):
	if i==0: plt.plot(t,0,c='g',alpha=1,label="valve event", marker='|', linestyle='None', markersize=10, markeredgewidth=1.5) 
	plt.axvline(t, c='g',alpha=0.1)
plt.plot(v("las_ss_error_thresh")/1000 + v("las_h_cmd")/1000,c="gray",alpha=0.7, label="contoller bound")
plt.plot(-v("las_ss_error_thresh")/1000 + v("las_h_cmd")/1000,c="gray",alpha=0.7)
plt.plot(t_,h, label="raw measurements", c='C1', alpha=0.6)
plt.plot(t_,hrf, label="balloon altitude",c='C0')
plt.ylabel("altitde (km)")
plt.xlabel("timestamp (UTC 2018)")
plt.legend(loc='uppper right', bbox_to_anchor=(1, 1),ncol=1)
plt.ylim([12,15])
plt.xlim([t_[0],t_[-1]])
plt.tight_layout()
plt.savefig("fig.png")

