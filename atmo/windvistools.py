import numpy as np 
import atmotools as at 
import pandas as pd 
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 

m = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=50,\
            llcrnrlon=-130,urcrnrlon=-65,resolution='l')
m.drawcoastlines()
m.drawcountries()
plt.show()

def plotWindOverMap(region,windobj):
	""" Region is a list of format [lon_center, lat_center, lon_width, lat_width]
		Windobj is an dict of wind data of the format output by makeWindArray in atmotools
	"""
	pass
	m = Basemap(llcrnrlon=-120,llcrnrlat=22,urcrnrlon=-60,urcrnrlat=55,
        projection='lcc',lat_1=33,lat_2=45,lon_0=-95,resolution='l')





