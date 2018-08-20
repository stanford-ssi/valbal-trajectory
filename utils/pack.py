import pygrib
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
import scipy.interpolate

from datetime import tzinfo, timedelta, datetime

class UTC(tzinfo):
    def utcoffset(self, dt):
        return timedelta(0)
    def tzname(self, dt):
        return "UTC"
    def dst(self, dt):
        return timedelta(0)
utc = UTC()

tree = None

varlist = ['U component of wind','V component of wind']
vv = ['uvel','vvel']

def read_file(cur):
    F = 'proc/%d.bin' % cur
    if os.path.exists(F):
        print('skipping',F)
        return
    fname = 'data/%d.grb2' % cur
    print('opening',fname)
    f = pygrib.open(fname) # open

    print(list(zip(varlist,vv)))
    assert len(varlist) == len(vv)

    d = f.select()
    #for x in d:
    #    print(x)
    #exit()
    grbs = []
    for x in varlist:
        if type(x) == int:
            grbs.append([d[x-1]])
        else:
            try: grbs.append(f.select(name=x))
            except:
                print('Not including',x)
                grbs.append(None)
            #print(grbs[-1])
    #exit()

    lats, lons = grbs[0][0].latlons()
    print(lats, lons)
    print(lats.shape, lons.shape)
    out = np.zeros((lats.shape[0], lats.shape[1], 8, 2), dtype=np.int16)
    for grb in grbs:
        i = 0
        for row in grb:
            if row['typeOfLevel'] == 'isobaricInhPa' and row.level >= 50 and row.level <= 350:
                j = 0 if row.name.startswith('U') else 1
                if i == 4:
                    print(row.values[42,120])
                    print(lats[42,120], lons[42,120])
                out[:,:,i,j] = row.values*100
                i += 1
    out.flatten().tofile('proc/%d.bin' % cur)
    del f
 

if __name__ == '__main__':
    files = sorted(map(lambda x: int(x.split('.')[0]), os.listdir('data')))
    
    for cur in files: read_file(cur)
