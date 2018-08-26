import urllib.request
import os
from datetime import tzinfo, timedelta, datetime

class UTC(tzinfo):
    def utcoffset(self, dt):
        return timedelta(0)
    def tzname(self, dt):
        return "UTC"
    def dst(self, dt):
        return timedelta(0)
utc = UTC()

launch_time = datetime(2017, 12, 9, 21, 39, tzinfo=utc)
#launch_time = datetime(2018, 5, 12, 20, 0, tzinfo=utc)
#launch_time = datetime(2018, 6, 14, 17, 0, tzinfo=utc)
print((launch_time-datetime(1970,1,1,tzinfo=utc)).total_seconds())
landing_time = datetime(2017, 12, 14, 23, 13, tzinfo=utc)
print((landing_time-launch_time).total_seconds())
#landing_time = datetime(2018, 5, 16, 21, 0, tzinfo=utc)
#landing_time = datetime(2018, 6, 18, 17, 0, tzinfo=utc)

times = ['0000_000',
         '0000_001',
         '0000_002',
         '0000_003',
         '0000_003',
         '0000_003',
         '0600_000',
         '0600_001',
         '0600_002',
         '0600_003',
         '0600_003',
         '0600_003',
         '1200_000',
         '1200_001',
         '1200_002',
         '1200_003',
         '1200_003',
         '1200_003',
         '1800_000',
         '1800_001',
         '1800_002',
         '1800_003',
         '1800_003',
         '1800_003']

times = ['0000_000',
         '0000_000',
         '0000_000',
         '0000_003',
         '0000_003',
         '0000_003',
         '0600_000',
         '0600_000',
         '0600_000',
         '0600_003',
         '0600_003',
         '0600_003',
         '1200_000',
         '1200_000',
         '1200_000',
         '1200_003',
         '1200_003',
         '1200_003',
         '1800_000',
         '1800_000',
         '1800_000',
         '1800_003',
         '1800_003',
         '1800_003']

assert len(times) == 24
it = launch_time
done = set()

while it < (landing_time + timedelta(hours=3)):
    url = "https://nomads.ncdc.noaa.gov/data/namanl/%04d%02d/%04d%02d%02d/namanl_218_%04d%02d%02d_%s.grb2" % (it.year, it.month, it.year, it.month, it.day, it.year, it.month, it.day, times[it.hour])
    url = "https://nomads.ncdc.noaa.gov/data/gfsanl/%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%s.grb2" % (it.year, it.month, it.year, it.month, it.day, it.year, it.month, it.day, times[it.hour])
    # 
    sp = times[it.hour].split('_')
    t = datetime(it.year, it.month, it.day, int(int(sp[0])/100+int(sp[1])))
    if url not in done:
        dt = int((t-datetime(1970,1,1)).total_seconds())
        fname = "../ignored/raw/GFS_anl_0deg5"+str(dt)+".grb2"
        if not os.path.exists(fname):
            print("Downloading",t)
            done.add(url)
            with urllib.request.urlopen(url) as u:
                with open(fname,"wb") as f:
                    f.write(u.read())
            
    it += timedelta(hours=1)
