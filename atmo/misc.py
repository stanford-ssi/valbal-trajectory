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

path = "../ignored/raw/gfs_anl_0deg5/gfsanl_4_20171212_0600_000.grb2"
grib = gb.open(path)
for g in grib:
	print(g)
	if g.shortName == "t":
		print(g.level,"t")

print()