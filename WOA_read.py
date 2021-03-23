#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++11/02/2021++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
#imports
from time import time
from copy import deepcopy
import os
import numpy as np
from sklearn.utils import Bunch
from netCDF4 import Dataset
from trapezoidal import regression_points_calc
import pickle
import sqlite3
from random_sampler import *

def read_WOA(woa_file,decimal_Latitude_range="14,40", decimal_Longitude_range="-35,-9"):#3D->3D
    """
requires:
    woa_file    file adress of file with World Ocean Atlas data.
                .csv with 2 first columns as latitude,longitude,
                rest in the row - values in deeper depths at
                coordinate.
                Depths:
                    0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,
                    75,80,85,90,95,100,125,150,175,200,225,250,
                    275,300,325,350,375,400,425,450,475,500,550,
                    600,650,700,750,800,850,900,950,1000,1050,
                    1100,1150,1200,1250,1300,1350,1400,1450,1500,
                    1550,1600,1650,1700,1750,1800,1850,1900,1950,
                    2000,2100,2200,2300,2400,2500,2600,2700,2800,
                    2900,3000,3100,3200,3300,3400,3500,3600,3700,
                    3800,3900,4000,4100,4200,4300,4400,4500,4600,
                    4700,4800,4900,5000,5100,5200,5300,5400,5500
ensures:3D numpy matrix(array), array with latitudes, array with
    longitudes, array with depths and size of grid
"""
    import gzip
    nulls = -9999
    bathimetry_spacing = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500]

    with gzip.open(woa_file,'rb') as handle:
        data = [[nulls if j=='' else float(j) for j in i[:-1].decode('utf-8').split(',')] for i in handle.readlines()[2:]]
    bunch = Bunch()
    bunch['ygrid'] = np.array(sorted(list(set([i[0] for i in data]))))
    bunch['xgrid'] = np.array(sorted(list(set([i[1] for i in data]))))
    bunch['depth'] = np.array(bathimetry_spacing, dtype='float32')


    longitude_min, longitude_max = [float(i) for i in decimal_Longitude_range.split(',')]
    latitude_min, latitude_max = [float(i) for i in decimal_Latitude_range.split(',')]

    lon_left = bunch.xgrid.searchsorted(longitude_min)
    lon_right = bunch.xgrid.searchsorted(longitude_max)
    lat_bottom = bunch.ygrid.searchsorted(latitude_min)
    lat_top = bunch.ygrid.searchsorted(latitude_max)

    bunch.xgrid = bunch.xgrid[lon_left:lon_right]
    bunch.ygrid = bunch.ygrid[lat_bottom:lat_top]


    out = np.array([[[[nulls]*len(bunch.depth) for x in bunch.xgrid] for y in bunch.ygrid]])
    for i in data:
        if not latitude_min<=i[0]<=latitude_max:
            continue
        if not longitude_min<=i[1]<=longitude_max:
            continue
        y = bunch.ygrid.searchsorted(i[0])
        x = bunch.xgrid.searchsorted(i[1])
        for z,e in enumerate(i[2:]):
            out[0,y,x,z] = e

    bunch['grid'] = abs(bunch.ygrid[0]-bunch.ygrid[1])#double check
    bunch['data'] = out
    return bunch


