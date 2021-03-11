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

#-----------------------------------------------------------------------

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


#---------------------------------
def read_copernicus_phys(nc_file_phys):
    """
requires:
    file name of netCDF4 with physical variables from copernicus
ensures:
    Bunch object
"""
    print("Physics, File: {}".format(nc_file_phys))
    #load data
    fh = Dataset(nc_file_phys)
    #extraction of data
    b1 = Bunch()
    b1['time'] = fh.variables['time'][:].data
    b1['depth'] = fh.variables['depth'][:].data
    b1['ygrid'] = fh.variables['latitude'][:].data
    b1['xgrid'] = fh.variables['longitude'][:].data
    b1['grid'] = abs(b1['ygrid'][0]-b1['ygrid'][1])
    
    #b1['northV'] = fh.variables['vo'][:].data
    b1['temperature'] = fh.variables['thetao'][:].data
    #b1['eastV'] = fh.variables['uo'][:].data
    b1['salinity'] = fh.variables['so'][:].data
    fh.close()
    return b1


def read_copernicus_PP(nc_file_PP):
    """
requires:
    file name of netCDF4 with primary production from copernicus
ensures:
    Bunch object
"""
    print('Primary Production, File: {}'.format(nc_file_PP))
    #load data
    fh = Dataset(nc_file_PP)
    #extraction of data
    b2 = Bunch()
    b2['time'] = fh.variables['time'][:].data
    b2['depth'] = fh.variables['depth'][:].data
    b2['ygrid'] = fh.variables['latitude'][:].data
    b2['xgrid'] = fh.variables['longitude'][:].data
    b2['grid'] = abs(b2['ygrid'][0]-b2['ygrid'][1])
    
    b2['O2'] = fh.variables['o2'][:].data
    b2['pH'] = fh.variables['ph'][:].data
    b2['chla'] = fh.variables['chl'][:].data
    b2['netPP'] = fh.variables['nppv'][:].data
    b2['phytoplankyonC'] = fh.variables['phyc'][:].data
    fh.close()
    return b2


def read_copernicus_Zoo(nc_file_zoo):
    """
requires:
    file name of netCDF4 with zooplankton from copernicus
ensures:
    Bunch object
"""
    print('Zooplankton, File: {}'.format(nc_file_zoo))
    #load data
    fh = Dataset(nc_file_zoo)
    #extraction of data
    b3 = Bunch()
    b3['time'] = fh.variables['time'][:].data
    b3['ygrid'] = fh.variables['latitude'][:].data
    b3['xgrid'] = fh.variables['longitude'][:].data
    b3['grid'] = abs(b3['ygrid'][0]-b3['ygrid'][1])
    #epipelagic
    b3['depth_epi'] = fh.variables['depth_epi'][:].data
    b3['mass_conc_epi'] = fh.variables['mnkc_epi'][:].data
    #upper mesopelagic
    b3['depth_umeso'] = fh.variables['depth_umeso'][:].data
    b3['mass_conc_umeso'] = fh.variables['mnkc_umeso'][:].data
    b3['mass_conc_umeso_mig'] = fh.variables['mnkc_ummeso'][:].data
    #lower mesopelagic
    b3['depth_lmeso'] = fh.variables['depth_lmeso'][:].data
    b3['mass_conc_lmeso'] = fh.variables['mnkc_lmeso'][:].data
    b3['mass_conc_lmeso_mig'] = fh.variables['mnkc_lmmeso'][:].data
    b3['mass_conc_lmeso_Hmig'] = fh.variables['mnkc_lhmmeso'][:].data
    fh.close()
    return b3


def read_copernicus(nc_file_phys, nc_file_PP, nc_file_zoo):
    """
requires:
    nc_file_phys   common file name of netCDF4 with physical variables from copernicus
    nc_file_PP     common file name of netCDF4 with primary production from copernicus
    nc_file_zoo    common file name of netCDF4 with zooplankton from copernicus
ensures:
    list with Bunch objects with all data compiled
"""
    ncfiles = (nc_file_phys, nc_file_PP, nc_file_zoo)
    data = [0,0,0]
    out = [0,0,0]
    for i,ncfile in enumerate(ncfiles):
        if i==0:
            wfunc = read_copernicus_phys
        elif i==1:
            wfunc = read_copernicus_PP
        elif i==2:
            wfunc = read_copernicus_Zoo
        path = '\\'.join(ncfile.split('\\')[:-1])
        if path=='':
            list_files = os.listdir()
        else:
            list_files = os.listdir(path)
        ncfileregex = [j for j in list_files if j.startswith(ncfile.split('\\')[-1])]
        data[i] = [0]*len(ncfileregex)
        for ii,file in enumerate(ncfileregex):
            data[i][ii] = wfunc(os.path.join(path,file))

        print("pack it home")
        out[i] = Bunch()
        out[i]['grid'] = data[i][0].grid
        out[i]['xgrid'] = data[i][0].xgrid
        out[i]['ygrid'] = data[i][0].ygrid
        if i<2:
            out[i]['depth'] = data[i][0].depth
        
        for key in data[i][0]:
            if key not in ['grid','depth','xgrid','ygrid']:
                out[i][key] = np.concatenate([data[i][jj][key] for jj in range(len(data[i]))])
        data[i] = 0
    return out


def averaged_copernicus(data):
    """
requires:
    data    list with 3 bunches of copernicus data
ensures:
    all datasets averaged on time dimension
"""
    for i,bunch in enumerate(data):
        for key in bunch:
            if key not in ['grid','depth','xgrid','ygrid','time']:
                data[i][key] = np.mean(bunch[key], axis=0)


def standardize_grid(data):
    """
requires:
    data    list with 3 bunches of copernicus data
ensures:
    all bunches in data have the same grid size filled with values
"""
    b2 = convert2higherResolutionGrid(data[1], data[0], d_tag=['O2',
                                                               'chla',
                                                               'netPP',
                                                               'pH',
                                                               'phytoplankyonC'])
    b3 = convert2higherResolutionGrid(data[2], data[0], d_tag=['depth_epi',
                                                               'mass_conc_epi',
                                                               'depth_umeso',
                                                               'mass_conc_umeso',
                                                               'mass_conc_umeso_mig',
                                                               'depth_lmeso',
                                                               'mass_conc_lmeso',
                                                               'mass_conc_lmeso_mig',
                                                               'mass_conc_lmeso_Hmig'])
    data[1] = b2
    data[2] = b3


#---------------------------------
def read_georeferenced_observations(obs):
    """species data observations
requires:
    obs    list(list()) each observations entry as:
               ["genus", "species", longitude, latitude, depth]
ensures:
    np.array ready for Bunch
"""
    dt = np.dtype([('species', np.unicode_,36),
                   ('dd long', np.float64),
                   ('dd lat', np.float64),
                   ('m depth', np.int16)])
    out = np.array([("_".join(x[0:2]), *x[2:]) for x in obs], dtype=dt)
    return out


def read_sql_georeferenced_observations(db_path = 'D:\\PhD\\articles\\article DB',db = 'samples.db'):
    """from SQL database file, collect georeferenced entries
requires:
    db_path    path of folder where db file is
               default: 'D:\\PhD\\articles\\article DB'
    db         db file name
               default: 'samples.db'
ensures: Numpy array with 'species' as 'genus_species',
    'dd long' float(),
    'dd lat' float() and 'depth' as int()
"""
    conn = sqlite3.connect(os.path.join(db_path,db))
    cursor = conn.cursor()
    cursor.execute('''
SELECT DISTINCT genus,sp,Longitude,Latitude,depth_min,depth_Max
FROM sample inner join species using (id_s) inner join net using (id_n,id_st,id_c) inner join station using (id_st,id_c) inner join cruise using (id_c)
WHERE id_c=1 AND Longitude<=-9 AND Latitude<=40 AND Longitude>=-35 AND Latitude>=14
ORDER BY genus,sp''')
    Lobs = cursor.fetchall()
    cursor.close()
    conn.close()
    Lobs = [[*i[:-2],round((i[-2]+i[-1])/2)] for i in [i for i in Lobs if type(i[-2])==int and type(i[-1])==int]]
    print('Filtering: Longitude<=-9 AND Latitude<=40 AND Longitude>=-35 AND Latitude>=14')
    return read_georeferenced_observations(Lobs)


#---------------------------------
def read_gebco(filename,decimal_Latitude_range="14,40", decimal_Longitude_range="-35,-9"):
    """GEBCO ('GEBCO_2014_1D.nc')
into matrix of depths
requires: filename[str()]
ensures:np.matrix(2D array), array with latitudes, array with
    longitudes and size of grid
    """
    #load data
    fh = Dataset(filename)
    #extraction of data
    out = Bunch()
    lon_left = fh.variables['lon'][:].searchsorted(float(decimal_Longitude_range.split(',')[0]))
    lon_right = fh.variables['lon'][:].searchsorted(float(decimal_Longitude_range.split(',')[1]))
    lat_bottom = fh.variables['lat'][:].searchsorted(float(decimal_Latitude_range.split(',')[0]))
    lat_top = fh.variables['lat'][:].searchsorted(float(decimal_Latitude_range.split(',')[1]))
    out['xgrid'] = fh.variables['lon'][lon_left:lon_right].data
    out['ygrid'] = fh.variables['lat'][lat_bottom:lat_top].data
    out['elevation'] = fh.variables['elevation'][lat_bottom:lat_top,lon_left:lon_right].data
    fh.close()
    out['grid'] = abs(out.ygrid[0]-out.ygrid[1])
    return out


#---------------------------------
def convert2higherResolutionGrid(Data, data, d_tag=['data']):#d_tag turn into a kwargs,args?
    """returns a higher resolution grid of Data
requires:
    Data    Bunch() or dict() with 'xgrid', 'ygrid', 'data'.
            Must have greater distance between gridpoints then _data_
    data    Bunch() or dict() with 'xgrid', 'ygrid', 'data'.
            Must have lesser distance between gridpoints then _Data_
    Grid    distance between gridpoints of _Data_ Bunch
    d_tag   Label(s) of the data to recalculate to new resolution.
            default = ['data']
ensures:
    Bunch() object with gridpoint system as _data_, with data from _Data_
"""
    Grid = Data.grid
    wk = Bunch()
    wk['x_grid'] = Data['xgrid'] + (Grid/2)
    wk['y_grid'] = Data['ygrid'] + (Grid/2)
    newgrid = deepcopy(Data)
    newgrid['x_grid'] = wk['x_grid'].searchsorted(data['xgrid'])
    newgrid['y_grid'] = wk['y_grid'].searchsorted(data['ygrid'])
    newgrid['xgrid'] = data['xgrid']
    newgrid['ygrid'] = data['ygrid']
    newgrid.grid = data.grid
    x,y = np.meshgrid(newgrid['x_grid'], newgrid['y_grid'])#double check
    for tag in d_tag:
        if tag not in ['grid','depth','xgrid','ygrid','time']:
            if len(Data[tag].shape)==2:
                newgrid[tag] = Data[tag][y,x]#double check
            elif len(Data[tag].shape)==3:
                newgrid[tag] = Data[tag][:,y,x]#double check
            else:
                raise IndexError('Number of dimensions of data are off for this algorythm')
    del newgrid['x_grid'],newgrid['y_grid']
    return newgrid


def calc_press_increase(B):
    """ calculates the pressure increase
relative to surface
"""
    return 1/(1 + B/10)


def calc_Q10(B, Q=2):#(pressure^temperature[metabolism--> Oxygen requirement]{needs standard value})-O2[availability]
    """https://en.wikipedia.org/wiki/Q10_(temperature_coefficient)
From the function:
<math>Q_{10}=\left( \frac{R_2}{R_1} \right )^{10 \,^\circ \mathrm{C}/(T_2-T_1) }</math>
Comes the reactivite increase:
<math>R_2 = R_1 ~Q_{10}^{(T_2-T_1)/10 \,^{\circ} \mathrm{C}}</math>
With the premisses here that:
R_1 = 1
Q_{10} = Q
T_1 = 0 celsius
T_2 = measured temperature
"""
    B['Q10'] = Q**(B.temperature/10)


def Celsius2Kelvin(Celsius):
    """turns degree Celsius temperatures into degree Kelvin
requires:    temperature in degrees Celsius(Int() or Float())
ensures:     temperature in degrees Kelvin
"""
    Celsius += 273


def createData(step_analysis=0,
               folder_c='D:\\PhD\\GIS-DBs\\copernico',
               folder_g='D:\\PhD\\GIS-DBs\\GEBCO',
               filename_g='GEBCO_2019.nc',
               filename_c1='global-reanalysis-phy-001-030-monthly_',
               filename_c2='global-reanalysis-bio-001-029-monthly_',
               filename_c3='global-reanalysis-bio-001-033-weekly_'):
    filename_gebco = os.path.join(folder_g,filename_g)
    filename1 = os.path.join(folder_c,filename_c1)
    filename2 = os.path.join(folder_c,filename_c2)
    filename3 = os.path.join(folder_c,filename_c3)
    if step_analysis==0:
        print('0. Collecting raw data')
        data = read_copernicus(filename1,filename2,filename3)
        with open("data.pickle",'wb') as handle:
            pickle.dump(data,handle)
        print("done")
    if step_analysis<=1:
        print('1. averaging time scale...', end="")
        with open("data.pickle",'rb') as handle:
            data = pickle.load(handle)
        averaged_copernicus(data)
        with open("dataa.pickle",'wb') as handle:
            pickle.dump(data,handle)
        print("done")
    if step_analysis<=2:
        print('2. expanding grid...', end="")
        with open("dataa.pickle",'rb') as handle:
            data = pickle.load(handle)
        depths = len(data[0].depth)+len(data[1].depth)
        for i,_ in enumerate(data[:2]):
            for key in data[i]:
                if key not in ['grid','depth','xgrid','ygrid','time']:
                    j = np.ones((depths,*data[i][key].shape[-2:]), dtype=np.float64)
                    for y,_ in enumerate(data[i].ygrid):
                        for x,_ in enumerate(data[i].xgrid):
                            if i==0:
                                wk = regression_points_calc(data[1].depth,np.c_[data[0].depth.tolist(),data[0][key][:,y,x].tolist()].tolist())
                            else:
                                wk = regression_points_calc(data[0].depth,np.c_[data[1].depth.tolist(),data[1][key][:,y,x].tolist()].tolist())
                            j[:,y,x] = [i[1] for i in wk]
                    data[i][key] = j
            data[i].depth = [i[0] for i in wk]

        standardize_grid(data)
        with open("dataaa.pickle",'wb') as handle:
            pickle.dump(data,handle)
        print("done")
    if step_analysis<=3:
        print('3. inserting all into 1 Bunch...', end="")
        with open("dataaa.pickle",'rb') as handle:
            data = pickle.load(handle)
        Celsius2Kelvin(data[0].temperature)
        Data = Bunch()
        Data['xgrid'] = data[0].xgrid
        Data['ygrid'] = data[0].ygrid
        Data['depth'] = np.array(data[0].depth)
        Data['grid'] = data[0].grid
        #print('gebco')
        b4 = read_gebco(filename_gebco,
                        decimal_Latitude_range="13.99,41.01",
                        decimal_Longitude_range="-35.01,-8.99")
        xx = b4.xgrid.searchsorted(Data.xgrid)
        yy = b4.ygrid.searchsorted(Data.ygrid)
        x,y = np.meshgrid(xx, yy)#double check
        del xx,yy
        d = np.expand_dims(Data.depth,axis=-1)
        d = np.expand_dims(d,axis=-1)
        b4.elevation = (d+np.stack([b4.elevation[y,x] for i in range(len(Data.depth))]))*(-1)#pelagicity#double check
        xx,yy = np.meshgrid(np.arange(len(Data.xgrid)),np.arange(len(Data.ygrid)))#double check
        j = np.ones((len(Data.depth),len(Data.ygrid),len(Data.xgrid)), dtype=np.float64)
        j[:,yy,xx] = calc_press_increase(d)#double check
        b4['pressure'] = j
        data.append(b4)
        for key in data[2]:
            if key not in ['grid','depth','xgrid','ygrid','time']:
                data[2][key] = np.stack([data[2][key] for i in range(len(Data.depth))])
        Data['coverages'] = np.stack([item for List in [[i[key] for key in i if key not in ['grid','depth','xgrid','ygrid','time']] for i in data] for item in List])
        del data
        with open("dataaaa.pickle",'wb') as handle:
            pickle.dump(Data,handle)
        print("done")
    if step_analysis<=4:
        print('4. retrieving final data...', end='')
        with open("dataaaa.pickle",'rb') as handle:
            Data = pickle.load(handle)
        print("done")
    return Data


def create_data_bunch(data,mask,test,train):
    """creates data bunch ready for SDM analysis
"""
    d = Bunch()
    d["coverages"] = data.coverages
    d["x_left_lower_corner"] = min(data.xgrid)
    d["y_left_lower_corner"] = min(data.ygrid)
    d["grid_size"] = data.grid
    d["Nx"] = len(data.xgrid)
    d["Ny"] = len(data.ygrid)
    d["Nz"] = len(data.depth)
    d["xgrid"] = data.xgrid
    d["ygrid"] = data.ygrid
    d["zgrid"] = data.depth
    d["mask"] = mask
    d["test"] = test
    d["train"] = train
    return d


def Load_cephalopods_macaronesia(step_analysis=4):
    Data = createData(step_analysis)
    water_mask = Data.coverages[-2]>=0


    dt = np.dtype([('species', np.unicode_,36),
                   ('dd long', np.float64),
                   ('dd lat', np.float64),
                   ('m depth', np.int16)])
    test = []
    train = []
    obs = read_sql_georeferenced_observations()
    for i in set(obs['species'].tolist()):
        print(np.count_nonzero(obs['species']==i),i)
        if np.count_nonzero(obs['species']==i)<10:
            continue
        oobs = obs[np.where(obs['species']==i)]
        n = round(len(oobs)*0.3)
        its,tr = np.array(WithoutReposition2(n,oobs))
        tst = np.array(negativeSampling(its,oobs))
        if len(test)==0:
            test = tst
            train = tr
        else:
            test = np.r_[test,tst]
            train = np.r_[train,tr]
    test = np.array(test,dtype=dt)
    train = np.array(train,dtype=dt)
    d = create_data_bunch(Data,water_mask,test,train)
    return d


if __name__=='__main__':
    print("hello")





