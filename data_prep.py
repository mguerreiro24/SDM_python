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
def read_copernicus_phys(nc_file_phys):
    """
requires:
    file name of netCDF4 with physical variables from copernicus
    source: ftp://my.cmems-du.eu/Core/GLOBAL_REANALYSIS_PHY_001_030/global-reanalysis-phy-001-030-monthly
ensures:
    Bunch object
"""
    print("Physics, File: {}".format(nc_file_phys))
    #load data
    fh = Dataset(nc_file_phys)
    #extraction of data
    b1 = Bunch()
    b1['time'] = convert1950HHto1970ss(fh.variables['time'][:].data)
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
    source: ftp://my.cmems-du.eu/Core/GLOBAL_REANALYSIS_BIO_001_029/global-reanalysis-bio-001-029-monthly/
ensures:
    Bunch object
"""
    print('Primary Production, File: {}'.format(nc_file_PP))
    #load data
    fh = Dataset(nc_file_PP)
    #extraction of data
    b2 = Bunch()
    b2['time'] = convert1950HHto1970ss(fh.variables['time'][:].data)
    b2['depth'] = fh.variables['depth'][:].data
    b2['ygrid'] = fh.variables['latitude'][:].data
    b2['xgrid'] = fh.variables['longitude'][:].data
    b2['grid'] = abs(b2['ygrid'][0]-b2['ygrid'][1])
    
    b2['O2'] = fh.variables['o2'][:].data
#    b2['pH'] = fh.variables['ph'][:].data
#    b2['chla'] = fh.variables['chl'][:].data
#    b2['netPP'] = fh.variables['nppv'][:].data
#    b2['phytoplankyonC'] = fh.variables['phyc'][:].data
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


def read_copernicus_Zoo_n(nc_file_zoo1, nc_file_zoo2):#2 files
    """
requires:
    file names of netCDF4 with zooplankton from copernicus
        nc_file_zoo1    with variables: mnkc_epi,
                                        mnkc_umeso,
                                        mnkc_mumeso,
                                        mnkc_lmeso,
                                        mnkc_mlmeso,
                                        mnkc_hmlmeso,
                                        #zooc,
                                        #zeu,
                                        #npp
        nc_file_zoo2    with variables: pelagic_layer_depth,
                                        #T,
                                        #U,
                                        #V
    source1: ftp://my.cmems-du.eu/Core/GLOBAL_MULTIYEAR_BGC_001_033/cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i/
    source2: ftp://my.cmems-du.eu/Core/GLOBAL_MULTIYEAR_BGC_001_033/cmems_mod_glo_bgc_my_0.083deg-lmtl-Fphy_PT1D-i/
ensures:
    Bunch object
"""
    print('Zooplankton, File: {}'.format(nc_file_zoo1))
    #load data
    fh = Dataset(nc_file_zoo1)
    #extraction of data
    b3 = Bunch()
    b3['time'] = fh.variables['time'][:].data
    b3['ygrid'] = fh.variables['latitude'][:].data
    b3['xgrid'] = fh.variables['longitude'][:].data
    b3['grid'] = abs(b3['ygrid'][0]-b3['ygrid'][1])
    #epipelagic
    b3['mass_conc_epi'] = fh.variables['mnkc_epi'][:].data
    #upper mesopelagic
    b3['mass_conc_umeso'] = fh.variables['mnkc_umeso'][:].data
    b3['mass_conc_umeso_mig'] = fh.variables['mnkc_ummeso'][:].data
    #lower mesopelagic
    b3['mass_conc_lmeso'] = fh.variables['mnkc_lmeso'][:].data
    b3['mass_conc_lmeso_mig'] = fh.variables['mnkc_lmmeso'][:].data
    b3['mass_conc_lmeso_Hmig'] = fh.variables['mnkc_lhmmeso'][:].data
    fh.close()
    print('Zooplankton, File: {}'.format(nc_file_zoo2))
    fh = Dataset(nc_file_zoo2)
    #epipelagic
    b3['depth_epi'] = fh.variables['pelagic_layer_depth'][:,0].data
    #upper mesopelagic
    b3['depth_umeso'] = fh.variables['pelagic_layer_depth'][:,1].data
    #lower mesopelagic
    b3['depth_lmeso'] = fh.variables['pelagic_layer_depth'][:,2].data
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
        #insert sorted by time
        out[i]['time'] = np.concatenate([data[i][jj]['time'] for jj in range(len(data[i]))])
        sorting_indexes = out[i]['time'].argsort()
        out[i]['time'] = out[i]['time'][sorting_indexes]
        #concatenating on time dimension (dimension 0)
        for key in data[i][0]:
            if key not in ['grid','depth','xgrid','ygrid','time']:
                out[i][key] = np.concatenate([data[i][jj][key] for jj in range(len(data[i]))])
                out[i][key] = out[i][key][sorting_indexes]
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
    b2 = convert2higherResolutionGrid(data[1], data[0], d_tag=['O2'])#,
##                                                               'chla',
##                                                               'netPP',
##                                                               'pH',
##                                                               'phytoplankyonC'])
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
               ["genus", "species", longitude, latitude, UnixTime, depth]
ensures:
    np.array ready for Bunch
"""
    dt = np.dtype([('species', np.unicode_,36),
                   ('dd long', np.float64),
                   ('dd lat', np.float64),
                   ('s 1970', np.int32),
                   ('m depth', np.int16)])
    out = np.array([("_".join(x[0:2]), *x[2:]) for x in obs], dtype=dt)
    return out


def read_sql_georeferenced_observations(db_path = 'D:\\PhD\\articles\\article DB',
                                        db = 'samples.db',
                                        test_Lat_range='14,40',
                                        test_Lon_range='-35,-9'):
    """from SQL database file, collect georeferenced entries
requires:
    db_path           path of folder where db file is
                        default: 'D:\\PhD\\articles\\article DB'
    db                db file name
                        default: 'samples.db'
    test_Lat_range    [str] csv of min and max latitude for data
                      that will be used for testing and mapping
                        default: '14,40'
    test_Lon_range    [str] csv of min and max longitude for data
                      that will be used for testing and mapping
                        default: '-35,-9'
ensures: Numpy array with
    'species' as 'genus_species',
    'dd long' float(),
    'dd lat' float() ,
    's 1970' as int() and 
    'depth' as int()
"""
    conn = sqlite3.connect(os.path.join(db_path,db))
    cursor = conn.cursor()
    cursor.execute('''
SELECT DISTINCT genus,sp,Longitude,Latitude,UnixTime,depth_min,depth_Max
FROM sample inner join species using (id_s) inner join net using (id_n,id_st,id_c) inner join station using (id_st,id_c) inner join cruise using (id_c)
WHERE id_c=1 OR id_c>3
ORDER BY genus,sp''')# OR id_c>3
    Lobs = cursor.fetchall()
    cursor.close()
    conn.close()
    Lobs = [[*i[:-2],round((i[-2]+i[-1])/2)] for i in [i for i in Lobs if type(i[-2])==int and type(i[-1])==int]]
    latm,latM = [float(t) for t in test_Lat_range.split(',')]
    lonm,lonM = [float(t) for t in test_Lon_range.split(',')]
    Lobs = [i for i in Lobs if lonm<=i[2]<=lonM and latm<=i[3]<=latM]
    return read_georeferenced_observations(Lobs)


def train_test():
    """preparing observation data for train and test of model
requires:
ensures:test,train datasets [Bunch()]
"""
    dt = np.dtype([('species', np.unicode_,36),
                   ('dd long', np.float64),
                   ('dd lat', np.float64),
                   ('s 1970', np.int32),
                   ('m depth', np.int16)])
    test = []
    train = []
    #collecting data(sensitive point)
    obs = read_sql_georeferenced_observations(test_Lat_range='-2,49',test_Lon_range='-45,-5')
    #point of selecting ratio and sampling for the 2 categories
    for i in set(obs['species'].tolist()):
        #print(np.count_nonzero(obs['species']==i),i)
        if np.count_nonzero(obs['species']==i)<9:
            continue
        oobs = obs[np.where(obs['species']==i)]
        n = round(len(oobs)*0.7)#percentage for training
        its,tr = WithoutReposition2(n,oobs)
        tst = negativeSampling(its,oobs)
        tr = np.array(tr)
        tst = np.array(tst)
        if len(test)==0:
            test = tst
            train = tr
        else:
            test = np.r_[test,tst]
            train = np.r_[train,tr]
    test = np.array(test,dtype=dt)
    train = np.array(train,dtype=dt)
    return test,train


def train_test_PA(db_path = 'D:\\PhD\\articles\\article DB',
                  db = 'samples.db'):
    """preparing [presence-]absence data for train and test of model
requires:
ensures:test_A,train_A -> absence datasets [Bunch()]
"""
    dt = np.dtype([('species', np.unicode_,36),
                   ('dd long', np.float64),
                   ('dd lat', np.float64),
                   ('s 1970', np.int32),
                   ('m depth', np.int16)])
    conn = sqlite3.connect(os.path.join(db_path,db))
    cursor = conn.cursor()
    cursor.execute('''
SELECT DISTINCT Longitude,Latitude,UnixTime,depth_min,depth_Max
FROM net inner join station using (id_st,id_c) inner join cruise using (id_c)
WHERE id_c=1 OR id_c>3
''')# OR id_c>3
    Lobs = cursor.fetchall()
    cursor.close()
    conn.close()
    Lobs = [[*i[:-2],round((i[-2]+i[-1])/2)] for i in [i for i in Lobs if type(i[-2])==int and type(i[-1])==int]]

    test,train = train_test()
    tt = np.r_[test,train]
    tt = np.array(sorted(tt,key=lambda x:x[0]))
    absent = []
    for specie in set(tt['species'].tolist()):
        l = tt[np.where(tt['species']==specie)]
        ll = list(zip(l['dd long'],l['dd lat'],l['s 1970'],l['m depth']))
        for o in Lobs:
            if o not in ll:
                absent.append([*specie.split('_'),*o])
    l_test = [1]*len(test)
    l_train = [1]*len(train)
    obs = read_georeferenced_observations(absent)
    #point of selecting ratio and sampling for the 2 categories
    for i in set(obs['species'].tolist()):
        oobs = obs[np.where(obs['species']==i)]
        n = round(len(oobs)*0.7)#percentage for training
        its,tr = WithoutReposition2(n,oobs)
        tst = negativeSampling(its,oobs)
        tr = np.array(tr)
        tst = np.array(tst)

        test = np.r_[test,tst]
        train = np.r_[train,tr]
    test = np.array(test,dtype=dt)
    train = np.array(train,dtype=dt)
    la_test = [0]*(len(test) - len(l_test))
    la_train = [0]*(len(train) - len(l_train))
    label_test = np.array(l_test + la_test)
    label_train = np.array(l_train + la_train)
    return (test,label_test),(train,label_train)


def create_community_bunch(species_name, train, test):
    """Create a bunch with information about a community

    This will use the test/train record arrays to extract the
    data specific to the given species name.
requires:
    species_name   list(str()): list(str(genus_species))
    train          np.array dtype=[('species', '<U36'), ('dd long', '<f8'), ('dd lat', '<f8'), ('s 1970', '<i4'), ('m depth', '<i2')])
    test           np.array dtype=[('species', '<U36'), ('dd long', '<f8'), ('dd lat', '<f8'), ('s 1970', '<i4'), ('m depth', '<i2')])
ensures:
    list(Bunch(species))
    coverages      
                       Copernicus - Physics,
                       Copernicus - Biochemical[O2],
                       Copernicus - Zooplankton,
                       GEBCO - distance from seaBed,
                       GEBCO - pressure increase by 10 m
    """
    f1 = "mercatorglorys12v1_gl12_mean_"
    f2 = "mercatorfreebiorys2v4_global_mean_"
    f3 = "global-reanalysis-bio-001-033-weekly_"
    folder_c = 'D:\\PhD\\GIS-DBs\\copernico\\greater range'
##    folder_c = 'D:\\PhD\\GIS-DBs\\copernico\\greater_range'
    file1 = os.path.join(folder_c,f1)
    file2 = os.path.join(folder_c,f2)
    file3 = os.path.join(folder_c,f3)
    folder_g='D:\\PhD\\GIS-DBs\\GEBCO'
    filename_g='GEBCO_2019.nc'
    fileg = os.path.join(folder_g,filename_g)
    list_B = [Bunch(name=' '.join(species.split("_"))) for species in species_name]
    
    #species_name = species_name.encode('ascii')
    points = dict(test=test, train=train)
    for si,species in enumerate(species_name):
        for label, pts in points.items():
            # choose points associated with each desired species (indexed to 'list_B')
            pts = pts[pts['species'] == species]
            list_B[si]['pts_%s' % label] = pts

    # determine coverage values for each of the training & testing points
    w_data = [{label:[0,0,0,0,0] for label in points.keys()} for sp in species_name]
    data = [0,0,0,0,0]
    for i,ntcd in enumerate((file1,file2,file3)):
        if i==0:
            wfunc = read_copernicus_phys
        elif i==1:
            wfunc = read_copernicus_PP
        elif i==2:
            wfunc = read_copernicus_Zoo
        path = '\\'.join(ntcd.split('\\')[:-1])
        #prepping up path of files to load (coverage files)
        if path=='':
            list_files = os.listdir()
        else:
            list_files = os.listdir(path)

        #getting the files names
        ncfileregex = [j for j in list_files if j.startswith(ntcd.split('\\')[-1])]
        #preparing to load data variables
        data[i] = [0]*len(ncfileregex)
        for label,_ in points.items():
            for si,_ in enumerate(species_name):
                w_data[si][label][i] = [0]*len(ncfileregex)

        #
        for ii,file in enumerate(ncfileregex):
            data[i][ii] = wfunc(os.path.join(path,file))#data loaded here
            for label,_ in points.items():
                for si,_ in enumerate(species_name):
                    w_data[si][label][i][ii] = Bunch()
                    w_data[si][label][i][ii]['time'] = data[i][ii].time#time info

            #getting positional data on loaded dataset
            grid = data[i][ii].grid
            xgrid = data[i][ii].xgrid + grid/2
            ygrid = data[i][ii].ygrid + grid/2
            #getting the points from the raw data of copernicus
            for label,_ in points.items():
                for si,_ in enumerate(species_name):
                    pts = list_B[si]['pts_%s' % label]
                    ix = np.searchsorted(xgrid, pts['dd long'])
                    iy = np.searchsorted(ygrid, pts['dd lat'])
                    if 'depth' in data[i][ii]:#getting positional data on loaded dataset
                        zgrid = data[i][ii].depth
                        iz = np.searchsorted(zgrid, pts['m depth'])
                        for key in data[i][ii]:
                            if key not in ['grid','depth','xgrid','ygrid','time']:
                                w_data[si][label][i][ii][key] = data[i][ii][key][:,iz,iy,ix]
                    else:
                        for key in data[i][ii]:
                            if key not in ['grid','depth','xgrid','ygrid','time']:
                                w_data[si][label][i][ii][key] = data[i][ii][key][:,iy,ix]
            data[i][ii] = 0
        #
        for label,_ in points.items():
            for si,_ in enumerate(species_name):
                pts = list_B[si]['pts_%s' % label]
                work_set = sorted(w_data[si][label][i], key=lambda x:x.time[0])
                time = np.array([t.time for t in work_set])

                #get closest date|apply to depth
                iw = np.argmin((np.abs(time-pts['s 1970'])),axis=0)
                w_data[si][label][i] = Bunch()
                for key in work_set[0]:
                    if key not in ['grid','depth','xgrid','ygrid','time']:
                        w_data[si][label][i][key] = np.stack([work_set[hh][key][0][ih] for ih,hh in enumerate(iw)], axis=0)
        del work_set

    #2 last entries are dependent of the GEBCO database of bathimetries
    setD = read_gebco(fileg,"-2.0,50.0","-36,-0.1")
    grid = setD.grid
    xgrid = setD.xgrid + grid/2
    ygrid = setD.ygrid + grid/2

    for label,_ in points.items():#finishing all data sets into output
        for si,_ in enumerate(species_name):
            # build last ones dependent on GEBCO
            pts = list_B[si]['pts_%s' % label]
            ix = np.searchsorted(xgrid, pts['dd long'])
            iy = np.searchsorted(ygrid, pts['dd lat'])
            #distance benthos
            w_data[si][label][3] = (pts['m depth']+setD['elevation'][iy, ix])*-1
            #pressure increase
            w_data[si][label][4] = calc_press_increase(pts['m depth'])
            #temperature in kelvin transformation
            Celsius2Kelvin(w_data[si][label][0]['temperature'])
            s = []#data coverages for bunch
            for h,_ in enumerate(w_data[si][label][:3]):#first 3 sets
                for key in w_data[si][label][h]:
                    if key not in ['grid','depth','xgrid','ygrid','time']:
                        s.append(w_data[si][label][h][key])
            for h in w_data[si][label][3:]:#last ones
                s.append(h)
                
            list_B[si]['cov_%s' % label] = np.stack(s).T
    return list_B


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
    print('Bathymetry, File: {}'.format(filename))
    out = Bunch()
    out['grid'] = abs(fh.variables['lon'][0].data-fh.variables['lon'][1].data)
    yy = fh.variables['lon'][:].data# + out.grid/2.0
    xx = fh.variables['lat'][:].data# + out.grid/2.0
    lon_left = yy.searchsorted(float(decimal_Longitude_range.split(',')[0]))
    lon_right = yy.searchsorted(float(decimal_Longitude_range.split(',')[1]))
    lat_bottom = xx.searchsorted(float(decimal_Latitude_range.split(',')[0]))
    lat_top = xx.searchsorted(float(decimal_Latitude_range.split(',')[1]))
    del xx,yy
    out['xgrid'] = fh.variables['lon'][lon_left:lon_right].data
    out['ygrid'] = fh.variables['lat'][lat_bottom:lat_top].data
    out['elevation'] = fh.variables['elevation'][lat_bottom:lat_top,lon_left:lon_right].data
    fh.close()
    
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


def convert1950HHto1970ss(hours):
    """Converts a hours since 1950-01-01 00:00:00 time stamp into a
seconds since 1970-1-1 00:00:00 time stamp
requires:number of hours since 1950-01-01 00:00:00
ensures:seconds since 1970-1-1 00:00:00
"""
    return (hours*60**2)-631152000


#------------------------------------------------------------------------------
def createData(step_analysis=0,
               folder_c='D:\\PhD\\GIS-DBs\\copernico\\macaronesia',
               folder_g='D:\\PhD\\GIS-DBs\\GEBCO',
               filename_g='GEBCO_2019.nc',
               filename_c1='global-reanalysis-phy-001-030-monthly_',
               filename_c2='global-reanalysis-bio-001-029-monthly_',
               filename_c3='global-reanalysis-bio-001-033-weekly_'):#data for map and test
    def Biomass(data):
        """for createData function, processing zooplankton data step
    """
        shape = data[0].temperature.shape
        biom = np.zeros(shape)
        
        data[2]['biomass'] = biom
        for d in range(0,shape[0]):
            for y in range(0,shape[1]):
                for x in range(0,shape[2]):
                    #epipelagic
                    if data[2].depth_epi[y,x]>=data[0].depth[d]:
                        try:
                            if data[2].mass_conc_lmeso_Hmig[y,x]>10**4 or data[2].mass_conc_umeso_mig[y,x]>10**4 or data[2].mass_conc_epi[y,x]>10**4:
                                raise OverflowError
                            data[2].biomass[d,y,x] = data[2].mass_conc_epi[y,x] + data[2].mass_conc_lmeso_Hmig[y,x] + data[2].mass_conc_umeso_mig[y,x]
                        except OverflowError:
                            try:
                                if data[2].mass_conc_umeso_mig[y,x]>10**4 or data[2].mass_conc_epi[y,x]>10**4:
                                    raise OverflowError
                                data[2].biomass[d,y,x] = data[2].mass_conc_epi[y,x] + data[2].mass_conc_umeso_mig[y,x]
                            except OverflowError:
                                if data[2].mass_conc_epi[y,x]<10**4:
                                    data[2].biomass[d,y,x] = data[2].mass_conc_epi[y,x]
                    #upper mesopelagic
                    elif data[2].depth_umeso[y,x]>=data[0].depth[d]:
                        try:
                            if data[2].mass_conc_umeso[y,x]>10**4 or data[2].mass_conc_lmeso_mig[y,x]>10**4:
                                raise OverflowError
                            data[2].biomass[d,y,x] = data[2].mass_conc_umeso[y,x] + data[2].mass_conc_lmeso_mig[y,x]
                        except OverflowError:
                            if data[2].mass_conc_umeso[y,x]<10**4:
                                data[2].biomass[d,y,x] = data[2].mass_conc_umeso[y,x]
                    #lower mesopelagic
                    elif data[2].depth_lmeso[y,x]>=data[0].depth[d]:
                        if data[2].mass_conc_lmeso[y,x]<10**4:
                            data[2].biomass[d,y,x] = data[2].mass_conc_lmeso[y,x]
        del data[2]['depth_epi'],data[2]['mass_conc_epi'],data[2]['depth_umeso'],data[2]['mass_conc_umeso'],data[2]['mass_conc_umeso_mig'],data[2]['depth_lmeso'],data[2]['mass_conc_lmeso'],data[2]['mass_conc_lmeso_mig'],data[2]['mass_conc_lmeso_Hmig']
        return data

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
        print('2. expanding grid...')
        with open("dataa.pickle",'rb') as handle:
            data = pickle.load(handle)
        depths = len(data[0].depth)+len(data[1].depth)
        print('higher resolution depth')
        for i,_ in enumerate(data[:2]):
            for key in data[i]:
                if key not in ['grid','depth','xgrid','ygrid','time']:
                    print(key)
                    j = np.ones((depths,*data[i][key].shape[-2:]), dtype=np.float64)
                    for y,_ in enumerate(data[i].ygrid):
##                        print(y)
                        for x,_ in enumerate(data[i].xgrid):
##                            print(y,x)
                            if i==0:
                                wk = regression_points_calc(data[1].depth,np.c_[data[0].depth.tolist(),data[0][key][:,y,x].tolist()].tolist())
                            else:
                                wk = regression_points_calc(data[0].depth,np.c_[data[1].depth.tolist(),data[1][key][:,y,x].tolist()].tolist())
                            j[:,y,x] = [i[1] for i in wk]
                    data[i][key] = j
            data[i].depth = [i[0] for i in wk]

        print('higher resolution grid')
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
##        for key in data[2]:#expanding to depth
##            if key not in ['grid','depth','xgrid','ygrid','time']:
##                data[2][key] = np.stack([data[2][key] for i in range(len(Data.depth))])
        data = Biomass(data)
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
    test,train = train_test()
    d = create_data_bunch(Data,water_mask,test,train)
    return d


def biomass(filename):
    """receievs old occurrences data, appends new field of biomass
requires:
    filename    name of the file[str].
                    _*tsv format*_ with headers:
                    -Longitude,
                    -Latitude,
                    -Bathymetry,
                    -temperature,
                    -salinity,
                    -O2,
                    -pressure,
                    -benthos_distance,
                    -epipelagic_depth,
                    -epipelagic_mass,
                    -mesopelagic_u_depth,
                    -mesopelagic_u_static_mass,
                    -mesopelagic_u_migratory_mass,
                    -mesopelagic_l_depth,
                    -mesopelagic_l_static_mass,
                    -mesopelagic_l_migratory_mass,
                    -mesopelagic_l_HighlyMigratory_mass
ensures:
    new file    filename+_new+.tsv
    In the above mentioned, last column has the biomass at that
        depth according to the values obtained in copernicus, by:
            -Occurrences shallower than the epipelagic depth will
                have values equal to the sum of the biomass found in
                the epipelagic, plus the migratory upper mesopelagic
                and the highly migratory lower mesopelagic;
            -Occurrences shallower than the upper mesopelagic depth
                but deeper than the epipelagic will have values
                equal to the sum of the biomass found in the
                upper mesopelagic and the migratory lower
                mesopelagic;
            -Occurrences shallower than the lower mesopelagic depth
                but deeper than the upper mesopelagic will have
                values equal to the sum of the biomass found in the
                lower mesopelagic;
            -Deeper occurrences will have 0.
"""
    with open(filename,'r') as handle:
        data = [i[:-1].split("\t") for i in handle.readlines()]

    #processing here!
    bm = 0
    for i,e in enumerate(data):
        if i==0:
            data[i].append('Biomass')
        else:
            if e[8]=='':
                data[i].append('')
                continue
            if float(e[2])<=float(e[8]):
                try:
                    bm = float(e[9]) + float(e[12]) + float(e[16])
                except:
                    bm = ''
            elif float(e[2])<=float(e[10]):
                try:
                    bm = float(e[11]) + float(e[15])
                except:
                    bm = ''
            elif float(e[2])<=float(e[13]):
                try:
                    bm = e[14]
                except:
                    bm = ''
            else:
                bm = 0.0
            data[i].append(str(bm))
    with open('{}_new.{}'.format(*filename.split('.')),'w') as handle:
        handle.write("{}\n".format('\n'.join(list(map('\t'.join,data)))))
if __name__=='__main__':
    print("hello")
    biomass('data\\Hreversa.tsv')
    #d = Load_cephalopods_macaronesia()
    #obs = read_sql_georeferenced_observations(test_Lat_range='-2,49',test_Lon_range='-45,-5')
##    test,train = train_test_PA()
##    species_name = 'Abraliopsis_atlantica'
    #len(test[0][np.where(test[1]==1)])
##    d = createData(3)
    




