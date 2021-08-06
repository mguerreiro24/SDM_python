#https://www.earthinversion.com/utilities/Writing-NetCDF4-Data-using-Python/
import netCDF4
import numpy as np
import pickle

def nc_out(DEPTH):
    with open("dataaaa.pickle",'rb') as handle:
        Data = pickle.load(handle)
    variables = ["temp",
                 "sal",
                 "o2",
##                 "ed",
##                 "em",
##                 "mud",
##                 "musm",
##                 "mumm",
##                 "mld",
##                 "mlsm",
##                 "mlmm",
##                 "mlhmm",
                 "biomass",
                 "bd",
                 "pr"]
    for vi,v in enumerate(variables):
        ncfile = netCDF4.Dataset('data/new_{}{}.nc'.format(v,DEPTH),mode='w',format='NETCDF4_CLASSIC')
        #creating dims
        lat_dim = ncfile.createDimension('lat', len(Data.ygrid)) # latitude axis
        lon_dim = ncfile.createDimension('lon', len(Data.xgrid)) # longitude axis
##        depth_dim = ncfile.createDimension('depth', len(Data.depth)) # depth axis


        #creating atributes
        ncfile.title='Macaronesia data'
        ncfile.subtitle="environmental data Macaronesia"
        ncfile.description="average"


        #vars dims
        lat = ncfile.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
##        depth = ncfile.createVariable('depth', np.float32, ('depth',))
##        depth.units = 'm'
##        depth.long_name = 'depth in meters'
        #writting data
        nlats = len(Data.ygrid)
        nlons = len(Data.xgrid)
##        ndepths = len(Data.depth)
        lat[:] = Data.ygrid
        lon[:] = Data.xgrid
##        depth[:] = Data.depth
        #Creating Variables
        if vi==0:#temperature
##            temp = ncfile.createVariable('temp',np.float64,('depth','lat','lon'))
            temp = ncfile.createVariable('temperature',np.float64,('lat','lon'))
            temp.units = 'K'
            temp.standard_name = 'sea_water_potential_temperature'
        elif vi==1:#salinity
            temp = ncfile.createVariable('salinity',np.float64,('lat','lon'))
            temp.units = 'PSU'
            temp.standard = 'salinity'
        elif vi==2:#Oxygen
            temp = ncfile.createVariable('O2',np.float64,('lat','lon'))
            temp.units = 'mmol m-3'
            temp.standard = 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
        elif vi==3:#Biomass
            temp = ncfile.createVariable('Biomass',np.float64,('lat','lon'))
            temp.units = 'g m-2'
            temp.standard = 'mass_concentration_of_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==3:#epi depth
##            temp = ncfile.createVariable('epipelagic_depth',np.float64,('lat','lon'))
##            temp.units = 'm'
##            temp.standard = 'sea_water_epipelagic_layer_depth'
##        elif vi==4:
##            temp = ncfile.createVariable('epipelagic_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_epipelagic_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==5:
##            temp = ncfile.createVariable('mesopelagic_u_depth',np.float64,('lat','lon'))
##            temp.units = 'm'
##            temp.standard = 'sea_water_upper_mesopelagic_layer_depth'
##        elif vi==6:
##            temp = ncfile.createVariable('mesopelagic_u_static_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_upper_mesopelagic_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==7:
##            temp = ncfile.createVariable('mesopelagic_u_migratory_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_upper_migrant_mesopelagic_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==8:
##            temp = ncfile.createVariable('mesopelagic_l_depth',np.float64,('lat','lon'))
##            temp.units = 'm'
##            temp.standard = 'sea_water_lower_mesopelagic_layer_depth'
##        elif vi==9:
##            temp = ncfile.createVariable('mesopelagic_l_static_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_lower_mesopelagic_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==10:
##            temp = ncfile.createVariable('mesopelagic_l_migratory_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_lower_migrant_mesopelagic_micronekton_expressed_as_wet_weight_in_sea_water'
##        elif vi==11:
##            temp = ncfile.createVariable('mesopelagic_l_HighlyMigratory_mass',np.float64,('lat','lon'))
##            temp.units = 'g m-2'
##            temp.standard = 'mass_concentration_of_lower_highly_migrant_mesopelagic_micronekton_expressed_as_wet_weight_in_sea_water'
        elif vi==5:#13:
            temp = ncfile.createVariable('pressure',np.float64,('lat','lon'))
            temp.units = 'atm/surface'
            temp.standard = 'pressure increase by 10 m'
        elif vi==4:#12:
            temp = ncfile.createVariable('benthos_distance',np.float64,('lat','lon'))
            temp.units = 'm'
            temp.standard = 'distance to Seabed'

        #variables
        temp[:,:] = Data.coverages[vi,DEPTH]

        # close the Dataset.
        ncfile.close()#; print('Dataset is closed!')


if __name__=='__main__':
##    Data = netCDF4.Dataset('D:\PhD\GIS-DBs\copernico\macaronesia\global-reanalysis-bio-001-033-weekly_1614005680347.nc')
    for i in range(0,125):
        print(i)
        nc_out(i)
##    Dat = netCDF4.Dataset('data\newtemp0.nc')
##    with open("dataaaa.pickle",'rb') as handle:
##        Data = pickle.load(handle)
##    for i,e in enumerate(Data.depth):
##        print(i,e)
