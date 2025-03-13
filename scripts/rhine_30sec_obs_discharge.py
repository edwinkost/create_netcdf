

import pcraster as pcr
import netCDF4 as nc
import numpy as np

cloneMap   = 
ncFileName = 

pcr.setclone(cloneMap)
cloneMap = pcr.boolean(1.0)

# latitudes and longitudes
latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))[::-1]
longitudes = np.unique(pcr.pcr2numpy(pcr.xcoordinate(cloneMap), vos.MV))

# netcdf format and zlib setup 
format = 'NETCDF3_CLASSIC'
zlib   = False

# netCDF attributes
attributeDictionary = {}
attributeDictionary['institution'] = "Utrecht University"
attributeDictionary['title'      ] = "GRDC daily observation data"
attributeDictionary['description'] = "GRDC daily observation data mapped at 30sec resolution of PCR-GLOBWB drainage network."


# create the netcdf file
rootgrp = nc.Dataset(ncFileName, 'w', format = format)

# create dimensions - time is unlimited, others are fixed
rootgrp.createDimension('time', None)
rootgrp.createDimension('lat', len(latitudes))
rootgrp.createDimension('lon', len(longitudes))

# set the attribute information
for k, v in list(attributeDictionary.items()): setattr(rootgrp,k,v)

# create time variable and its unit and reference
date_time = rootgrp.createVariable('time', 'f4', ('time',))
date_time.standard_name = 'time'
date_time.long_name = 'Days since 1901-01-01'
date_time.units     = 'days since 1901-01-01'
date_time.calendar  = 'standard'

# create latitude and longitude variables
lat = rootgrp.createVariable('lat','f4',('lat',))
lat.long_name = 'latitude'
lat.units = 'degrees_north'
lat.standard_name = 'latitude'

lon = rootgrp.createVariable('lon','f4',('lon',))
lon.standard_name = 'longitude'
lon.long_name = 'longitude'
lon.units = 'degrees_east'

# set the values for latitude and longitudes
lat[:] = latitudes
lon[:] = longitudes

# create the variable "observed_discharge"
standardVarName = "obs_discharge"
shortVarName = standardVarName
longVarName  = "observed_discharge"
varUnits     = "m3.s-1"
var = rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) , fill_value = 1e20, zlib = zlib)
var.standard_name = standardVarName
var.long_name = longVarName
var.units = varUnits

# write to an empty netcdf file 
rootgrp.sync()
rootgrp.close()


# assign discharge values from a station
#
# step 1: read GRDC nc discharge file (note the grdc discharge file must contain only the dates 1 Jan 1990 to 31 Dec 2010, e.g. using: "cdo selyear,1990/2010 input.nc output.nc")
grdc_discharge_file = 
grdc_data        = nc.Dataset(grdc_discharge_file)
grdc_time_series = np.array(grdc_data.variables["runoff_mean"][:])
#
# step2: assign grdc_time_series to our netcdf file
rootgrp = nc.Dataset(ncFileName,  'a')
# - indices for latitude and longitude
i_lat = 30
i_lon = 30
rootgrp.variables[shortVarName][:,i_lat,i_lon] = grdc_time_series
rootgrp.sync()
rootgrp.close()


