#!/usr/bin/python

import xarray as xr
import numpy as np
#print("xarray version : ", xr.__version__)

# Open the CLASSIC output file that needs to be transformed to be on GMT:
hhoriginal = 'rmLeaf.nc'
variableName = 'rmLeaf'

# Open up a land mask.
landmask = '/raid/rd40/data/CTEM/OBSERVATION_DATASETS/CTEM_t63_landmask_trimmed.nc'

outfile = 'test.nc'

#======================

delt = 1800. # Model timestep for physics in seconds.
shortSteps = 86400. / delt # Number of physics timesteps in one day 


# Open using xarray.
hhdataorig = xr.open_dataset(hhoriginal)
mask = xr.open_dataset(landmask)

# Grab the longitudes
lons = hhdataorig.lon
lats = hhdataorig.lat

for lon in lons.data:
  for lat in lats.data:

      # Now adjust this longitude's values for the difference from Greenwich.
      # timeZone is the number of delt timesteps
      timeZone = (round(lon / 15. * 60.)) / (delt / 60.)
      
      if (lon >= 180.) or (lon < 0): # The longitudes run from 0 to 360 or from -180 to 180
          timeZone = shortSteps - timeZone
      else: # Because of how cshift works, this needs to be negative.
          timeZone = -timeZone
      
      timeZone = int(round(timeZone))
      if (timeZone != 0):      
          maskval = mask['landmask'].sel(lon=lon, lat=lat, method='nearest')
          if (maskval.data.flatten()[0] == 1):             

              # Now that we know how many timesteps we need to shift, we can apply the roll.
              newvals = hhdataorig[variableName].sel(lon=lon, lat=lat).roll(time=timeZone, roll_coords=False)
              # Replace the old values with the rolled values
              hhdataorig[variableName].loc[{'lon':lon, 'lat':lat}] = newvals

hhdataorig.to_netcdf(outfile)
