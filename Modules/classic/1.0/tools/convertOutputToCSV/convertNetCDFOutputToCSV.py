#!/usr/bin/env python3
import xarray as xr
import sys
import csv
import pandas as pd

# This script takes in the CLASSIC output netCDF and writes it 
# out as a csv file in a format suitable for AMBER. This script 
# is intended only for site-level simulations. 

# J. Melton. Jun 2019 

# Read in the file info from the command line
inNcFile = str(sys.argv[1])

# Read in the file using xarray
inNC = xr.open_dataset(inNcFile)

# Convert to a pandas dataframe 
df = inNC.to_dataframe()

# Output to csv.
df.to_csv(sys.argv[1].replace('.nc', '.csv'))

