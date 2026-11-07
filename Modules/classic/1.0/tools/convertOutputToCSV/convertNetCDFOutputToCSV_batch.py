#!/space/hall1/sitestore/eccc/crd/ccrp/mib001/usr_conda/envs/py3-base/bin/python
import xarray as xr
import sys
import csv
import pandas as pd
import os, fnmatch
# This script takes in the CLASSIC output netCDF and writes it
# out as a csv file in a format suitable for AMBER. This script
# is intended only for site-level simulations.

# J. Melton. Jun 2019
# Read in the file info from the command line
input_dir = str(sys.argv[1])


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


for filename in find_files(input_dir, '*.nc'):
    if "rsFile" in filename:
        continue
    # Read in the file using xarray
    inNC = xr.open_dataset(filename)

    # Convert to a pandas dataframe
    df = inNC.to_dataframe()

    # Output to csv.
    df.to_csv(filename.replace('.nc', '.csv'))
