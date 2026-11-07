Convert output files to CSV {#convertToCSV}
========

# Purpose

It can be convenient in some instances to have the model output files in CSV format. This little script takes in the netCDF CLASSIC output file and converts to CSV.

# Important information

This script is designed for site-level simulations only, i.e. 1D output files.

# Dependencies

The script uses python with the xarray library (http://xarray.pydata.org/en/stable/). Available via http://xarray.pydata.org/en/stable/installing.html and can be installed using conda. Additional python libraries include sys,  csv, and pandas. All are generally included with an Anaconda installation.

# Execute

To run the program, use the following command:

./convertNetCDFOutputToCSV.py inFile.nc

and it will create a CSV format outut file with the same name (inFile.csv) in the same directory

# Example output

For a monthly sensible heat flux file:

        lat,lon,time,hfss
        -35.66,148.15,2001-01-31 00:00:00,76.36956653287922
        -35.66,148.15,2001-02-28 00:00:00,65.62731883147316
        -35.66,148.15,2001-03-31 00:00:00,52.27449460017114
