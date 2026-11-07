Modify global scale half-hourly output files to GMT {#shiftToGMT}
========

# Purpose

The half-hourly files use the timestamps at each cell's location. Some users require the files to be all referenced to GMT. This script adjusts the model outputs at each longitude to be GMT. This means a shift and for the first and last days of the time series some data is taken from the previous/next days.

# Important information

For this to work properly your run should be outputting half-hourly outputs for the full year, i.e. day 1 to 365 as it does a shift for the data that is up to 12 hours (longitude dependent).

# Dependencies

The script uses python with the xarray library (http://xarray.pydata.org/en/stable/) and numpy. Available via http://xarray.pydata.org/en/stable/installing.html and can be installed using conda.

# Execute

To run the program, edit the paths to the original half-hourly output file (hhoriginal), the CLASSIC land mask (landmask; this could removed from the script if you don't have a landmask, it just speeds the script up), the variable name (variableName) and the output file (outfile)

