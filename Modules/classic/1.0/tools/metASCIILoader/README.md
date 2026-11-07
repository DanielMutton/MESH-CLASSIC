ASCII to NetCDF met file loader {#asciiMet}
========

The metASCIILoader program takes in ASCII met files and converts them to NetCDF format. It is located in tools/metASCIILoader.

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/metASCIILoader [file.MET] [longitude] [latitude] [containsLeaps]`

Where *longitude* and *latitude* are the desired real values for the local file and the *file.MET* is your local met file. *containsLeaps* is a logical switch. True if your input file contains leap years, i.e. years that contain Feb 29th, false if not. It is expecting either 'true' or 'false'.

The provided longitude and latitude values will be written as a property in the netcdf file.

See @ref forcingData for expected units and variables.

# Structure

The original met file looks like this:

        Minutes, Hour, Day, Year, Shortwave, Longwave, Precipitation, Temperature, Humidity, Wind, Pressure
        [...]
          6  0    1  1901     0.00   264.07    0.0000E+00     2.78   4.286E-03    3.02    99308.05
          6 30    1  1901     0.00   259.21    0.0000E+00     4.06   4.607E-03    1.33    99228.51
          7  0    1  1901     0.00   259.21    0.0000E+00     4.43   4.724E-03    1.29    99230.83
          7 30    1  1901     0.00   259.21    0.0000E+00     4.79   4.841E-03    1.26    99233.15
          8  0    1  1901    52.61   259.21    0.0000E+00     5.16   4.958E-03    1.23    99235.47
          8 30    1  1901   106.40   259.21    0.0000E+00     5.53   5.075E-03    1.19    99237.79
          9  0    1  1901   148.44   259.21    0.0000E+00     5.90   5.192E-03    1.16    99240.11
          9 30    1  1901   177.86   259.21    0.0000E+00     6.27   5.309E-03    1.13    99242.43
         10  0    1  1901   195.90   259.21    0.0000E+00     6.64   5.426E-03    1.09    99244.75
         10 30    1  1901   205.17   259.21    0.0000E+00     7.01   5.543E-03    1.06    99247.07
         11  0    1  1901   208.82   259.21    0.0000E+00     7.38   5.660E-03    1.03    99249.38
         11 30    1  1901   209.69   259.21    8.2154E-05     7.75   5.777E-03    0.99    99251.70
         12  0    1  1901   209.75   259.21    0.0000E+00     8.12   5.894E-03    0.96    99254.02
         12 30    1  1901   209.69   252.32    0.0000E+00     8.49   6.011E-03    0.93    99256.34
         [...]

It is important to ensure there are no extra columns or trailing white spaces after the pressure field.