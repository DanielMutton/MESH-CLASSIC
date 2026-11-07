ASCII to NetCDF CO2 file loader {#makeGHGfiles}
========

The ghgASCIILoader program takes in ASCII GHG files and converts them to NetCDF format.

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/ghgASCIILoader  [ghg name, e.g. CO2] [input file]`

Where ghg name is either 'CO2' or 'CH4' and the input file is an ASCII text file with year and GHG mole fraction. Expected units for both are ppmv.

# Structure

The text GHG file is expected to look like this:

          1700  276.59
          1701  276.62
          1702  276.65
          1703  276.67
          1704  276.70
          [...]
