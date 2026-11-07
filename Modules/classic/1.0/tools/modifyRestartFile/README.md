Modify netcdf model initialization files for whole fields or point locations {#modifyRS}
========

# Purpose

Aids in manipulating netcdf initialization files. The tool is intended for regional or global runs where it is less straight-forward to deal with a netcdf file.

1. See values of a given variable for a given grid cell from your restart file,

2. Print values of all variables for a given grid cell from your restart file to an ASCII file,

3. Change the values of a given variable for a given grid cell in your restart file (at multiple levels, if needed),

4. Change the values (and reduce the levels, if needed) of a given variable (at multiple levels, if needed) **everywhere** in your domain over land with domain-wide constant values


# Dependencies

The script uses bash shell scripting and python

# Execute

To run the program, use the following command:

`./modifyRestartFile.sh`

This will prompt with a usage message. The usage is:

 `./modifyRestartFile.sh [file to be modified]`

 After that follow prompts.
