Create netcdf model initialization files from depricated INI (and CTM) file formats or fortran namelists {#initTool}
========

The initFileConverter program takes in ASCII CLASS initialization files (.INI) and, if required, CTEM initialization files (.CTM) and converts them to the newer NetCDF format. **It is also possible to take in a fortran namelist format. This is the preferred method for site-level users.** 

Two example namelist files are included in this directory. example.nml is for a single tile whereas exampleTwoTile.nml has two tiles in the grid cell.

The expected format of the INI file is:
<pre><code>
      Testing Site,          !this line ignored in initFileConverter
  J. Melton                  !this line ignored in initFileConverter
  CRD                        !this line ignored in initFileConverter
     12.56     19.69     40.00     40.00     50.00   -1.0    1    1
   0.000   0.063   0.023   0.711   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
 200.000 125.000  85.000 100.000          30.000  40.000  30.000  30.000
   0.650   0.500   0.500   0.500           1.050   0.600   1.000   1.000
 100.000 100.000 100.000 100.000           5.000   5.000   5.000   5.000
   0.100   1.728   1.000
 4.0E-02 3.0E-01 2.0E+03 1.0E-05       1    **10**
      47.7      46.7      42.0
      24.7      27.2      31.8
       1.1       0.7       0.3
     20.85     26.60     27.82     13.74      0.00      0.00
     0.043     0.194     0.244     0.000     0.000     0.000     0.000
    0.0000    0.0000      0.00     0.000    0.0000     1.000
   0.100   0.100
   0.250   0.350
   3.750   4.100
         1       365         1       365 !this line ignored in initFileConverter
      2009      2010      2009      2010 !this line ignored in initFileConverter
</code></pre>

**Special Note: The converter expects a soil colour index in the file (highlighted above). Please see the CLASSIC manual for more information.**

* Additionally it is often best to make the snow in the canopy (SCAN), liquid in the canopy (RCAN) zero as they can often lead to instabilities when running the model at a new site for the first time (see CLASSIC manual for more on this).
* Any variables that are new to the initfile are set to default values (usually 0) if not in the INI or CTM files. They can, however, be set in the example namelist files. 
* The INI file read-in is not setup to read in > 4 CLASS PFTs or > 9 CTEM PFTs. It is desirable to move to the namelist technique for those cases, or edit the initFileConverter.f90 for your needs.

and the expected CTM format is:
<pre><code>
 CTEM INITIALIZATION FILE. CTEM's 9 PFTs ARE                                    !this line ignored in initFileConverter
     NDL     NDL     BDL     BDL     BDL    CROP    CROP   GRASS   GRASS        !this line ignored in initFileConverter
     EVG     DCD     EVG DCD-CLD DCD-DRY      C3      C4      C3      C4!Note 2 !this line ignored in initFileConverter
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
   1.000   0.000   0.000   0.000   1.000   0.429   0.571   0.629   0.371
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.27561E-03 8.57120E-02 1.27454E-02 4.19591E-04 4.67570E-02
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 4.09982E-03 4.96495E-02
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 9.47277E-02 1.89120E-01 4.18575E-03 0.00000E+00 0.00000E+00
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 8.07029E-02 4.43796E-02 2.66396E-03 9.20991E-02 1.21723E+00
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 3.91125E-02 1.32837E-01 1.43084E+00 8.01093E-02 6.45154E-01 0.00000E+00
 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.98504E-01 1.18263E-01 1.13012E+00 7.50368E-01 5.94591E+00 0.00000E+00
       4       4       4       4       1       2       2       1       1
       0       0       0       0       0       7       7       0       7
   0.000   0.000   0.000   0.009   0.012   0.119
   0.076   0.076   0.052   0.024   0.000   0.000
    0.50    !this line ignored in initFileConverter
    0.12    !this line ignored in initFileConverter
   1        !this line ignored in initFileConverter
 0.04661 0.20073 0.31897 0.45505 0.61088 0.72254 0.79003 0.83401                                  !this line ignored in initFileConverter
 0.00000 0.00000 0.00000 0.00000 0.00000 0.00015 0.00015 0.00015 0.00015 0.00015 0.00000 0.00000  !this line ignored in initFileConverter
</code></pre>

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/initFileConverter [file.INI or file.nml] [file.CTM]

The suffixes are case-insensitive. The namelist can have the suffix .txt or .nml.

The initialization file created will be placed in the same folder as the INI file and will have the same file prefix (file.nc).
