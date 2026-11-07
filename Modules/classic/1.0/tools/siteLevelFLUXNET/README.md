# FLUXNET site automation scripts

The tools in this directory were created with the intent of allowing easy use of CLASSIC for site-level runs on FLUXNET sites. The use case for each script will be covered briefly in this document.

## comparative\_plot\_generator.py

This script takes the output from site-level runs and plots them against observational data to give an idea of model performance over specific variables. The variables plotted are defined near the top of the script. When run from the command line, the most common format is
```
python3 comparative_plot_generator.py {path_to_CLASSIC_outputs_folder} -o {path_to_plot_output_folder}
```
Due to the large number of libraries used in this script, it is recommended that you run it from the CLASSIC singularity container.

## prep\_jobopts.sh

Takes the `template_job_options_file.txt` from configurationFiles and makes a copy for each FLUXNET site. Each site has unique attributes that must be specified in the file, such as start and end dates, as well as the location of meteorological files.

## process\_outputs.sh

Runs the `convertOutputToCSV` script on all netCDF outputs for the FLUXNET sites, then runs the `comparative_plot_generator.py` (so the user doesn't have to do so manually).

## run\_sites.sh

Sequentially runs the FLUXNET sites through CLASSIC. Note: only sites with a directory present in `inputFiles/FLUXNETsites` will be run. This script generally takes a large amount of time.

## table.py

Generates an informative table in .html format detailing site-level information. This can be useful for drawing insights on CLASSIC output.

## yaml_reader.py

Useful tool for reading site-level information from individual `siteinfo.yaml` files (contains FLUXNET metadata).
