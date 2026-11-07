# This script creates custom job options files for each fluxnet site in the inputFiles/FLUXNETsites directory
# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
inputs=$rootdir/inputFiles/FLUXNETsites
# We'll need the yaml reader to access the site_info.yaml for each site
yaml_reader=$script_location/yaml_reader.py
# Iterate through the site directories
for f in $inputs/*; do
  if [ -d $f ]; then
    # Get the name of the site (using its directory)
    current=${f##*/}
    echo "Creating job options file for $current..."
    # Get metadata about the site for later use
    start_year=$(python3 $yaml_reader $f/siteinfo.yaml start)
    end_year=$(python3 $yaml_reader $f/siteinfo.yaml end)
    ZRFH=$(python3 $yaml_reader $f/siteinfo.yaml ZRFH)
    ZRFM=$(python3 $yaml_reader $f/siteinfo.yaml ZRFM)
    jof=$f/job_options_file.txt
    yes | cp $rootdir/configurationFiles/template_job_options_file.txt $jof

    # Now we use sed (stream editor) to tailor the job options file with the info we gained
    # from site_info.yaml
    # Replace start year and end year
    sed -i "/readMetStartYear/s| = [0-9]* | = $start_year | ; /readMetEndYear/s| = [0-9]* | = $end_year | ; /metLoop/s| = [0-9]* | = 1 |" $jof

    # Replace all the meteorological files
    sed -i "/metFileFss/s|'.*'|'$f/metVar_sw.nc'| ; /metFileFdl/s|'.*'|'$f/metVar_lw.nc'| ; /metFilePre/s|'.*'|'$f/metVar_pr.nc'| ; /metFileTa/s|'.*'|'$f/metVar_ta.nc'| ; /metFileQa/s|'.*'|'$f/metVar_qa.nc'| ; /metFileUv/s|'.*'|'$f/metVar_wi.nc'| ; /metFilePres/s|'.*'|'$f/metVar_ap.nc'|" $jof

    # Use the proper CO2 file
    sed -i "/CO2File/s|'.*'|'inputFiles/TRENDY_CO2_1700_2018.nc'|" $jof

    # Replace the init and restart files
    sed -i "/init_file/s|'.*'|'$f/${current}_init.nc'| ; /rs_file_to_overwrite/s|'.*'|'$f/rsfile.nc'|" $jof

    # Turn off fires and land use, and competition
    sed -i "/lnduseon/s|=\s*\..*\.\s*,|= \.false\. ,| ; /fixedYearLUC/s|=\s*[-0-9]*\s*,|= -9999 ,| ; /dofire/s|=\s*\..*\.\s*,|= \.false\. ,| ; /PFTCompetition/s|=\s*\..*\.|= \.false\.|" $jof

    # Replace the IDISP, IZREF, ZRFH, and ZRFM
    sed -i "/IDISP/s|= [0-9\.]*\s*,|= 1 ,| ; /IZREF/s|= [0-9\.]*\s*,|= 1 ,| ; /ZRFH/s|= [0-9\.]*,|= $ZRFH ,| ; /ZRFM/s|= [0-9\.]*,|= $ZRFM ,|" $jof

    # Replace output_directory
    sed -i "/output_director/s|'.*'|'$rootdir/outputFiles/FLUXNETsites/$current'|" $jof

    mkdir -p outputFiles/FLUXNETsites/${current}
  fi

done
