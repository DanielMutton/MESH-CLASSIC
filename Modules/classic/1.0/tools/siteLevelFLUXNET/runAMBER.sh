# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}

# Setup configure.R with all the local path locations.
sed -i "/amber.gitrepo.path <-/s|\".*\"|\"$rootdir/tools/siteLevelFLUXNET/AMBER\"| ; /mod.csv.path/s|\".*\"|\"$rootdir/outputFiles/FLUXNETsites\"| ; /ref.csv.path <-/s|\".*\"|\"$rootdir/inputFiles/observationalDataFLUXNET\"| ; /outputDir <-/s|\".*\"|\"$rootdir/outputFiles/AMBER\"|" $rootdir/tools/siteLevelFLUXNET/AMBER/configure.R

Rscript $rootdir/tools/siteLevelFLUXNET/AMBER/configure.R 2>/dev/null

for file in outputFiles/AMBER/*.tex; do
  pandoc -s $file -o ${file%.tex}.html
done
