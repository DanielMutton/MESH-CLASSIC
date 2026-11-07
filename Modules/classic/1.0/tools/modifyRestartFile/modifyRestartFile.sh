#!/bin/bash

# This script modifies rsFile (needed as argument) for a given grid cell and for given variables

# IT REQUIRES THE USE OF ANACONDA PYTHON

# Vivek Arora.


if [ $# -lt 1 ]; then
  echo " "
  echo "Usage: "
  echo "$0    rs_file_name.nc "
  echo " "
  echo "This script allows to ... "
  echo " "
  echo "a) see values of a given variable for a given grid cell from your restart file,"
  echo " "
  echo "b) print values of all variables for a given grid cell from your restart file to an ASCII file,"
  echo " "
  echo -e "c) change the values of a given variable for a given grid cell in your restart file (at\n   multiple levels, if needed), and"
  echo " "
  echo -e "d) change the values (and reduce the levels, if needed) of a given variable (at multiple levels,\n   if needed) EVERYWHERE in your domain over land with domain-wide constant values."
  echo " "
  echo " "
  exit
fi

#----------------------------------------------------------------------------------------------
# TYPICALLY USER INPUT NOT REQUIRED BELOW THIS LINE

clear

orig_rs=$1


# if string contains forward slash that means file names points to a non-local file
# in some folder. in this case get the filename of rsFile ignoring the whole path
# for later use

if [[ "$1" == *\/* ]]; then  # if path contains /
   rsfile_name_without_path=`echo $1 | rev | awk -F'/' '{print $1}' | rev`
else # if it's a local file in his folder
   rsfile_name_without_path=$1
fi




mkdir -p temp_dir
rm -fr temp_dir/*
cp $1 temp_dir/. # cp restrat file to our temp_dir so that we cna work on it
cd temp_dir

d1=`date +%Y_%m_%d_%H_%M_%S`
output_rs="rsFile_${d1}.nc"
output_rs3="rsFile_${d1}.nc3"
output_tmp="rsFile_tmp.nc"

# cp orig+rs to temp_dir


if [[ "$1" == *\/* ]]; then  # if path contains /
   rsfile_name_without_path=`echo $1 | rev | awk -F'/' '{print $1}' | rev`
   cp $orig_rs .
else # if it's a local file in this folder
   cp ../$orig_rs .
fi


# remove history
ncatted -hO -a history,global,d,, $orig_rs # remove history


# spit out the names of all variables in the restart file

vars_in_rs_file=`ncdump -h ${rsfile_name_without_path} | grep 'int\|float\|double' | awk '{print $2}' | awk -F'(' '{print $1}'`
ncdump -h ${rsfile_name_without_path} | grep 'int\|float\|double\|long_name' | paste - - -d '|' | grep -v "history"  | grep -v "owner of" > vars.txt


# remove dimensions from variable list
cat vars.txt | grep -v "double lat\|double lon\|double tile\|double layer\|double icp1\|double iccp1\|double ic\|double icc\|double months\|double slope\|int nmtest" | sort -k2  > vars2.txt



touch -f vars3.txt


print_var_names(){
  echo " "
  echo "Following variables are in the restart file for CLASS-CTEM other than the dimensional"
  echo "variables like lat, lon, tile, month, icp, icc, etc."
  echo " "

  if [ $first_time -eq 1 ]; then
    sleep 2
  else
    sleep 0
  fi

  count=1
  while read line; do
    part1=`echo $line | awk -F'|' '{print $1}' | sed 's/(/ (/g'`
    part2=`echo $line | awk -F'|' '{print $2}' | sed  's/^[ \t]*//;s/[ \t]*$//'` # sed removes leading and trailing spaces
    long_name=`echo $part2 | awk -F'"' '{print $2}' | awk -F'"' '{print $1}'`

    #echo "$count | $line "| sed 's/(/ (/g'
    echo "$count | $part1 | $long_name" | awk -F'|' '{printf "%3d | %-45s | %s\n", $1, $2,$3}'
    echo "$count | $part1 | $long_name" | awk -F'|' '{printf "%3d | %-45s | %s\n", $1, $2,$3}' >> vars3.txt
    count=`expr $count + 1`
  done < vars2.txt
}

choose_variable () {
   echo "---------------------------------------------------------------------------------- "
   echo -e "\e[1;7mSelect a variable from the above list to see its values for a given grid cell ...\e[0m"
   echo " "
   read  -p "Type here variable name you selected :" varselected

   echo " "
   echo "You selected:" $varselected
   echo " "

   # check if typed variable exists

   exist=0
   for v in $vars_in_rs_file; do
     if [ "x${varselected}" = "x${v}" ]; then
       exist=1
     fi
   done

   if [ $exist -eq 0 ]; then
      echo " "
      echo " "
      echo -e "\e[1;7mDid you mistype? Your selected variable doesn't exist. Note that variables \e[0m"
      echo -e "\e[1;7mare case-senitive. Please try again or press Ctrl-C to exit. \e[0m"
      echo " "
      echo " "
      choose_variable
   fi
}


what_would_u_like () {
   echo " "
   echo "---------------------------------------------------------------------------------- "
   echo " "
   echo -e "\e[1;7mWould you like to \n
1) see [and modify] values of a given variable for a given grid cell [or everyhere on land in your domain] or
2) print out values of all variables for a given grid cell\n\e[0m"
   echo " "
   read  -p "Choose 1 or 2 :" whatwouldulike

   if [ "x$whatwouldulike" = "x1" ]; then
     uselss=1
   elif [ "x$whatwouldulike" = "x2" ]; then
     uselss=1
   else
     echo " "
     echo "Did you mistype? Only 1 and 2 are valid options."
     echo " "
     what_would_u_like
   fi

}

what_would_u_like


if [ "x$whatwouldulike" = "x1" ]; then # user selected to see [and modify] values for a given grid cell
  first_time=1
  print_var_names
  choose_variable
fi


# find lat and lon range and min and max indices

min_lat=`ncdump -v lat ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lat = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print $1}'`
max_lat=`ncdump -v lat ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lat = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print $NF}'`
no_of_lats=`ncdump -v lat ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lat = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print NF}'`

min_lon=`ncdump -v lon ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lon = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print $1}'`
max_lon=`ncdump -v lon ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lon = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print $NF}'`
no_of_lons=`ncdump -v lon ${rsfile_name_without_path} | sed -e '1,/data:/ d' | sed 's/lat = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -F',' '{print NF}'`




get_lat_lon_from_user () {

  echo " "
  echo "---------------------------------------------------------------------------------- "
  echo -e "\e[1;7mNow enter latitude/longitude (as floating point numbers) of this grid cell or its lat/lon indices (as integers) ...\e[0m"
  echo " "
  echo "e.g. 10.50/286.75 (for actual latitude/longitude)"
  echo " "
  echo "or 48/102 (for lat/lon indices) "
  echo " "
  echo "NOTE THAT in your restart file : Latitudes go from $min_lat [1] to $max_lat [$no_of_lats], and "
  echo "                                 Longitudes go from $min_lon [1] to $max_lon [$no_of_lons] "
  echo " "


  echo " "
  read  -p "Type here latitude/longitude (as floating point numbers) or lat/lon indices (as integers) : " latloncombo

  lat1=`echo $latloncombo | awk -F'/' '{print $1}'`
  lon1=`echo $latloncombo | awk -F'/' '{print $2}'`



  # spit out a python script that we will use
  # this script find the nearest lat/lon and its index

  echo "
from netCDF4 import Dataset
import numpy as np
import sys


def geo_idx(dd, dd_array):
   \"\"\"
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    \"\"\"
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx


args=sys.argv

if len(sys.argv) != 4:
   print \"Need 3 arguments\"
   print \"Usage: anaconda python script.py 1)filename 2)longitude 3)latitude\"
   quit()


infile=args[1]
in_lon = float(args[2])
in_lat = float(args[3])

nc_f = infile
nci = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                             # and create an instance of the ncCDF4 class

#nci = netCDF4.Dataset(infile)
lats = nci.variables['lat'][:]
lons = nci.variables['lon'][:]

lat_idx = geo_idx(in_lat, lats)
lon_idx = geo_idx(in_lon, lons)


print lons[lon_idx], lon_idx
print lats[lat_idx], lat_idx
" > script.py



  echo $lat1, $lon1

  if [[ "$lat1" =~ ^[0-9]+$ ]]; then # if lat1 is an integer
     lat_index=$lat1
     lat_integer=1
     lat=`ncdump -v lat $rsfile_name_without_path | sed -e '1,/data:/ d' | sed 's/lat = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -v L=$lat_index -F',' '{print $L}'`
  else
     # call our python script to find nearest lat and it's index
     anaconda python script.py $rsfile_name_without_path $lon1 $lat1 > py_out.txt
     lat=`tail -1 py_out.txt| awk '{print $1}'`
     lat_index=`tail -1 py_out.txt| awk '{print $2+1}'`
     lat_integer=0
  fi


  if [[ "$lon1" =~ ^[0-9]+$ ]]; then # if lon1 is an integer
     lon_index=$lon1
     lon_integer=1
     lon=`ncdump -v lon $rsfile_name_without_path | sed -e '1,/data:/ d' | sed 's/lon = //g' | sed 's/}//g' | sed 's/;//g' | tr '\n' ' ' | awk -v L=$lon_index -F',' '{print $L}'`
  else
     lon=`head -1 py_out.txt| awk '{print $1}'`
     lon_index=`head -1 py_out.txt| awk '{print $2+1}'`
     lon_integer=0
  fi


  # check if one is integer and the other is not

  if [ $lat_integer -eq 1 -a $lon_integer -eq 0 -o $lat_integer -eq 0 -a $lon_integer -eq 1 ]; then
    echo " "
    echo "I'm confused ... the number types of lat ($lat1) and lon ($lon1) are not the same."
    echo "One is an integer and the other float"
    echo " "
    exit
  fi

  # remove leading and trailing  spaces from lat and lon
  # so that ncks works

  lat=`echo $lat | sed  's/^[ \t]*//;s/[ \t]*$//'`
  lon=`echo $lon | sed  's/^[ \t]*//;s/[ \t]*$//'`

  echo "--------------------------------------------------------------- "
  echo "Nearest lat and lon and their indices for your chosen grid cell"
  echo " "
  echo $lat [$lat_index], $lon [ $lon_index]
  echo " "

}


get_lat_lon_from_user


values_ever_changed_in_rsfile=0

see_and_modify_values_for_a_grid_cell () {

    if [ $first_time -ne 1 ]; then
       print_var_names
       choose_variable
    fi

    echo "--------------------------------------------------------------- "
    echo "Dimensions of your chosen variables and its values at your selected grid cell "
    echo " "

    #  Now extract values for the chosen variable and the chosen grid cell

    rm -f out1.nc
    ncks -v $varselected -d lon,${lon},${lon} -d lat,${lat},${lat} ${rsfile_name_without_path} out1.nc
    ncdump out1.nc | grep "${varselected}(tile"
    ncdump out1.nc | sed -e '1,/data:/ d' | grep -v "}"

    no_of_tiles=`ncdump out1.nc | sed -e '1,/data:/ d' | grep tile | awk -F'= ' '{print $2}' | awk -F' ;' '{print $1}'`
    # check if no_of_tiles empty or not
    if [ -z $no_of_tiles ]; then # if empty
      no_of_tiles=0
    fi
    # this script currently can't handle more than 1 tile
    if [ $no_of_tiles -gt 1 ]; then
       echo " "
       echo " "
       echo "It seems number of tiles is more than 1. This script can't go any further."
       echo " "
       exit
    fi



    ask_question_1 () {
      echo " "
      echo "--------------------------------------------------------------- "
      echo -e "\e[1;7mChange values for this grid cell ...\e[0m"
      echo " "
      read  -p "Would you like to change the values of [$varselected] for the grid cell you chose (y/n):" want_change
    }

    check_answer_1() {
       if [ "$want_change" != "y" -a "$want_change" != "n" ]; then
           echo " "
           echo "Only y or n are acceptable answers. "
           echo " "
           ask_question_1
           check_answer_1
       else
           useless=1
       fi
    }

    ask_question_1
    check_answer_1

    # if user wants variable to be not changed we exit, else we ask for values he/she wants
    how_to_continue() {

         if [ "$want_change" == "n" ]; then
           echo " "
           echo "---------------------------------------------------- "
           echo " "
         elif [ "$want_change" == "y" ]; then
           echo " "
           echo " "
           echo "--------------------------------------------------------------- "
           read  -p "Type your value(s) for [$varselected] (for mutiple values use comma to separate them) : " values_typed

           no_of_values_typed=`echo $values_typed | awk -F',' '{print NF; exit}'`

           echo " "
           echo " It seems you provided [$no_of_values_typed] value[s] that you want to replace [$varselected] with."
           echo " "
           read  -p "Are you happy with your values (y or n)? Type n to change them : " happy_with_vals_y_n

           if [ "$happy_with_vals_y_n" == "n" ]; then
             echo " "
             how_to_continue
           elif [ "$happy_with_vals_y_n" != "n" -a "$happy_with_vals_y_n" != "y" ]; then
             echo " "
             echo "Only y or n are acceptable answers. "
             echo " "
             sleep 0.5
             how_to_continue
           elif [ "$happy_with_vals_y_n" == "y" ]; then

             # find if the variable is 3 (tile,lat,lon) or 4 (tile, layer/pft, lat, lon) dimensional
             # or just 2 dimensional (lat,lon) like GC

             dimensions=`ncdump out1.nc | grep "${varselected}(tile" | awk -F'(' '{print $2}' | awk -F')' '{print $1}'`
             no_of_dims=`echo $dimensions| awk -F',' '{print NF; exit}'`

             #echo "no_of_dims = $no_of_dims"

             if [ $no_of_dims -eq 3 ]; then # SDEP (tile, lat, lon) like variable

               # for a 3D variable like SDEP  # of values type should be 1. I.E. there are no levels

               if [ $no_of_values_typed -ne 1 ]; then
                   echo " "
                   echo -e "\e[1;7mIt seems you provided $no_of_values_typed values and only 1 is needed ...\e[0m"
                   echo " "
                   continue
               fi

               change_values1 () {
                 echo " "
                 echo -e "\e[1;7mChange values for your grid cell or EVERYWHERE OVER LAND (based on GC < -0.5) ...\e[0m"
                 echo " "
                 read  -p "You can choose to change to your values for your [chosen (c)] grid cell or [everywhere over land (e)] on globe/domain (c or e) : " chosen_or_everywhere

                 if [ "$chosen_or_everywhere" == "c" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values for CHOSEN grid cell ...\e[0m"
                    ncap2 -s "$varselected(:,$lat_index-1,$lon_index-1)=$values_typed" $rsfile_name_without_path ${output_tmp}
                 elif [ "$chosen_or_everywhere" == "e" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values EVERYHERE (over land) in the file ...\e[0m"
                    #ncap2 -s "where(GC < -0.5) $varselected(:,:,:)=$values_typed" $rsfile_name_without_path ${output_tmp}
                    ncap2 -s "var_tmp= $varselected(:,:,:); where (GC < -0.5) var_tmp=${values_typed}; $varselected(:,:,:)=var_tmp;"  $rsfile_name_without_path ${output_tmp}
                 else
                   echo " "
                   echo -e "\e[1;7mIt seems you mistyped ...\e[0m"
                   echo -e "\e[1;7mYour options were c (for chosen grid cell) or e (for everywhere on the globe/domain) ...\e[0m"
                   echo " "
                   sleep 1
                   change_values1
                 fi
               }
               change_values1

             elif [ $no_of_dims -eq 4 ]; then # FCAN (tile, icp1, lat, lon) like variable

               # find no of levels | layers, pft, icp, icp1, icc, iccp1
               no_of_levels=`ncdump out1.nc | sed -e '1,/data:/ d' | grep "tile\|pft\|icp\|ic\|icp1\|icc\|iccp1\|layer\|months\|slope" | awk -F'= ' '{print $2}' | awk -F' ;' '{print $1}' | awk -F',' '{print NF; exit}'`

               #echo "no_of_levels= $no_of_levels"

               if [ $no_of_levels -gt $no_of_values_typed ]; then
                   echo " "
                   echo -e "\e[1;7mIt seems # of levels in $varselected in original file are $no_of_levels and you want to cut them down to $no_of_values_typed.\e[0m"
                   echo " "
                   echo "The script will now cut down the levels for this variable EVERYHWERE and not just at the chosen grid cell."
                   echo " "
                   read  -p "Do you want to continue cutting down the dimensions (y or n) : " cut_dimension_y_n

                   if [ "$cut_dimension_y_n" == "n" ]; then
                     echo " "
                     echo "Well then let's go back and start again!!"
                     echo "---------------------------------------------------- "
                     echo " "
                     sleep 1
                     continue
                   elif [ "$cut_dimension_y_n" == "y" ]; then
                     echo " "
                     k1=`expr $no_of_values_typed - 1`
                     ncks -d layer,0,${k1} ${rsfile_name_without_path} tmp1.nc
                     mv tmp1.nc ${rsfile_name_without_path}
                   fi

               elif [ $no_of_levels -lt $no_of_values_typed ]; then
                  echo " "
                  echo "--------------------------------------------------------------- "
                  echo "NCO doesn't allow to increase the levels for a given variable. "
                  echo "--------------------------------------------------------------- "
                  echo " "
                  sleep 1
                  continue
               fi

               # now change the values at the given grid cell

               change_values2 () {
                 echo " "
                 echo -e "\e[1;7mChange values for your grid cell or EVERYWHERE OVER LAND (based on GC < -0.5) ...\e[0m"
                 echo " "
                 read  -p "You can choose to change to your values for your [chosen (c)] grid cell or [everywhere over land (e)] on globe/domain (c or e) : " chosen_or_everywhere

                 if [ "$chosen_or_everywhere" == "c" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values for CHOSEN grid cell ...\e[0m"
                    for n in $(seq $no_of_values_typed); do
                         v1=`echo $values_typed | awk -v v=$n -F',' '{print $v}'`
                         echo "Changing value $n of $varselected to $v1"
                         ncap2 -s "$varselected(:,$n-1,$lat_index-1,$lon_index-1)=$v1" $rsfile_name_without_path tmp1.nc
                         mv tmp1.nc $rsfile_name_without_path
                    done
                 elif [ "$chosen_or_everywhere" == "e" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values EVERYHERE OVER LAND in the file ...\e[0m"
                    for n in $(seq $no_of_values_typed); do
                         v1=`echo $values_typed | awk -v v=$n -F',' '{print $v}'`
                         echo "Changing value $n of $varselected to $v1"
                         # ncap2 -s "where(GC < -0.5) $varselected(:,$n-1,:,:)=$v1" $rsfile_name_without_path tmp1.nc
                         n2=`expr $n - 1`
                         ncap2 -s "var_tmp=$varselected(:,${n2},:,:); where (GC < -0.5) var_tmp=${v1}; $varselected(:,${n2},:,:)=var_tmp;"  $rsfile_name_without_path tmp1.nc
                         mv tmp1.nc $rsfile_name_without_path
                    done
                 else
                   echo " "
                   echo -e "\e[1;7mIt seems you mistyped ...\e[0m"
                   echo -e "\e[1;7mYour options were c (for chosen grid cell) or e (for everywhere on the globe/domain) ...\e[0m"
                   echo " "
                   change_values2
                 fi
               }

               change_values2
               mv $rsfile_name_without_path ${output_tmp}


             elif [ $no_of_dims -eq 0 ]; then # GC (lat, lon) like variable

               # for a 2D variable like GC # of values type should be 1. I.E. there are no levels

               if [ $no_of_values_typed -ne 1 ]; then
                   echo " "
                   echo -e "\e[1;7mIt seems you provided $no_of_values_typed values and only 1 is needed ...\e[0m"
                   echo " "
                   echo " "
                   sleep 2
                   continue
                   #how_to_continue

               fi

               change_values3 () {
                 echo " "
                 echo -e "\e[1;7mChange values for your grid cell or EVERYWHERE OVER LAND (based on GC < -0.5) ...\e[0m"
                 echo " "
                 read  -p "You can choose to change to your values for your [chosen (c)] grid cell or [everywhere over land (e)] on globe/domain (c or e) : " chosen_or_everywhere

                 if [ "$chosen_or_everywhere" == "c" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values for CHOSEN grid cell ...\e[0m"
                    ncap2 -s "$varselected($lat_index-1,$lon_index-1)=$values_typed" $rsfile_name_without_path ${output_tmp}
                 elif [ "$chosen_or_everywhere" == "e" ]; then
                    echo " "
                    echo -e "\e[1;7mChanging values EVERYHERE in the file ...\e[0m"
                    ncap2 -s "where(GC < -0.5) $varselected=$values_typed" $rsfile_name_without_path ${output_tmp}
                 else
                   echo " "
                   echo -e "\e[1;7mIt seems you mistyped ...\e[0m"
                   echo -e "\e[1;7mYour options were c (for chosen grid cell) or e (for everywhere on the globe/domain) ...\e[0m"
                   echo " "
                   change_values3
                 fi
               }
               change_values3

             else
               echo " "
               echo " "
               echo -e "\e[1;7mSURPRISING but it seems this script doesn't know how to handle ${varselected}.\e[0m"
               echo " "
               echo " "
               exit
             fi # how many levels

           fi # happy_with_values
           mv ${output_tmp} $rsfile_name_without_path  # get ready to modify the file again if needed
           values_ever_changed_in_rsfile=1
         fi # want_change values

    } # how_to_continue


    how_to_continue

    first_time=0

} # see_and_modify_values_for_a_grid_cell


print_values_for_a_grid_cell () {

    sleep 0.5
    echo " "
    echo "Now printing values of all restart variables  to a text file ..."
    echo " "
    sleep 0.5


    # prefix of input rsfile
    prefix=`echo $rsfile_name_without_path | awk -F'.nc' '{print $1}'`

    output_txt="rsFile_vars_from_${prefix}_for_${lat_index}_${lon_index}.txt"
    output_nc="rsFile_vars_from_${prefix}_for_${lat_index}_${lon_index}.nc"
    rm -f ../$output_txt ../$output_nc
    touch -f ../$output_txt


    # throw all vars to a netcdf file as well for this grid cell
    rm -f out1.nc
    ncks -d lon,${lon},${lon} -d lat,${lat},${lat} ${rsfile_name_without_path} out1.nc
    nccopy -k classic out1.nc ../$output_nc


    count=1
    while read line; do
      v=`echo $line | awk -F'|' '{print $2}' | awk '{print $2}'`
      long_name=`echo $line | awk -F'|' '{print $3}'`
      echo "$count | [${v}] ------------------------------------------------------"

      #  Now extract values for the chosen variable and the chosen grid cell

      rm -f out1.nc
      ncks -v $v -d lon,${lon},${lon} -d lat,${lat},${lat} ${rsfile_name_without_path} out1.nc
      #ncdump out1.nc | grep "${v}(tile"
      #ncdump out1.nc | sed -e '1,/data:/ d' | grep -v "}"

      echo -e "\n----------- $v ($long_name} " >> ../${output_txt}
      ncdump out1.nc | grep "${v}(tile" >> ../${output_txt}
      ncdump out1.nc | sed -e '1,/data:/ d' | grep -v "}" >> ../${output_txt}

      count=`expr $count + 1`
    done < vars3.txt

    echo " "
    echo "Output txt file for grid cell [${lat_index}/${lon_index}] thrown to ............. ${output_txt}"
    echo "Output netcdf for grid cell [${lat_index}/${lon_index}] thrown to ............. ${output_nc}"
    echo " "
} # print_values_for_a_grid_cell


ask_if_they_want_to_do_it_again () {
   echo " "
   echo -e "\e[1;7mDo you want to see [and modify] values of another variable [for another grid cell]?\e[0m"
   echo " "
   read  -p "Choose (y) to continue and (n) to exit : " continue_again_1
}

ask_if_user_wants_grid_cell_to_change () {
   echo " "
   echo -e "\e[1;7mDo you want to choose a new grid cell ..\e[0m"
   echo " "
   read  -p "Choose (y) to change or press ENTER to continue with LAT $lat [$lat_index], LON $lon [$lon_index]: " change_grid_cell
}


do_it_again=1

if [ "x$whatwouldulike" = "x1" ]; then # user selected to see [and modify] values for a given grid cell

  while [ $do_it_again -eq 1 ]; do


       see_and_modify_values_for_a_grid_cell
       ask_if_they_want_to_do_it_again
       if [ "$continue_again_1" == "n" ]; then
         do_it_again=0
         echo ""
         echo "Good bye."
         echo ""
       elif [ "$continue_again_1" == "y" ]; then
         ask_if_user_wants_grid_cell_to_change
         if [ "$change_grid_cell" == "y" ]; then
           get_lat_lon_from_user
         elif [ "$change_grid_cell" == "n" -o "$change_grid_cell" == "" ]; then
           useless=1
         else
           echo " "
           echo "Only y or ENTER are acceptable answers. "
           echo " "
           ask_if_user_wants_grid_cell_to_change
         fi
       else
         echo " "
         echo "Only y or n are acceptable answers. "
         echo " "
         ask_if_they_want_to_do_it_again
       fi
  done

  # finally having made changes to the variables requested by user move $output_tmp to
  # $output_rs

  if [ $values_ever_changed_in_rsfile -eq 1 ]; then

      # also convert output rsFile to netcdf3 for easy viewing in Cview
      nccopy -k classic ${rsfile_name_without_path} ../${output_rs3}
      mv ${rsfile_name_without_path} ../${output_rs}
      echo " "
      echo "Output rsFile thrown to ............. ${output_rs} [NetCDF4]"
      echo "Output rsFile thrown to ............. ${output_rs3} [NetCDF3]"
      echo " "
  fi

elif [ "x$whatwouldulike" = "x2" ]; then # user selected to print values of all restart vars for  agiven grid cell
  first_time=1
  print_var_names
  print_values_for_a_grid_cell
fi
