#!/usr/bin/env python3
import sys
import os
import shutil

# USE: regtest.py $RUNPATH test_name

def main():

  # get the path where the netcdf folders are located
  output_directory = sys.argv[1] + '/classic_checksums/' + sys.argv[2]
  internal_dirs = os.listdir(output_directory)
  files = []
  netdirs = []

  # put all.csv files in all the netcdf folders into an array
  for odir in [x for x in internal_dirs if 'netcdf' in x]:
    netdir = output_directory + '/{}'.format(odir) + '/checksums'
    if os.path.isdir(netdir):
      for file in os.listdir(netdir):
        files.append(netdir + '/' + file)
      netdirs.append(netdir)

  # read the lines from all files into a list
  lines = []
  with open(output_directory + '/checksums.csv', 'w') as f:
    for file in files:
      with open(file) as g:
        for line in g.readlines():
          lines.append(line)
    # lines must be sorted by the first, second, and last element
    # (for easy comparison .csv file comparison)
    splitLines = [line.split(',') for line in lines]
    parsedLines = []
    for line in splitLines:
      newline = []
      for i, item in enumerate(line):
        if i < 2:
          newline.append(int(item))
        else:
          newline.append(item)
      parsedLines.append(newline)
    parsedLines = sorted(parsedLines, key=lambda x : x[2])
    parsedLines = sorted(parsedLines, key=lambda x : x[1])
    parsedLines = sorted(parsedLines, key=lambda x : x[0])

    # turn the integer values in the line back into string values
    completeLines = []
    for line in parsedLines:
      newline = []
      for i, item in enumerate(line):
        if i < 2:
          newline.append(str(item))
        else:
          newline.append(item)
      completeLines.append(",".join(newline))

    # now write the parsed, sorted lines into the new .csv file
    for line in completeLines:
      f.write(line)


if __name__ == '__main__':
    main()
