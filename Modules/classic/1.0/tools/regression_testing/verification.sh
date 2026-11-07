#!/bin/bash
testcase=$1
bdir=/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums
rf=${CI_BUILD_REF:0:8}

if [ ! -f $bdir/$testcase/checksums.csv ]; then
  echo "Error: No checksums.csv file found for $testcase. Aborting pipeline."
  exit 1
elif [ -z "$(cat $bdir/$testcase/checksums.csv)" ]; then
  echo "Error: checksums.csv file is empty for $testcase. Aborting pipeline."
  exit 1
fi

mkdir -p $bdir/checksums/$rf
mv $bdir/$testcase/checksums.csv $bdir/checksums/$rf/${testcase}_checksums.csv

casediff=$( diff $bdir/checksums/$rf/${testcase}_checksums.csv $bdir/checksums/baseline/${testcase}_checksums.csv )

if [ "$casediff" != "" ]; then
  diff $bdir/checksums/$rf/${testcase}_checksums.csv $bdir/checksums/baseline/${testcase}_checksums.csv
  echo "$testcase checksum verification has failed"
  exit 1
fi
echo "$testcase checksum verification passed"
exit 0
