#!/bin/bash
JOB_ID=$1
echo "Job ID: $JOB_ID"
sleep 20
let sleepcount=0
while [ $sleepcount -lt 1080 ] ;
do
  status=$( jobst -c $HDNODE 2>/dev/null | grep $1 | tr -s '|' ' ' | cut -d ' ' -f 3 )
  status=${status#*|}
  if [ "$status" == "Q" ]; then
    echo "Queued!"
  elif [ "$status" == "R" ]; then
    echo "Running!"
  elif [ -z "$status" ]; then
    break
  else
    echo "Unexpected jobstat result:"
    echo $( jobst -c $HDNODE 2>/dev/null | grep $1 )
  fi
  let sleepcount=$sleepcount+1
  sleep 5
done

[ $sleepcount -eq 1080 ] && echo "Error: Sleep count timed out. The job must be stuck in the queue. Checksums will not be checked." && exit 1
exit 0
