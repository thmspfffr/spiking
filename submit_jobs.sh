#!/bin/sh

echo -n "How many jobs do you want to submit? "
read NJOBS

for i in $( seq 1 $NJOBS); do
    let var1=10*$i;
    echo 'Start Job ' $i 'wait for: ' $var1 's'
    python2.7 run_model.py &
    sleep $var1
done

