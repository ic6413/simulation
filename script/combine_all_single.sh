#!/bin/bash

cd output/single_all

START=0
END=10000
STEPDIFF=1

# delete headers
for i in $(eval echo "{$START..$END..$STEPDIFF}")
do
sed "-e 1,9d; s/^/$i /" dump.all.single.$i > tmp$i
done

cat $(eval echo "tmp{$START..$END..$STEPDIFF}") > dump.all.single.allstep

rm -r tmp*