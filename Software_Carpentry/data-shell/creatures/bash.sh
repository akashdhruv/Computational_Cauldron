#!/bin/bash
for filename in *.dat;
do
  head -n 5 $filename
done
