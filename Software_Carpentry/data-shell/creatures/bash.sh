#!/bin/bash
for filename in "$2";
do
  head -n 5 $filename
done
