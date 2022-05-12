#!/bin/csh 

echo $argument > tmp
sed -e "s/=/ = /" tmp > fmp

awk '{print $1}' fmp > tmp
setenv argument `cat tmp`

awk '{print $3}' fmp > tmp
setenv variable `cat tmp`

rm -f tmp fmp
