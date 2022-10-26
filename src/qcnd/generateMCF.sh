#!/bin/bash

if [ $# != 6 ]; then
    echo "Parameters error"
    exit 1
fi

./pargen $1 $2 $3 $4 $5 $6 &&
echo "par generated!" &&
./netgen/netgen < netgen-$1-$2-$3-$4-$5-$6.par > netgen-$1-$2-$3-$4-$5-$6.dmx &&
echo "dmx generated!" &&
./qfcgen netgen-$1-$2-$3-$4-$5-$6.dmx &&
echo "qfc generated!" &&
python3 convertertomat.py netgen-$1-$2-$3-$4-$5-$6 &&
echo "mat generated!"


if [ -d "../data" ];
then
    mv netgen-$1-$2-$3-$4-$5-$6.* ../data
else
    mkdir ../data
    mv netgen-$1-$2-$3-$4-$5-$6.* ../data
fi