#!/bin/bash

if [ $# != 6 ]; then
    echo "Parameters error"
    exit 1
fi

./pargen $1 $2 $3 $4 $5 $6
./netgen/netgen < netgen-$1-$2-$3-$4-$5-$6.par > netgen-$1-$2-$3-$4-$5-$6.dmx
./qfcgen netgen-$1-$2-$3-$4-$5-$6.dmx
