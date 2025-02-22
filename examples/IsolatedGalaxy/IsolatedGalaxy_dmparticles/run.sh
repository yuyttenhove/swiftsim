#!/bin/bash

if [ ! -e 3e11-star-only-DM-halo-galaxy.hdf5 ]
then
    echo "Fetching initial conditons for the isolated galaxy with an external potential ..."
    ./getIC.sh
fi 

../../swift --external-gravity --self-gravity --stars --threads=4 isolated_galaxy.yml 2>&1 | tee output.log


echo "Make plots of conservation of total angular momentum" 
if command -v python3 &>/dev/null; then
    python3 angularmomentum.py 
else
    python angularmomentum.py 
fi

echo "Make plots of change of vertical and radial profile"
if command -v python3 &>/dev/null; then
    python3 profilefit.py 
else
    python profilefit.py 
fi
