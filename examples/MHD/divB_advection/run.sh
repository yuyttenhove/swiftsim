#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e divB.hdf5 ]
then
    echo "Generating initial conditions for the divB advection example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift  --drift-all --hydro --threads=4 divB.yml 2>&1 | tee output.log

python3 plotSolution.py 11
