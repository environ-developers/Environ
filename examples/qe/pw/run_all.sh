#!/bin/sh

set -eu

examples=" \
    sccs \
    sscs \
    pbc \
    slab \
    helmholtz \
    helmholtz_linpb \
    mott_schottky \
    solvent_aware \
    field_aware \
"

for dir in $examples; do
    cd "$dir"
    ./run_example.sh
    cd ../
done
