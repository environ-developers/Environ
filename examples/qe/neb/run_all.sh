#!/bin/sh

set -eu

examples=" \
    collinear_proton_transfer
"

for dir in $examples; do
    cd "$dir"
    ./run_example.sh
    cd ../
done
