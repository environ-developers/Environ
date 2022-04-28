#!/bin/sh

set -eu

examples=" \
    SiO2 \
"

for dir in $examples; do
    cd "$dir"
    ./run_example.sh
    cd ../
done