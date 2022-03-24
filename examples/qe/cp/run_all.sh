#!/bin/sh

set -eu

examples=" \
    water
"

for dir in $examples; do
    cd "$dir"
    ./run_example.sh
    cd ../
done
