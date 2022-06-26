#!/bin/sh

set -eu

examples=" \
    lanczos \
    davidson \
"

for dir in $examples; do
    cd "$dir"
    ./run_example.sh
    cd ../
done