#!/bin/sh

set -eu

for dir in example*; do
    cd "$dir"
    ./run_example.sh
    cd ../
done
