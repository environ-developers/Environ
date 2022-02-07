#!/bin/sh

set -eu

for dir in */; do
    cd "$dir"
    ./run_all.sh
    cd ../
done
