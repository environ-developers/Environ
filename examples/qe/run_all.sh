#!/bin/sh

set -eu

for dir in pw td cp; do
    cd "$dir"
    ./run_all.sh
    cd ../
done
