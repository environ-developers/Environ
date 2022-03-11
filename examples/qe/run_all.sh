#!/bin/sh

set -eu

for dir in pw neb tddfpt cp; do
    cd "$dir"
    ./run_all.sh
    cd ../
done
