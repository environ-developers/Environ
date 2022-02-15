#!/bin/sh

for dir in example*/; do
    (cd "$dir" || exit; ./run_example.sh)
done
