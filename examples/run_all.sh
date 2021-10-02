#!/bin/sh

for i in $(seq 10); do

    if [ "$i" -lt 10 ]; then
        cd example0"$i" || exit
    else
        cd example"$i" || exit
    fi

    ./run_example.sh

    cd ../
done
