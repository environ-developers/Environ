#!/bin/bash

# TODO: add check for long tests

runtime() {
    grep "PWSCF *:" "$1" | grep -Eo "[0-9]{2}\.[0-9]{2}" | tail -1
}

# LOAD TEST DIRECTORIES
if [[ "$*" ]]; then
    testdirs="$*"
else
    mapfile -t testdirs <<<"$(find . -name "pw_*")"
fi

# START PROCESS
for dir in "${testdirs[@]}"; do
    printf "\n%10s\n" "$dir" # HEADER

    # ------------------------------- INPUT VALIDATION
    if [ -z "$dir" ]; then
        echo "Missing directory"
        continue
    fi

    if ! compgen -G "$dir/benchmark*" >/dev/null; then
        echo "Missing benchmarks"
        continue
    fi

    if ! compgen -G "$dir/test\.out*" >/dev/null; then
        echo "Missing tests"
        continue
    fi
    # ------------------------------------------------

    mapfile -t files <<<"$(find "$dir" -name "benchmark*")"
    printf "%50s\t%11s\n" "Test" "delta_t (s)"
    for file in "${files[@]}"; do
        inp=$(echo "$file" | cut -d "=" -f 2-) # JOB NAME
        bench_t=$(runtime "$file")
        test_t=$(runtime "$(find "$dir" -name "test.out*$inp*")")
        delta="$(echo "$bench_t" "$test_t" | awk '{printf "%.2f", $2 - $1}')"
        printf "%50s\t%7s\n" "$inp" "$delta"
    done
done
