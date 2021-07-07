#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 2.0
#     
#     Environ 2.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#     
#     Environ 2.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Edan Bainglass (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Runtime comparison for tests
#
#----------------------------------------------------------------------------------------

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
