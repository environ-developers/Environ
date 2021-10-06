#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2017 Quantum ESPRESSO (www.quantum-espresso.org)
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
# Authors: Oliviero Andreussi (Department of Physics, UNT)
#		   Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Compute dependencies for the PWscf directory tree
#
#----------------------------------------------------------------------------------------

# remove unwanted dependencies
no_deps() {
    for no_dep in $*; do
        echo "/@$no_dep@/d" >>removedeps.tmp
    done
    sed -f removedeps.tmp make.depend >tmp
    mv tmp make.depend
    rm removedeps.tmp
}

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd "$(echo "$0" | sed 's/\(.*\)\/.*/\1/')" || exit # extract pathname
TOPDIR=$(pwd)

dirs="utils FFTs src"

for dir in $dirs; do

    DEPENDS=

    # generate dependencies file (only for directories that are present)
    if test -d "$TOPDIR"/../"$dir"; then
        cd "$TOPDIR"/../"$dir" || exit

        case "$dir" in
        FFTs) DEPENDS="$DEPENDS ../utils" ;;
        src) DEPENDS="$DEPENDS ../utils ../FFTs" ;;
        esac

        "$TOPDIR"/moduledep.sh "$DEPENDS" >make.depend

        # list of all system modules
        sysdeps="iso_c_binding ifcore"

        # list of all external library modules or include files
        libdeps="mpi omp_lib mkl_dfti mkl_dfti.f90 fftw3.f03 fftw3.f"

        # list of all cuda-related modules
        cudadeps="cudafor cufft flops_tracker"

        no_deps "$sysdeps $libdeps $cudadeps"

        # check for missing dependencies
        if grep @ make.depend; then
            notfound=1
            echo WARNING: dependencies not found in directory "$dir"
        else
            echo directory "$dir" : ok
        fi
    else
        echo directory "$dir" : not present in "$TOPDIR"
    fi
done

if test "$notfound" = ""; then
    echo all dependencies updated successfully
fi
