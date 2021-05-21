#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd $(echo $0 | sed 's/\(.*\)\/.*/\1/') # extract pathname
TOPDIR=$(pwd)

dirs=" UtilXlib FFTXlib src "

for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=$(echo ${dir%/})

    DEPENDS=

    # generate dependencies file (only for directories that are present)
    if test -d $TOPDIR/../$DIR; then
        cd $TOPDIR/../$DIR

        if test "$DIR" = "FFTXlib"; then
            DEPENDS="$DEPENDS ../UtilXlib"
        fi

        if test "$DIR" = "src"; then
            DEPENDS="$DEPENDS ../UtilXlib ../FFTXlib"
        fi

        $TOPDIR/moduledep.sh $DEPENDS >make.depend

        # handle special cases: modules for C-fortran binding,hdf5, MPI, FoX, libxc
        sed '/@iso_c_binding@/d' make.depend >tmp
        mv tmp make.depend
        sed '/@hdf5@/d' make.depend >tmp
        mv tmp make.depend
        sed '/@mpi@/d' make.depend >tmp
        mv tmp make.depend
        sed '/@fox_dom@/d;/@fox_wxml@/d;/@m_common_io@/d' make.depend >tmp
        mv tmp make.depend
        sed '/@xc_version.h@/d;/@xc_f03_lib_m@/d' make.depend >tmp
        mv tmp make.depend

        if test "$DIR" = "FFTXlib"; then
            # more special cases: modules for FFTs, GPU, OpenMP
            sed '/@omp_lib@/d' make.depend >tmp
            mv tmp make.depend
            sed '/@mkl_dfti/d' make.depend >tmp
            mv tmp make.depend
            sed '/@fftw3.f/d;s/@fftw.c@/fftw.c/' make.depend >tmp
            mv tmp make.depend
            sed '/@cudafor@/d;/@cufft@/d;/@flops_tracker@/d' make.depend >tmp
            mv tmp make.depend
        fi

        if test "$DIR" = "UtilXlib"; then
            sed '/@ifcore@/d' make.depend >tmp
            mv tmp make.depend
            sed '/@cudafor@/d' make.depend >tmp
            mv tmp make.depend
        fi

        # check for missing dependencies
        if grep @ make.depend; then
            notfound=1
            echo WARNING: dependencies not found in directory $DIR
        else
            echo directory $DIR : ok
        fi
    else
        echo directory $DIR : not present in $TOPDIR
    fi
done
if test "$notfound" = ""; then
    echo all dependencies updated successfully
fi
