#!/usr/bin/env sh

###################################################################################################

# check if library is already modified
already_modified() {
    test "$(grep env_ "$lib"/*.f90)" && return 0 || return 1
}

# apply 'env_' prefix to object
modify() {
    find "$lib" -maxdepth 1 \( -name "*.f90" -o -name "*.c" \) -print0 | xargs -0 sed -i "s/\<$1\>/env_$1/gi"

    if test "$lib" = UtilXlib && test ! "$1" = errore; then
        find FFTXlib -maxdepth 1 \( -name "*.f90" -o -name "*.c" \) -print0 | xargs -0 sed -i "s/\<$1\>/env_$1/gi"
    fi

    # clean bad mods
    if test "$lib" = FFTXlib && test "$1" = allowed; then
        sed -i 's/not env_allowed/not allowed/g' "$lib"/*.f90
        sed -i 's/Max env_allowed fft dimension/Max allowed fft dimension/' "$lib"/*.f90
    fi
}

# revert modifications
revert_modification() {

    # undo info broadcasting
    if test "$lib" = FFTXlib; then
        sed -i \
            -e '/USE env_mp, ONLY: env_mp_bcast/d' \
            -e '/CALL env_mp_bcast(info, dfft%root, dfft%comm)/d' "$lib"/scatter_mod.f90
    fi

    # remove env_ prefix
    sed -i 's/env_//g' "$lib"/*.f90 2>/dev/null

}

# remove extra whitespace including tabs
strip() {
    sed -e "s/\t/ /g" -e 's/  */ /g' -e 's/^ //'
}

# return unique list of object names
clean_results() {
    rev | strip | cut -d ' ' -f 1 | rev | uniq
}

# get files to modify
get_files() {
    find "$lib" -maxdepth 1 \( -name "*.f90" -o -name "*.c" \) -print0
}

# get objects to modify
get_objects() {
    get_files | xargs -0 grep -hiE "$1"
}

# remove current library files
clean_library() {
    rm "$lib"/*.f90 "$lib"/*.c 2>/dev/null
}

# remove references to discarded modules
discard_modules() {
    if test "$lib" = FFTXlib; then
        sed -i '/SUBROUTINE invfft_b/,/END SUBROUTINE invfft_b/d' FFTXlib/*.f90 2>/dev/null
        sed -i '/elif defined(__FFTW)/,/USE fft_scalar_fftw/d' FFTXlib/*.f90 2>/dev/null
        sed -i '/fft_smallbox/d' FFTXlib/*.f90 2>/dev/null

    elif test "$lib" = UtilXlib; then
        sed -i '/error_handler/,/error_handler/d' UtilXlib/error_handler.f90 2>/dev/null
    fi

}

# prefix env_ to modules, routines, and interfaces
apply_prefix() {
    modules="$(get_objects '^\s*module ' | awk 'BEGIN {IGNORECASE = 1} {if ($2!="procedure") print $0}' | clean_results)"
    derived_types="$(get_objects '^\s*type(\s*,\s*\w+)*\s*(:{2})?\s*[^(]\w+' | clean_results)"
    subroutines="$(get_objects '^\s*subroutine ' | cut -d '(' -f 1 | clean_results)"
    functions="$(get_objects '^\s*(pure|recursive|integer|logical|(complex|real)\s*\(\w+\)\s*)?\s*function ' | grep -ioE 'function \w+' | clean_results)"
    interfaces="$(get_objects '^\s*(abstract)?\s*interface\s+\w+' | clean_results)"

    if test "$lib" = UtilXlib; then
        c_routines="$(get_objects 'cclock|scnds' | cut -d '(' -f 1 | clean_results)"
    fi

    for object in \
        $modules \
        $derived_types \
        $subroutines \
        $functions \
        $interfaces \
        $c_routines; do
        modify "$object"
    done
}

# fix needed by PyE-C
broadcast_info() {
    if test "$lib" = FFTXlib; then
        sed -i \
            -e '/USE env_fft_param/a \        USE env_mp, ONLY: env_mp_bcast' \
            -e '/info < 0/i \  CALL env_mp_bcast(info, dfft%root, dfft%comm)' FFTXlib/scatter_mod.f90
    fi
}

# fix mp buffer deallocation
add_allocated_checks() {
    if test "$lib" = UtilXlib; then
        sed -i 's/\(DEALLOCATE(mp_buff_r, mp_buff_i)\)/IF (ALLOCATED (mp_buff_r) .AND. ALLOCATED (mp_buff_i)) \1/' UtilXlib/mp_base.f90
        sed -i 's/\(DEALLOCATE(mp_buff_r_d, mp_buff_i_d)\)/IF (ALLOCATED (mp_buff_r_d) .AND. ALLOCATED (mp_buff_i_d)) \1/' UtilXlib/mp_base_gpu.f90
    fi
}

# copy library files over from QE
copy_files() {

    if test "$lib" = FFTXlib; then

        find "$1/$lib" -maxdepth 1 -name '*.f90' \
            ! -name '*smallbox*' \
            ! -name 'fft_test.f90' \
            ! -name 'fftw_interfaces.f90' \
            ! -name 'fft_scalar.FFTW.f90' \
            -exec cp -t "$lib" {} +

    elif test "$lib" = UtilXlib; then

        # fortran files
        find "$1/$lib" -maxdepth 1 -name '*.f90' \
            ! -name 'clib_wrappers.f90' \
            ! -name 'device_helper.f90' \
            ! -name 'divide.f90' \
            ! -name 'export_gstart_2_solvers.f90' \
            ! -name 'fletcher32_mod.f90' \
            ! -name 'hash.f90' \
            ! -name '*mem*' \
            ! -name 'mp_bands_util.f90' \
            ! -name 'set_mpi_comm_4_solvers.f90' \
            ! -name 'thread_util.f90' \
            -exec cp -t "$lib" {} +

        # c files
        cp "$1/$lib"/cptimer.c "$lib"

    fi
}

# check files exist to be modified
check_files() {
    if test "$(find "$lib" -maxdepth 1 \( -name "*.f90" -o -name "*.c" \) 2>/dev/null)"; then
        return 0
    else
        printf "\nNo files found\n\n" && return 1
    fi
}

check_path() {
    if test -d "$1/$lib"; then
        return 0
    else
        printf "\n$lib not found at %s\n\n" "$1" && return 1
    fi
}

# print script usage message
help_message() {
    printf "\nupdate -h for help\n\n"
}

###################################################################################################

test ! "$1" && help_message && exit

while getopts "hcmrl:p:v:" opt; do

    case "$opt" in
    h)
        echo
        echo "usage: update_libs.sh <flag> <option>"
        echo
        echo "To update library, run"
        echo
        echo "   update_libs.sh -l <FFTXlib|UtilXlib> -p <qeroot> -m -v ###"
        echo
        echo "-c             clean current library state (rm *.f90)"
        echo
        echo "-l <library>   library to update"
        echo
        echo "-p <qe-path>   copy files from <qe-path>/FFTXlib to Environ/FFTXlib"
        echo
        echo "-m             apply 'env_' prefix modification"
        echo
        echo "-r             revert 'env_' prefix modification"
        echo
        echo "-v ###         update version in README.md to ###"
        echo
        ;;
    c)
        printf "\nClean library (*.f90) (y|n)? "
        read -r o
        test "$o" = y && clean_library && echo && exit
        ;;
    m)
        check_files || exit
        already_modified && exit
        discard_modules
        apply_prefix
        broadcast_info
        add_allocated_checks
        ;;
    r)
        revert_modification && exit
        ;;
    l)
        lib="$OPTARG"
        if test "$lib" = fft; then
            lib=FFTXlib
        elif test "$lib" = util; then
            lib=UtilXlib
        else
            printf "\nBad library argument -> %s\n\n" "$lib"
            exit
        fi
        ;;
    p)
        test $lib || exit
        check_path "$OPTARG" || exit
        clean_library
        copy_files "$OPTARG"
        ;;
    v)
        # update library version number
        version=$(grep -oE '[0-9][0-9.]+' "$lib"/README.md)
        sed -i "s/$version/$OPTARG/" "$lib"/README.md
        ;;
    *)
        help_message
        ;;
    esac

done
