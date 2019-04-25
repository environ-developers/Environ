#!/bin/bash
#
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.1
#
#    Environ 1.1 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.1 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# PATCH SCRIPT FOR Environ
#
# Author: Oliviero Andreussi (Department of Physics, UNT)
#

export ENVIRON_VERSION="1.1"
export ENVIRON_PATCH=$(cd "${BASH_SOURCE%/*}"; pwd)
export ENVIRON_DIR="${ENVIRON_PATCH}/.."
export ENVIRON_TEST="${ENVIRON_DIR}/tests"
export ENVIRON_SRC="${ENVIRON_DIR}/src"

# Installation script assumes that Environ directory is placed
# in main QE directory. Modify the following line if this is not the case

export QE_DIR="$ENVIRON_DIR/.."

# Set QE source directories and check that they are present

export PW_SRC="${QE_DIR}/PW/src"
if [ ! -d $PW_SRC ]; then
   echo "Cannot find PW/src directory"
   echo "Searching in $PW_SRC"
   exit
fi

export CP_SRC="${QE_DIR}/CPV/src"
if [ ! -d $CP_SRC ]; then
   echo "Cannot find CPV/src directory"
   echo "Searching in $CP_SRC"
   exit
fi

export TD_SRC="${QE_DIR}/TDDFPT/src"
if [ ! -d $TD_SRC ]; then
   echo "Cannot find TDDFPT/src directory"
   echo "Searching in $TD_SRC"
   exit
fi

export XS_SRC="${QE_DIR}/XSpectra/src"
if [ ! -d $XS_SRC ]; then
   echo "Cannot find XSpectra/src directory"
   echo "Searching in $XS_SRC"
   exit
fi

if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$0  (-patch or -revert) {all|pw|cp|td|xs}"
 echo " -patch  : apply Environ patch "
 echo " -revert : revert code to original "
 echo " a second argument may be used to specify the code to patch"
 exit
fi

case "$1" in
    (-patch)
	if [ "$#" -eq 1 ]; then
 	   for i in pw cp td xs ; do
	       PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
	       if test -e "$PATCH_SCRIPT" ; then
		   echo "-- Applying patches to $i"
		   bash "$PATCH_SCRIPT"
	       fi
	   done
	else
	    case "$2" in
	    (all)
 		for i in pw cp td xs ; do
		    PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
		    if test -e "$PATCH_SCRIPT" ; then
			echo "-- Applying patches to $i"
			bash "$PATCH_SCRIPT"
		    fi
		done
		;;
	    (pw)
		PATCH_SCRIPT="${ENVIRON_PATCH}/pw-patch.sh"
		if test -e "$PATCH_SCRIPT" ; then
		    echo "-- Applying patches to pw"
		    bash "$PATCH_SCRIPT"
		fi
		;;
	    (cp)
		PATCH_SCRIPT="${ENVIRON_PATCH}/cp-patch.sh"
		if test -e "$PATCH_SCRIPT" ; then
		    echo "-- Applying patches to cp"
		    bash "$PATCH_SCRIPT"
		fi
		;;
	    (td)
 		for i in pw td ; do
		    PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
		    if test -e "$PATCH_SCRIPT" ; then
			echo "-- Applying patches to $i"
			bash "$PATCH_SCRIPT"
		    fi
		done
		;;
	    (xs)
		for i in pw xs ; do
		    PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
		    if test -e "$PATCH_SCRIPT" ; then
			echo "-- Applying patches to $i"
			bash "$PATCH_SCRIPT"
		    fi
		done
	    esac
	fi
	;;
    (-revert)
	if [ "$#" -eq 1 ]; then
 	   for i in pw cp td xs ; do
	       PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
	       if test -e "$PATCH_SCRIPT" ; then
		   echo "-- Reverting patches to $i"
		   bash "$PATCH_SCRIPT"
	       fi
	   done
	else
	    case "$2" in
	    (all)
 		for i in pw cp td xs ; do
		    REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
		    if test -e "$REVERT_SCRIPT" ; then
			echo "-- Reverting patches to $i"
			bash "$REVERT_SCRIPT" > /dev/null
		    fi
		done
		;;
	    (pw)
		REVERT_SCRIPT="${ENVIRON_PATCH}/pw-revert.sh"
		if test -e "$REVERT_SCRIPT" ; then
		    echo "-- Reverting patches to pw"
		    bash "$REVERT_SCRIPT"
		fi
		;;
	    (cp)
		REVERT_SCRIPT="${ENVIRON_PATCH}/cp-revert.sh"
		if test -e "$REVERT_SCRIPT" ; then
		    echo "-- Reverting patches to cp"
		    bash "$REVERT_SCRIPT"
		fi
		;;
	    (td)
 		for i in pw td ; do
		    REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
		    if test -e "$REVERT_SCRIPT" ; then
			echo "-- Reverting patches to ${i}"
			bash "$REVERT_SCRIPT"
		    fi
		done
		;;
	    (xs)
		for i in pw xs ; do
		    REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
		    if test -e "$REVERT_SCRIPT" ; then
			echo "-- Reverting patches to ${i}"
			bash "$REVERT_SCRIPT"
		    fi
		done
	    esac
	fi
esac
