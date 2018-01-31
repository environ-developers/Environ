# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.0
#
#    Environ 1.0 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# PATCH SCRIPT FOR Environ
#
# This script has been adapted from original patch scripts
# of plumed (www.plumed-code.org) by Layla Martin-Samos
#
#!/bin/bash

ENVIRON_VERSION="1.0"
ENVIRON_dir="$PWD/../../Environ"
patch_SCRIPT="$ENVIRON_dir/patches/qe-patch.sh"
revert_SCRIPT="$ENVIRON_dir/patches/qe-revert.sh"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  echo > /dev/null
}

function to_do_before_revert () {
  echo > /dev/null
}

function to_do_after_revert () {
  echo > /dev/null
}


if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$0  (-patch) (-revert)   "
 echo " -patch  : apply Environ patch "
 echo " -revert : revert code to original "
 exit
fi

case "$1" in
(-patch)
  echo "* I will try to patch Environ version $ENVIRON_VERSION ..."
  if test -e "Environ_PATCH" ; then
    echo "-- File Environ_PATCH exists"
    echo "-- I guess you have already patched Environ $(tail -1 Environ_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
  fi
  echo "#Please do not remove or modify this file"                    >  Environ_PATCH
  echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
  echo "$ENVIRON_VERSION"                                              >> Environ_PATCH

  to_do_before_patch

  if test -e "$patch_SCRIPT" ; then
    echo "-- Applying patches"
    bash "$patch_SCRIPT"
  fi

  to_do_after_patch

  echo "- DONE!"
;;

(-revert)
  echo "* I will try to revert Environ version $Environ_VERSION ..."
  if test ! -e Environ_PATCH ; then
    echo "-- File Environ_PATCH is not there"
    echo "-- I guess you never patched, so there is nothing to revert"
    echo "* ABORT"
    exit
  fi

  echo "-- Executing pre script"

  to_do_before_revert

  if test -e "$revert_SCRIPT" ; then
    echo "-- Reverting patches"
    bash "$revert_SCRIPT"
  fi

  echo "-- Executing post script"
  to_do_after_revert

  rm "Environ_PATCH"

  echo "* DONE!"

esac
