#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2017 Quantum ESPRESSO (www.quantum-espresso.org)
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 3.0
#
#     Environ 3.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#
#     Environ 3.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Oliviero Andreussi (Department of Physics, UNT)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# This file was adapted from the following Quantum ESPRESSO v6.1 file:
#
#     QE/test-suite/check_pseudo.sh
#
#----------------------------------------------------------------------------------------

if test "`which curl`" = "" ; then
   if test "`which wget`" = "" ; then
      echo "### wget or curl not found: will not be able to download missing PP ###"
   else
      DOWNLOADER="wget -O"
      # echo "wget found"
   fi
else
   DOWNLOADER="curl -o"
   # echo "curl found"
fi

inputs=`find $1* -type f -name "*.in" -not -name "test.*" -not -name "benchmark.*"`
pp_files=`for x in ${inputs}; do grep UPF ${x} | awk '{print $3}'; done`

for pp_file in ${pp_files} ; do
if ! test -f ${ESPRESSO_PSEUDO}/${pp_file} ; then
	#echo -n "Downloading ${pp_file} to ${ESPRESSO_PSEUDO} ... "
	${DOWNLOADER} ${ESPRESSO_PSEUDO}/${pp_file} ${NETWORK_PSEUDO}/${pp_file} 2> /dev/null
	if test $? != 0 ; then
		echo "FAILED, do it manually -- Testing aborted!"
		exit -1
	#else
		#echo "done."
	fi
#else
	#echo "No need to download ${pp_file}."
fi
done
