#!/bin/bash
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
# Copyright (C) 2001-2017 Quantum ESPRESSO group
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
# This file was adapted from the following Quantum ESPRESSO v6.1 file
#
#      QE/test-suite/check_pseudo.sh
#
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#

source ${ENVIRON_ROOT}/tests/ENVIRONMENT

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
