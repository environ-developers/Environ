#!/bin/bash
if [ -f tmp.txt ]; then
    rm tmp.txt
    fi
for i in *f90;
do
    modename=${i%%f90}o
    grep '^\s*USE' $i | grep -v mp | grep -v fft_base | \
	grep -v fft_types | grep -v scatter_mod | \
	grep -v fft_interfaces | \
	awk -v s=$modename '{printf("%s : %s.o\n",s,$2)}' | \
	sed s/','// | sed s/core_base/core_mod/ | \
	sed s/environ_base/environ_mod/ | \
	sed s/electrostatic_base/electrostatic_mod/ >> tmp.txt
done

cat -n tmp.txt | sort -uk2 | sort -nk1 | cut -f2- > make.depend
rm tmp.txt
