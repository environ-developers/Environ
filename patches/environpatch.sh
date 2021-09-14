#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
# PATCH script for Environ
#
#----------------------------------------------------------------------------------------

export ENVIRON_VERSION="2.0"

ENVIRON_PATCH=$(
	cd "${BASH_SOURCE%/*}" || exit
	pwd
)
export ENVIRON_PATCH

export ENVIRON_DIR="${ENVIRON_PATCH}/.."
export ENVIRON_TEST="${ENVIRON_DIR}/tests"
export ENVIRON_SRC="${ENVIRON_DIR}/src"

# other variables used in patch/revert scripts
export PATCHED=0
export REVERTED=0

# local variables
DOTS='........................................'
PW_DEP_DIRS="NEB/src \
			 PP/src \
			 PHonon/PH \
			 PHonon/Gamma \
			 PWCOND/src HP/src \
			 GWW/pw4gww \
			 GWW/head \
			 GWW/bse \
			 GWW/simple"

# Installation script assumes that Environ directory is placed
# in main QE directory. Modify the following line if this is not the case

export QE_DIR="$ENVIRON_DIR/.."

# Set QE source directories and check that they are present

export PW_SRC="${QE_DIR}/PW/src"
if [ ! -d "$PW_SRC" ]; then
	echo "Cannot find PW/src directory"
	echo "Searching in $PW_SRC"
	exit
fi

export CP_SRC="${QE_DIR}/CPV/src"
if [ ! -d "$CP_SRC" ]; then
	echo "Cannot find CPV/src directory"
	echo "Searching in $CP_SRC"
	exit
fi

export TD_SRC="${QE_DIR}/TDDFPT/src"
if [ ! -d "$TD_SRC" ]; then
	echo "Cannot find TDDFPT/src directory"
	echo "Searching in $TD_SRC"
	exit
fi

export XS_SRC="${QE_DIR}/XSpectra/src"
if [ ! -d "$XS_SRC" ]; then
	echo "Cannot find XSpectra/src directory"
	echo "Searching in $XS_SRC"
	exit
fi

if [ "$#" -eq 0 ]; then
	echo "USAGE :"
	echo "$0  (-patch or -revert) {all|pw|cp|td|xs}"
	echo " -patch  : apply Environ patch "
	echo " -revert : revert code to original "
	echo " a second argument may be used to specify the code to patch"
	exit
fi

# functions used in patch/revert scripts

fill_with_dots() {
	printf "%s" "$1${DOTS:${#1}}"
}

check_src_patched() {
	if test -e "Environ_PATCH"; then
		echo "  - src already patched!"
		PATCHED=1
	else
		PATCHED=0
	fi
}

patch_makefile() {
	if test "$(grep '# Environ patch' Makefile)"; then
		echo "  - $1Makefile already patched!"
		return
	else
		fill_with_dots "  - Patching $1Makefile"
		mod1="MODFLAGS+=\$(MOD_FLAG)../../Environ/src \$(MOD_FLAG)../../Environ/utils"
		mod2="QEMODS+=../../Environ/libs/*"

		sed -i.tmp -e '/^TLDEPS/a \
# Environ patch \
'"$mod1"' \
'"$mod2"' \
# end Environ patch' Makefile && rm Makefile.tmp

		printf " done! \n"
	fi
}

message() {
	fill_with_dots "  - $1 src"
}

check_src_reverted() {
	if test ! -e Environ_PATCH; then
		echo "  - src has not been patched!"
		REVERTED=1
	else
		REVERTED=0
	fi
}

revert_makefile() {
	if test "$(grep '# Environ patch' Makefile)"; then
		fill_with_dots "  - Reverting $1Makefile"
		sed -i.tmp -e '/# Environ patch/,/# end Environ patch/d' Makefile && rm Makefile.tmp
		printf " done! \n"
	else
		echo "  - $1Makefile has not been patched!"
		return
	fi
}

patch_makedeps() {

	case $1 in
	pw) DEPS="PW/src" ;;
	cp) DEPS="CPV/src" ;;
	tddfpt) DEPS="PW/src | TDDFPT/src" ;;
	xspectra) DEPS="PW/src | XSpectra/src" ;;
	*) DEPS="PW/src | TDDFPT/src | XSpectra/src | CPV/src" ;;
	esac

	sed -i.tmp -e '/cd $TOPDIR\/..\/$DIR/a \
	# Environ patch\
	case $DIR in\
	'"$DEPS"')\
		DEPENDS="$DEPENDS $LEVEL2/Environ/src $LEVEL2/Environ/utils"\
		;;\
	esac\
	# end Environ patch' $file && rm $file.tmp
}

case "$1" in
-patch)

	# patch to QE/install/makedeps.sh
	file="../install/makedeps.sh"

	if test "$(grep '# Environ patch' $file)"; then
		sed -i.tmp -e '/# Environ patch/,/# end Environ patch/d' $file && rm $file.tmp
	fi

	printf "\n* Patching install/makedeps.sh.........."
	patch_makedeps "$2"
	printf " done! \n\n"

	case "$2" in
	pw)
		PATCH_SCRIPT="${ENVIRON_PATCH}/pw-patch.sh"
		if test -e "$PATCH_SCRIPT"; then
			echo "* Applying patches to pw"
			. "$PATCH_SCRIPT"
		fi
		;;
	cp)
		PATCH_SCRIPT="${ENVIRON_PATCH}/cp-patch.sh"
		if test -e "$PATCH_SCRIPT"; then
			echo "* Applying patches to cp"
			. "$PATCH_SCRIPT"
		fi
		;;
	tddfpt)
		for i in pw td; do
			PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
			if test -e "$PATCH_SCRIPT"; then
				echo "* Applying patches to $i"
				. "$PATCH_SCRIPT"
			fi
		done
		;;
	xspectra)
		for i in pw xs; do
			PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
			if test -e "$PATCH_SCRIPT"; then
				echo "* Applying patches to $i"
				. "$PATCH_SCRIPT"
			fi
		done
		;;
	*)
		for i in pw cp td xs; do
			PATCH_SCRIPT="${ENVIRON_PATCH}/${i}-patch.sh"
			if test -e "$PATCH_SCRIPT"; then
				echo "* Applying patches to $i"
				. "$PATCH_SCRIPT"
			fi
		done
		;;
	esac

	cd "$ENVIRON_DIR" || exit # return to Environ root directory

	# apply Makefile patches to QE packages that rely on pw.x
	printf "\n* Patching Makefiles dependent on pw.x\n"
	for i in $PW_DEP_DIRS; do
		loc=../$i
		if test -e "$loc"; then
			(
				cd "$loc" || exit
				patch_makefile "$i/"
			)
		fi
	done

	;;
-revert)

	# revert patch to QE/install/makedeps.sh
	file="../install/makedeps.sh"
	if test "$(grep '# Environ patch' $file)"; then
		printf "\n* Reverting install/makedeps.sh........."
		sed -i.tmp -e '/# Environ patch/,/# end Environ patch/d' $file && rm $file.tmp
		printf " done! \n\n"
	else
		printf "\n* install/makedeps.sh has not been patched! \n\n"
	fi

	case "$2" in
	pw)
		REVERT_SCRIPT="${ENVIRON_PATCH}/pw-revert.sh"
		if test -e "$REVERT_SCRIPT"; then
			echo "* Reverting patches to pw"
			. "$REVERT_SCRIPT"
		fi
		;;
	cp)
		REVERT_SCRIPT="${ENVIRON_PATCH}/cp-revert.sh"
		if test -e "$REVERT_SCRIPT"; then
			echo "* Reverting patches to cp"
			. "$REVERT_SCRIPT"
		fi
		;;
	tddfpt)
		for i in pw td; do
			REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
			if test -e "$REVERT_SCRIPT"; then
				echo "* Reverting patches to ${i}"
				. "$REVERT_SCRIPT"
			fi
		done
		;;
	xspectra)
		for i in pw xs; do
			REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
			if test -e "$REVERT_SCRIPT"; then
				echo "* Reverting patches to ${i}"
				. "$REVERT_SCRIPT"
			fi
		done
		;;
	*)
		for i in pw cp td xs; do
			REVERT_SCRIPT="${ENVIRON_PATCH}/${i}-revert.sh"
			if test -e "$REVERT_SCRIPT"; then
				echo "* Reverting patches to $i"
				. "$REVERT_SCRIPT"
			fi
		done
		;;
	esac

	cd "$ENVIRON_DIR" || exit # return to Environ root directory

	# revert Makefile patches to QE packages that rely on pw.x
	printf "\n* Reverting Makefiles dependent on pw.x\n"
	for i in $PW_DEP_DIRS; do
		loc=../$i
		if test -e "$loc"; then
			(
				cd "$loc" || exit 1
				revert_makefile "$i/"
			)
		fi
	done

	;;
esac
