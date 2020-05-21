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
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#

default: all

all: doc libenviron

doc:
	if test -d Doc ; then \
        (cd Doc ; $(MAKE) || exit 1 ) ; fi

libenviron:
	if ! test -d ../PW ; then \
		( make libfft ) ; \
		( make libutil ) ; fi
	if test -d src ; then \
        ( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
        else $(MAKE) $(MFLAGS) ; fi ) ; fi ; \

libfft : 
	( cd libs/FFTXlib ; $(MAKE) TLDEPS= all || exit 1 )

libutil : 
	( cd libs/UtilXlib ; $(MAKE) TLDEPS= all || exit 1 )

clean :
	if test -d src ; then \
        ( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
        else $(MAKE) clean ; fi ) ; fi ;\
	( cd libs/FFTXlib ; $(MAKE) clean ) ; \
	( cd libs/UtilXlib ; $(MAKE) clean )

doc_clean:
	if test -d Doc ; then \
        (cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
