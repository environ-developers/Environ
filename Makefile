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
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#

default: all

all: libenviron doc

doc:
	if test -d Doc ; then \
        (cd Doc ; $(MAKE) || exit 1 ) ; fi

libenviron:
	if test -d src ; then \
        ( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
        else $(MAKE) $(MFLAGS) ; fi ) ; fi ; \

clean :
	if test -d src ; then \
        ( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
        else $(MAKE) clean ; fi ) ; fi ;\

doc_clean:
	if test -d Doc ; then \
        (cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
