# Makefile for Environ
# Adapted from PW main Makefile

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
