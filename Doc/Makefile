include ../make.inc

HELPDOC=$(QEDIR)/dev-tools/helpdoc

doc:  all
all:  defs

clean:
	@ rm -f input_xx.xsl
	@ rm -f INPUT_*.html INPUT_*.xml
	@ rm -f INPUT_Environ.txt
	@ rm -f $(QEDIR)/Doc/INPUT_Environ.*


defs: link_input_xx INPUT_Environ.html INPUT_Environ.txt link_on_main_doc

link_input_xx:
	@(if test ! -f input_xx.xsl; then \
	(if test -f $(QEDIR)/dev-tools/input_xx.xsl; then \
	(ln -sf $(QEDIR)/dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; \
	fi) ; fi)

INPUT_Environ.html: %.html: %.def
	$(HELPDOC) $<
INPUT_Environ.txt: %.txt: %.def
	$(HELPDOC) $<

link_on_main_doc:
	-@( cd $(QEDIR)/Doc ; ln -fs ../Environ/Doc/INPUT_Environ.html . ; \
	ln -fs ../Environ/Doc/INPUT_Environ.txt . ; \
	ln -fs ../Environ/Doc/INPUT_Environ.xml .)
