EXTRA_DIST = README.db refman.ps Doxyfile.in

docdir = $(datadir)/doc
doc_DATA = html latex man
DOX = @DOXYGEN@
HEADER_LIST = @top_srcdir@/*/*.cc @top_srcdir@/*/*.h


all: stamp-doc

html: stamp-doc

clean :	clean-local

stamp-doc: Doxyfile $(HEADER_LIST)
	$(DOX) Doxyfile
	@echo > timestamp > $@

install-docDATA:
	$(mkinstalldirs) $(docdir);
	(cp -r html $(docdir); cd $(docdir)/html; rm -rf CVS)
	(cp -r man $(docdir);cd $(docdir)/man; rm -rf CVS)

uninstall-docDATA:
	cd $(docdir) ; rm -rf html
	cd $(docdir) ; rm -rf man
	rmdir $(docdir)

clean-local:
	@rm -rf html latex man rtf stamp-doc
	@mkdir html latex man rtf
