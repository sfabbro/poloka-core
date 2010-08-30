CLEANFILES = *~ alltelinst.*
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* *.new

EXTRA_DIST = virtualinstrument.cc virtualinstrument.h slasub.c

include_HEADERS = \
	alltelinst.h \
	astroutils.h \
	fileutils.h \
	fitsexception.h \
	fitsimage.h \
	fitstoad.h \
	frame.h \
	gzstream.h \
	image.h \
	point.h \
	polokaexception.h \
	stringlist.h \
	vutils.h

noinst_HEADERS = wcs_defines.h wcscon.h

lib_LTLIBRARIES = libtoadsbase.la
libtoadsbase_la_SOURCES = \
	$(include_HEADERS) \
	$(noinst_HEADERS) \
	astroutils.cc \
	fileutils.cc \
	fitsimage.cc \
	fitstoad.cc \
	frame.cc \
	gzstream.cc \
	image.cc \
	signalhandlers.cc \
	stringlist.cc \
	trapfpe.c \
	vutils.cc \
	wcscon.c

libtoadsbase_la_CPPFLAGS = @CFITSIO_CFLAGS@
libtoadsbase_la_LIBADD =  @CFITSIO_LIBS@

telinst_sources = $(top_srcdir)/telinst/*.cc
classname = \(.*class *\)\([^ :]*\)\(.*\)

alltelinst.cc: Makefile $(telinst_sources)
	@echo "// This file has been automatically generated; do not edit" > $@.new
	@echo "/* Source files were:" $(telinst_sources) "*/" >> $@.new
	@for telinst in $(telinst_sources); do echo "#include \"$$telinst\""; done >> $@.new
	@echo "#ifndef USE_WCS" >> $@.new
	@echo "static AcceptorType AcceptorsArray [] = {" >> $@.new
	@grep TYPE_SELECTOR $(telinst_sources) | $(SED) -e 's/$(classname)/  \&\2::Acceptor,/' >> $@.new
	@echo "  &Unknown::Acceptor};" >> $@.new
	@echo "// IsOfKind template instantiation" >> $@.new
	@grep TYPE_TESTER $(telinst_sources) | $(SED) \
	 	-e 's/$(classname)/template bool IsOfKind< \2 >(const FitsHeader \&);/' \
		-e 's/$(classname)/template bool IsOfKind< \2 >(const string \&);/' >> $@.new
	@echo "#endif // USE_WCS" >> $@.new
	@diff $@.new $@ &> /dev/null || mv -f $@.new $@

alltelinst.h: alltelinst.cc
	@echo "// This file has been automatically generated; do not edit" > $@.new
	@echo "/* Source files were:" $(telinst_sources) "*/" >> $@.new
	@echo "#ifndef ALLTELINST__H" >> $@.new
	@echo "#define ALLTELINST__H" >> $@.new
	@grep TYPE_TESTER $(telinst_sources) | $(SED) -e 's/$(classname)/class \2;/' >> $@.new
	@echo "#endif" >> $@.new
	@diff $@.new $@ &> /dev/null || mv -f $@.new $@

# hack to produce alltelinst.h before the rest. (fitstoad.h should be before any sources)
fitsimage.cc: alltelinst.h
