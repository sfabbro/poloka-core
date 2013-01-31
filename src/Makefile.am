CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* \
	dbfileLex.c \
	dbfileParse.c \
	dbfileParse.h

AM_YFLAGS = -l -d

lib_LTLIBRARIES = libpoloka-core.la

EXTRA_DIST = virtualinstrument.cc virtualinstrument.h slasub.c
BUILT_SOURCES = alltelinst.cc alltelinst.h

noinst_HEADERS =  docmain.h
src_includedir = $(includedir)/poloka
src_include_HEADERS = \
	allreducedimage.h \
	apersestar.h \
	array2d.h \
	astroutils.h \
	basestar.h \
	ccd.h \
	cluster2.h \
	cluster.h \
	convolution.h \
	convtemplate.h \
	countedref.h \
	datacards.h \
	dbconfigexception.h \
	dbimage.h \
	detection.h \
	dicstar.h \
	dictfile.h \
	dimage.h \
	exceptions.h \
	fastfinder.h \
	fastifstream.h \
	fatpoint.h \
	fileutils.h \
	fitsexception.h \
	fitsimage.h \
	fitsslice.h \
	fitstoad.h \
	frame.h \
	gaussianfit.h \
	globalval.h \
	gtransfo.h \
	histo1d.h \
	histo2d.h \
	histo4d.h \
	histopeakfinder.h \
	hprofile1d.h \
	imageback.h \
	imagebinning.h \
	image.h \
	imageinterpolation.h \
	imagelist.h \
	imagematch.h \
	imagepair.h \
	imagesubtraction.h \
	imagesum.h \
	imageutils.h \
	intframe.h \
	intpoint.h \
	iohelpers.h \
	kernelfit.h \
	kernelfitter.h \
	listmatch.h \
	matchexception.h \
	matvect.h \
	nstarmatch.h \
	photoratio.h \
	pixelblock.h \
	point.h \
	polokaexception.h \
	poly1.h \
	poly2.h \
	quali_box.h \
	reducedimage.h \
	reducedutils.h \
	scorescollection.h \
	seeing_box.h \
	sestar.h \
	sextractor_box.h \
	simps2d.h \
	starlistexception.h \
	starlist.cc \
	starlist.h \
	starmatch.h \
	stringlist.h \
	sub.h \
	subimage.h \
	swarpstack.h \
	toadscards.h \
	transformedimage.h \
	usnoutils.h \
	vutils.h \
	wcscon.h \
	wcs_defines.h \
	wcsutils.h

libpoloka_core_la_SOURCES = \
	$(src_include_HEADERS) \
	dbfileParse.y \
	dbfileLex.l \
	allreducedimage.cc \
	apersestar.cc \
	astroutils.cc \
	basestar.cc \
	cluster.cc \
	convolution.cc \
	datacards.cc \
	dbimage.cc \
	detection.cc \
	dicstar.cc \
	dictfile.cc \
	dimage.cc \
	ellipticaper.cc \
	fastfinder.cc \
	fastifstream.cc \
	fileutils.cc \
	fitsimage.cc \
	fitsslice.cc \
	fitstoad.cc \
	frame.cc \
	gaussianfit.cc \
	globalval.cc \
	gtransfo.cc \
	gzstream.cc \
	histo1d.cc \
	histo2d.cc \
	histo4d.cc \
	histopeakfinder.cc \
	hprofile1d.cc \
	image.cc \
	imageback.cc \
	imageinterpolation.cc \
	imagematch.cc \
	imagepair.cc \
	imagesubtraction.cc \
	imagesum.cc \
	imageutils.cc \
	kernelfit.cc \
	kernelfitter.cc \
	listmatch.cc \
	matvect.cc \
	multimagstar.cc \
	nstarmatch.cc \
	poly1.cc \
	poly2.cc \
	quali_box.cc \
	photoratio.cc \
	reducedimage.cc \
	reducedutils.cc \
	scorescollection.cc \
	seeing_box.cc \
	sestar.cc \
	sextractor_box.cc \
	signalhandlers.cc \
	starmatch.cc \
	stringlist.cc \
	sub.cc \
	subimage.cc \
	swarpstack.cc \
	toadscards.cc \
	transformedimage.cc \
	usnoutils.cc \
	vutils.cc \
	wcscon.c \
	wcsutils.cc

libpoloka_core_la_CPPFLAGS = \
	@CFITSIO_CFLAGS@ @SEXTRACTOR_CFLAGS@ @ZLIB_CFLAGS@

libpoloka_core_la_LIBADD = \
	@CFITSIO_LIBS@ @LAPACK_LIBS@ @SEXTRACTOR_LIBS@ @ZLIB_LIBS@

telinst_sources = $(top_srcdir)/telinst/*.cc
classname = \(.*class *\)\([^ :]*\)\(.*\)

alltelinst.cc: Makefile $(telinst_sources)
	@echo "// This file has been automatically generated; do not edit" > $@.new
	@echo "/* Source files were:" $(telinst_sources) "*/" >> $@.new
	@for telinst in $(telinst_sources); do echo "#include \"$$telinst\""; done >> $@.new
	@echo "#ifndef USE_WCS" >> $@.new
	@echo "static AcceptorType AcceptorsArray [] = {" >> $@.new
	@$(GREP) TYPE_SELECTOR $(telinst_sources) | $(SED) -e 's/$(classname)/  \&\2::Acceptor,/' >> $@.new
	@echo "  &Unknown::Acceptor};" >> $@.new
	@echo "// IsOfKind template instantiation" >> $@.new
	@$(GREP) TYPE_TESTER $(telinst_sources) | $(SED) \
		-e 's/$(classname)/template bool IsOfKind< \2 >(const FitsHeader \&);/' \
		-e 's/$(classname)/template bool IsOfKind< \2 >(const string \&);/' >> $@.new
	@echo "#endif // USE_WCS" >> $@.new
	@diff $@.new $@ &> /dev/null || mv -f $@.new $@

alltelinst.h: alltelinst.cc
	@echo "// This file has been automatically generated; do not edit" > $@.new
	@echo "/* Source files were:" $(telinst_sources) "*/" >> $@.new
	@echo "#ifndef ALLTELINST__H" >> $@.new
	@echo "#define ALLTELINST__H" >> $@.new
	@$(GREP) TYPE_TESTER $(telinst_sources) | $(SED) -e 's/$(classname)/class \2;/' >> $@.new
	@echo "#endif" >> $@.new
	@diff $@.new $@ &> /dev/null || mv -f $@.new $@
