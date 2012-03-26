CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* \
	dbfileLex.c \
	dbfileParse.c \
	dbfileParse.h

AM_YFLAGS = -l -d

lib_LTLIBRARIES = libtoads.la

noinst_HEADERS =  docmain.h
src_includedir = $(includedir)/poloka
src_include_HEADERS = \
	allreducedimage.h \
	apersestar.h \
	basestar.h \
	ccd.h \
	cluster2.h \
	cluster.h \
	convolution.h \
	convtemplate.h \
	countedref.h \
	daophotio.h \
	daostar.h \
	datacards.h \
	dbconfigexception.h \
	dbimage.h \
	detection.h \
	dicstar.h \
	dictfile.h \
	dimage.h \
	ellipticaper.h \
	exceptions.h \
	fastfinder.h \
	fastifstream.h \
	fatpoint.h \
	gaussianfit.h \
	globalval.h \
	gtransfo.h \
	histopeakfinder.h \
	hprofile1d.h \
	imagebinning.h \
	imagelist.h \
	imagematch.h \
	imagepair.h \
	imagesubtraction.h \
	imagesum.h \
	imageutils.h \
	iohelpers.h \
	kernelfit.h \
	kernelfitter.h \
	lambert_dust.h \
	lambert_tag.h \
	listmatch.h \
	matchexception.h \
	multimagstar.h \
	nstarmatch.h \
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
	sub.h \
	subimage.h \
	swarpstack.h \
	toadscards.h \
	transformedimage.h \
	usnoutils.h \
	wcsutils.h

libtoads_la_SOURCES = \
	$(src_include_HEADERS) \
	dbfileParse.y \
	dbfileLex.l \
	allreducedimage.cc \
	apersestar.cc \
	basestar.cc \
	cluster.cc \
	convolution.cc \
	daophotio.cc \
	daostar.cc \
	datacards.cc \
	dbimage.cc \
	detection.cc \
	dicstar.cc \
	dictfile.cc \
	dimage.cc \
	ellipticaper.cc \
	fastfinder.cc \
	fastifstream.cc \
	gaussianfit.cc \
	globalval.cc \
	gtransfo.cc \
	histopeakfinder.cc \
	hprofile1d.cc \
	imagematch.cc \
	imagepair.cc \
	imagesubtraction.cc \
	imagesum.cc \
	imageutils.cc \
	kernelfit.cc \
	kernelfitter.cc \
	lambert_dust.cc \
	listmatch.cc \
	multimagstar.cc \
	nstarmatch.cc \
	quali_box.cc \
	reducedimage.cc \
	reducedutils.cc \
	scorescollection.cc \
	seeing_box.cc \
	sestar.cc \
	sextractor_box.cc \
	starmatch.cc \
	sub.cc \
	subimage.cc \
	swarpstack.cc \
	toadscards.cc \
	transformedimage.cc \
	usnoutils.cc \
	wcsutils.cc

# the $(top_builddir)/src_base is to include the makefile generated alltelinst.h
libtoads_la_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_builddir)/src_base \
	@CFITSIO_CFLAGS@ @SEXTRACTOR_CFLAGS@

libtoads_la_LIBADD = \
	$(top_builddir)/src_base/libtoadsbase.la \
	$(top_builddir)/src_utils/libtoadsutils.la \
	@CFITSIO_LIBS@ @LAPACK_LIBS@ @SEXTRACTOR_LIBS@
