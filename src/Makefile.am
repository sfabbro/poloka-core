CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* \
	dbfileLex.c \
	dbfileParse.c \
	dbfileParse.h

AM_YFLAGS = -l -d

lib_LTLIBRARIES = libpoloka-core.la

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
	fitsimage.h \
	fitsslice.h \
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
	virtualinstrument.h \
	vutils.h \
	wcscon.h \
	wcs_defines.h \
	wcsutils.h

libpoloka-core_la_SOURCES = \
	$(src_include_HEADERS) \
	dbfileParse.y \
	dbfileLex.l \
	allreducedimage.cc \
	apersestar.cc \
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
	listmatch.cc \
	multimagstar.cc \
	nstarmatch.cc \
	quali_box.cc \
	photoratio.cc \
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

libpoloka-core_la_CPPFLAGS = \
	@CFITSIO_CFLAGS@ @SEXTRACTOR_CFLAGS@ @Z_CFLAGS@

libpoloka-core_la_LIBADD = \
	@CFITSIO_LIBS@ @LAPACK_LIBS@ @SEXTRACTOR_LIBS@ @Z_LIBS@
