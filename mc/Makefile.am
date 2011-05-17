CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*

lib_LTLIBRARIES = libtoadsmc.la

mc_includedir = $(includedir)/poloka
mc_include_HEADERS = \
	addlcfakes.h \
	agaussian.h \
	dual_sub.h \
	fakelist.h \
	fakeobject.h \
	fakesnia.h \
	host.h \
	imagemanager.h \
	lcsim.h \
	mcimage.h \
	mcsub.h \
	myrdm.h \
	simsnstar.h \
	simulation.h

libtoadsmc_la_SOURCES = \
	$(mc_include_HEADERS) \
	addlcfakes.cc \
	agaussian.cc \
	dual_sub.cc \
	fakelist.cc \
	fakeobject.cc \
	fakesnia.cc \
	host.cc \
	imagemanager.cc \
	lcsim.cc \
	mcimage.cc \
	mcsub.cc \
	myrdm.cc \
	simsnstar.cc \
	simulation.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/psf \
	@SNFIT_CFLAGS@ @CFITSIO_CFLAGS@

libtoadsmc_la_LIBADD = \
	$(top_builddir)/src/libtoads.la \
	$(top_builddir)/psf/libtoadspsf.la \
	@SNFIT_LIBS@ @CFITSIO_LIBS@

LDADD = libtoadsmc.la
bin_PROGRAMS = genlcfakes make_lcsim mc_addsn make_dualsub make_mc_cat

genlcfakes_SOURCES = genlcfakes.cc
make_lcsim_SOURCES = make_lcsim.cc
mc_addsn_SOURCES = mc_addsn.cc
make_dualsub_SOURCES = make_dualsub.cc
make_mc_cat_SOURCES = make_mc_cat.cc
