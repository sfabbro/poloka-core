
lib_LTLIBRARIES = libtoadsmc.la

include_HEADERS = \
	agaussian.h \
	simsnstar.h \
	simulation.h \
	myrdm.h \
	mcimage.h \
	mcsub.h \
	dual_sub.h \
	fakelist.h \
	fakeobject.h \
	fakesnia.h \
	imagemanager.h \
	host.h \
	lcsim.h

libtoadsmc_la_SOURCES = \
	$(include_HEADERS) \
	simsnstar.cc \
	simulation.cc \
	myrdm.cc \
	agaussian.cc \
	mcimage.cc \
	mcsub.cc \
	dual_sub.cc \
	addlcfakes.cc \
	fakelist.cc \
	fakeobject.cc \
	fakesnia.cc \
	imagemanager.cc \
	host.cc \
	lcsim.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/psf \
	@SNFIT_CFLAGS@

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
