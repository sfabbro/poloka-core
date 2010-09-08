## this library is not compiled if sextractor is disable ( --disable-sex in configure)

if USESEXTRACTOR
#SNFIT_INC = -I/afs/in2p3.fr/throng/snovae/softsnls/snfit/v1r0/src

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src -I$(top_srcdir)/src_utils  $(SNFIT_INC) -I$(top_srcdir)/psf -I$(top_srcdir)/dao  -I$(top_builddir)/src_base -I$(top_builddir)/src -I$(top_builddir)/dao @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ -DFNAME
AM_LDFLAGS  = -L$(top_builddir)/src_base -L$(top_builddir)/src -L$(top_builddir)/dao -L$(top_srcdir)/mc -ltoadspsf -ltoadsbase -ltoads -L@libdir@ CFITSIO_LDFLAGS@ @CMT_LDFLAGS@


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3



####   matvec.cc_cppflags = $(optimize_it)

lib_LTLIBRARIES = libtoadsmc.la
libtoadsmc_la_SOURCES = \
#simsnstar_dict.cc \
#simulation_dict.cc \
simsnstar.cc simsnstar.h \
simulation.cc simulation.h \
myrdm.cc myrdm.h \
agaussian.cc agaussian.h \
mcimage.cc mcimage.h \
mcsub.cc mcsub.h \
dual_sub.cc dual_sub.h \
addlcfakes.cc \
fakelist.cc fakelist.h \
fakeobject.cc fakeobject.h \
imagemanager.cc imagemanager.h \
host.cc host.h

if ENABLE_MC
libtoadsmc_la_SOURCES += lcsim.cc lcsim.h fakesnia.cc fakesnia.h 
endif

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
host.h

if ENABLE_MC
include_HEADERS += lcsim.h
endif



CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

@prefix@/include/% : @abs_builddir@/%
	ln -fs @abs_builddir@/$* @prefix@/include/$*

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadsmc.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif

endif

