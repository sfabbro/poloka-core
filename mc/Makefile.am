AM_CPPFLAGS = -Wall -I../src -I../objio -I../dao @XML_CFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ -DFNAME
AM_LDFLAGS  = -L../src -L../dao -L../objio -L@libdir@ -ltoadio @XML_LIBS@ CFITSIO_LDFLAGS@ @CMT_LDFLAGS@

%_dict.cc : %.h
	../objio/dictgen -o $@ $<


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3



####   matvec.cc_cppflags = $(optimize_it)

lib_LTLIBRARIES = libtoadsmc.la
libtoadsmc_la_SOURCES = \
simsnstar_dict.cc \
simulation_dict.cc \
simsnstar.cc simsnstar.h \
simulation.cc simulation.h \
myrdm.cc myrdm.h \
agaussian.cc agaussian.h \
mcimage.cc mcimage.h \
mcsub.cc mcsub.h 

include_HEADERS = \
agaussian.h \
simsnstar.h \
simulation.h \
myrdm.h \
mcimage.h \
mcsub.h


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
