
AM_CPPFLAGS = -Wall -I../src -I../objio -I../dao @XML_CFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS  = -L../src -L../dao -L../objio -L@libdir@ -ltoadio @XML_LIBS@ CFITSIO_LDFLAGS@ @CMT_LDFLAGS@


lib_LTLIBRARIES = libtoadslc.la

libtoadslc_la_SOURCES =  \
	gausspsf.cc	 \
	photstar.cc      \
	vignet.cc	 \
	vignetphot.cc	 \
	refstar.cc       \
	lcio.cc	         \
	lightcurve.cc	 \
	simfitvignet.cc  \
	simfit.cc        \
	simfitphot.cc 	 

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadslc.*))

endif