EXTRA_DIST =

if USECERNLIB
 bin_PROGRAMS = l2tup 2im2tup im2tup 

 AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/dbutils -I$(top_srcdir)/objio @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
 AM_LDFLAGS = -L$(top_builddir)/src  -L$(top_builddir)/objio -L$(top_builddir)/cern_stuff -L$(top_builddir)/lapack_stuff -L@libdir@ -ltoads -ltoadio   -lcernstuff -llapackstuff @CMT_LDFLAGS@ @CERNLIB_LDFLAGS@ @XML_LIBS@ @SEX_LDFLAGS@  @CFITSIO_LDFLAGS@ -lg2c

 im2tup_SOURCES = im2tup.cc
 l2tup_SOURCES = l2tup.cc
 2im2tup_SOURCES = 2im2tup.cc


######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

all   : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})

endif

endif