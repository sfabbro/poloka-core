EXTRA_DIST =


if USECERNLIB

# if you miss CERNLIB_COMMAND, you probably did not update ../m4/cernlib.m4
 CERNLIB_LDFLAGS = `@CERNLIB_COMMAND@ packlib`

 bin_PROGRAMS = l2tup im2tup  2im2tup

 AM_CPPFLAGS = -I$(top_srcdir)/src_base -I$(top_srcdir)/dbutils @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
 AM_LDFLAGS = -L$(top_builddir)/src_base -R @libdir@ -L@libdir@ -ltoadsbase @CMT_LDFLAGS@ @CFITSIO_LDFLAGS@ ${CERNLIB_LDFLAGS}

 im2tup_SOURCES = im2tup.cc
 2im2tup_SOURCES = 2im2tup.cc
 l2tup_SOURCES = l2tup.cc

##  2im2tup_SOURCES = 2im2tup.cc


######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : %
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

all   : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})

endif

endif
