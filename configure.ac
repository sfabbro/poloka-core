## Initialization 
AC_PREREQ(2.59)
AC_INIT(poloka, 0.3.0, [poloka@in2p3.fr])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src_base/fitsimage.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(autobuild)
AM_INIT_AUTOMAKE

## Keep this one before libtool check
AM_DISABLE_STATIC

## Checks for programs to build poloka libraries
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_LEX
AC_PROG_LIBTOOL
# AC_PROG_SED is available only in autoconf 2.60 and later.
# With older autoconf, assume that 'sed' is ok.
ifdef([AC_PROG_SED], [AC_PROG_SED], [SED=sed AC_SUBST(SED)])
ifdef([AC_PROG_GREP], [AC_PROG_GREP], [GREP=grep AC_SUBST(GREP)])
AC_PROG_YACC

## Check for mandatory blas/lapack
AX_LAPACK(, [AC_MSG_ERROR([blas or lapack was not found. 
*** Please get blas and lapack installed, or use the options 
*** --with-blas or --with-lapack.])
])

## Check for mandatory cfitsio
AX_CFITSIO(, [AC_MSG_ERROR([cfitsio library could not be detected
*** Please install cfitsio or use the options:
*** --with-cfitsio-includedir and --with-cfitsio-libdir.])
])

## Check for mandatory SExtractor repacked library
AX_LIBSEX(, [AC_MSG_ERROR([SExtractor library could not be detected 
*** Please install libsex or use the options:
*** --with-libsex-includedir and --with-libsex-libdir.])
])

## Check for mandatory zlib (for gztream)
AX_CHECK_ZLIB

## Check for optional cernlib to build
CERNLIB=false
AC_MSG_CHECKING([whether to build cernlib tools])
AC_ARG_ENABLE(cernlib,
    AC_HELP_STRING([--enable-cernlib], [cernlib support (default=no)]),
   [enable_cernlib=$enableval],
   [enable_cernlib=no])
if test x"$enable_cernlib" = x"yes"; then
    AC_MSG_RESULT([yes])
    AC_LIB_CERNLIB([packlib], [HROPEN],
    		   [CERNLIB=true],
		   [AC_MSG_ERROR([packlib was not found.
		   *** Please install cernlib (packlib) in order to build the cernlib tools.])])
    CERNLIB_PACKLIBS="$CERNLIB_LIBS"
    AC_SUBST(CERNLIB_PACKLIBS)
else
    AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(CERNLIB, test x$CERNLIB = xtrue)

## Check for optional snfit to build mc module
MC=false
AC_MSG_CHECKING([whether to build the Monte Carlo module])
AC_ARG_ENABLE([mc],
   AS_HELP_STRING([--enable-mc], [build montecarlo module (default=no)]),
   [enable_mc=$enableval],
   [enable_mc=no])
if test x"$enable_mc" = "xyes"; then
    AC_MSG_RESULT([yes])
    AX_SNFIT([MC=true], [AC_MSG_ERROR([snfit library was not found. 
   *** Please install snfit or use the options:
   *** --with-snfit-includedir and --with-snfit-libdir.])
   ])
else
    AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(MC, test x$MC = xtrue)

## Check for optional imagemagick to output some vignettes
IMAGEMAGICK=false
AC_MSG_CHECKING([whether to build the imagemagick utilities])
AC_ARG_ENABLE([imagemagick],
      [AS_HELP_STRING([--enable-imagemagick], [build imagemagick tools (default=no)])],
      [enable_imagemagick=$enableval], 		   	   
      [enable_imagemagick=no])
if test x"$enable_imagemagick" = x"yes"; then
    AC_MSG_RESULT([yes])
    AX_IMAGEMAGICK([IMAGEMAGICK=true], 
		   [AC_MSG_ERROR([ImageMagick++ library was not found. *** Please install ImageMagick.])])
else
    AC_MSG_RESULT([no])
fi
AM_CONDITIONAL([IMAGEMAGICK], [test x$IMAGEMAGICK = xtrue])

AC_CONFIG_FILES([Makefile
		 poloka.pc
                 cern_stuff/Makefile
                 cern_utils/Makefile
                 datacards/Makefile
                 doc/Doxyfile
                 doc/Makefile
		 flat/Makefile
		 imagemagick_utils/Makefile
  		 lc/Makefile
  		 mc/Makefile
                 psf/Makefile
                 simphot/Makefile
                 src/Makefile
		 src_base/Makefile
		 src_utils/Makefile
                 telinst/Makefile
                 utils/Makefile])
AC_OUTPUT
