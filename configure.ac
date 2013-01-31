## Initialization 
AC_PREREQ(2.59)
AC_INIT(poloka, 0.8.7, [poloka team])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/fitsimage.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([foreign])

## Keep this one before libtool check
AM_DISABLE_STATIC

## Checks for programs to build poloka libraries
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_LEX
AC_PROG_F77
AC_PROG_LIBTOOL
AC_F77_LIBRARY_LDFLAGS

# AC_PROG_SED is available only in autoconf 2.60 and later.
# With older autoconf, assume that 'sed' is ok.
ifdef([AC_PROG_SED], [AC_PROG_SED], [SED=sed AC_SUBST(SED)])
ifdef([AC_PROG_GREP], [AC_PROG_GREP], [GREP=grep AC_SUBST(GREP)])
AC_PROG_YACC

## Check for mandatory lapack
PKG_CHECK_MODULES([LAPACK],
   		  [lapack],
		  [LAPACK_PC=lapack],
		  [AX_LAPACK([PC_LIBS="$PC_LIBS $LAPACK_LIBS"],
		  AC_MSG_ERROR([Cannot find lapack libraries]))])   
AC_SUBST([BLAS_PC])
AC_SUBST([LAPACK_PC])
AC_SUBST([PC_LIBS])

## Check for mandatory cfitsio
AX_CHECK_PKG_LIB(
	[cfitsio],
	[fitsio.h],
	[cfitsio],[ffopen],,
	[AC_MSG_ERROR([Could not find cfitsio library])])

## Check for mandatory SExtractor repacked library
AX_CHECK_PKG_LIB(
	[sextractor],
	[define.h],
	[sextractor],[makeit],,
	[AC_MSG_ERROR([Could not find SExtractor library])])

## Check for mandatory zlib (for gztream)
AX_CHECK_PKG_LIB(
	[zlib],
	[zlib.h],
	[z],[inflateEnd],,
	[AC_MSG_ERROR([Could not find zlib library])])


AC_CONFIG_FILES([Makefile
		 poloka-core.pc
                 datacards/Makefile
                 doc/Doxyfile
                 doc/Makefile
                 src/Makefile
                 telinst/Makefile
                 tools/Makefile])
AC_OUTPUT
