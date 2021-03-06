## Initialization
AC_PREREQ([2.68])
AC_INIT([poloka-core],[0.2.0],[poloka team])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([poloka/fitsimage.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([foreign])

## Checks for programs to build poloka-core library
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LEX
AC_PROG_YACC

# AC_PROG_SED is available only in autoconf 2.60 and later.
# With older autoconf, assume that 'sed' is ok.
ifdef([AC_PROG_SED], [AC_PROG_SED], [SED=sed AC_SUBST(SED)])
ifdef([AC_PROG_GREP], [AC_PROG_GREP], [GREP=grep AC_SUBST(GREP)])

# F77 needed to link with blas/lapack
AC_PROG_F77
LT_INIT([disable-static])

## Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([inttypes.h libintl.h limits.h malloc.h stddef.h stdlib.h string.h unistd.h])

## Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

## Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_SEARCH_LIBS([cos], [m], [],
	       [AC_MSG_ERROR([Cannot find the cos() function])])


## Check for mandatory blas
PKG_CHECK_MODULES([BLAS],
   		  [blas],
		  [BLAS_PC=blas],
		  [AX_BLAS([PC_LIBS="$PC_LIBS $BLAS_LIBS"],
		  AC_MSG_ERROR([Cannot find blas libraries]))])

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
AX_PKG_LIB(
	[cfitsio],
	[fitsio.h],
	[cfitsio], [ffopen], [],
	[AC_MSG_ERROR([Cannot find cfitsio library])])

## Check for mandatory SExtractor repacked library
AX_PKG_LIB(
	[sextractor],
	[sextractor/define.h],
	[sextractor], [makeit], [],
	[AC_MSG_ERROR([Cannot find SExtractor library])])

AC_CONFIG_FILES([Makefile
		 poloka-core.pc
                 datacards/Makefile
                 doc/Doxyfile
                 doc/Makefile
		 poloka/Makefile
                 telinst/Makefile
		 tools/Makefile
                 ])
AC_OUTPUT
