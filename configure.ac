#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

PACKAGE=repack2
VERSION=0.1.0

AC_PREREQ(2.57)
AC_INIT(repack2,0.1.0,regnault@in2p3.fr)
AC_CONFIG_SRCDIR([config.h.in])
#AM_CONFIG_HEADER([config.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST
#AM_DISABLE_STATIC

# default prefix for installation is in the current dir
AC_PREFIX_DEFAULT([`\pwd`/installarea])

# optional features

# Checks for programs.
AC_PROG_F77
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_PROG_YACC
AM_PROG_LEX
AC_PROG_MAKE_SET
AC_PATH_PROG(DOXYGEN,doxygen)
SWIG_PROG(1.3.20)
SWIG_ENABLE_CXX


CPPFLAGS="$CPPFLAGS -I/usr/local/include"

# Checks for libraries.
AM_PATH_PYTHON(2.0.0)
# FIXME: Replace `main' with a function in `-lcrypt':
AC_CHECK_LIB([crypt], [main])
# FIXME: Replace `main' with a function in `-ldl':
AC_CHECK_LIB([dl], [main])
# FIXME: Replace `main' with a function in `-lg2c':
AC_CHECK_LIB([g2c], [main])
# FIXME: Replace `main' with a function in `-lnsl':
AC_CHECK_LIB([nsl], [main])
AM_PATH_XML2(2.6.2)

# Use CMT to link packages, see cmt/configure.in
AC_ARG_ENABLE(cmt,
[AC_HELP_STRING([--enable-cmt],[use cmt to find external software (default is NO)])],
[CHECK_CMT],
[ # else put here the additionnal software to check
CHECK_CFITSIO(2.470)
CHECK_SEX
CHECK_LAPACK
])

# Check cernlib 
AC_ARG_ENABLE(cernlib,
[AC_HELP_STRING([--disable-cernlib],[Do not use cernlib, programs such as l2tup will no be compiled ])],,
[CHECK_CERNLIB])


# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE
AC_LANG(C)
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long double)

# Checks for library functions.
AC_FUNC_LSTAT
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor getcwd memchr memset pow rint setenv sqrt strchr strrchr strspn strstr])

CXXFLAGS="-g"
CPPFLAGS="-g"
CFLAGS="-g"
FFLAGS="-g"



AC_CONFIG_FILES([Makefile
                 bindings/Makefile
                 cern_stuff/Makefile
                 dao_stuff/Makefile
                 datacards/Makefile
                 doc/Makefile
                 doc/Doxyfile
                 objio/Makefile
                 src/Makefile
                 telinst/Makefile
                 utils/Makefile
		 cern_utils/Makefile
                 repack-config
                 repack.spec
		 cmt/requirements
		 cmt/version.cmt
])

AC_OUTPUT
