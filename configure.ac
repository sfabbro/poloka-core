#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(repack2, 0.1.0, regnault@in2p3.fr)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST
AM_DISABLE_STATIC

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
AC_CHECK_PROG(DOXYGEN,doxygen)

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcernstuff':
#AC_CHECK_LIB([cernstuff], [main])
# FIXME: Replace `main' with a function in `-lcrypt':
AC_CHECK_LIB([crypt], [main])
# FIXME: Replace `main' with a function in `-ldl':
AC_CHECK_LIB([dl], [main])
# FIXME: Replace `main' with a function in `-lg2c':
AC_CHECK_LIB([g2c], [main])
# FIXME: Replace `main' with a function in `-lnsl':
AC_CHECK_LIB([nsl], [main])
# FIXME: Replace `main' with a function in `-lpacklib':
#AC_CHECK_LIB([packlib], [main])
# FIXME: Replace `main' with a function in `-ltoadio':
#AC_CHECK_LIB([toadio], [main])
# FIXME: Replace `main' with a function in `-ltoads':
#AC_CHECK_LIB([toads], [main])
# FIXME: Replace `main' with a function in `-lxml2':
#AC_CHECK_LIB([xml2], [main])

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

AC_CONFIG_FILES([Makefile
                 bindings/Makefile
                 cern_stuff/Makefile
                 dao_stuff/Makefile
                 datacards/Makefile
                 dbutils/Makefile
                 doc/Makefile
                 doc/Doxyfile
                 objio/Makefile
                 src/Makefile
                 telinst/Makefile
                 utils/Makefile])
AC_OUTPUT
