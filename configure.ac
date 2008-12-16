#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
 
PACKAGE=poloka
VERSION=0.1.0



AC_PREREQ(2.57)
AC_INIT(poloka,0.1.0,regnault@in2p3.fr)
AC_CONFIG_SRCDIR([config.h.in])
#AM_CONFIG_HEADER([config.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST
AM_DISABLE_STATIC



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



CPPFLAGS="$CPPFLAGS -I/usr/local/include"

# Checks for libraries.
### AM_PATH_PYTHON(2.0.0)


# FIXME: Replace `main' with a function in `-lcrypt':
AC_CHECK_LIB([crypt], [main])
# FIXME: Replace `main' with a function in `-ldl':
AC_CHECK_LIB([dl], [main])
# FIXME: Replace `main' with a function in `-lg2c':
AC_CHECK_LIB([g2c], [main])
# FIXME: Replace `main' with a function in `-lnsl':
AC_CHECK_LIB([nsl], [main])


# disable objio to remove swig and xml dependence 
#AC_ARG_ENABLE(objio,
#[AC_HELP_STRING([--disable-objio],[Do not use persistence, all programs using it will no be compiled, no need for xml nor swig ])],
#[objio=false],
#[ objio=true])
#AM_CONDITIONAL(OBJIOENABLE, test x$objio = xtrue)

#if test  x$objio = xtrue ; then
#AM_PATH_XML2(2.6.2)
#SWIG_PROG(1.3.20)
#SWIG_ENABLE_CXX
#fi

# Check sex 
AC_ARG_ENABLE(sex,
[AC_HELP_STRING([--disable-sex],[Do not use sextractor, all programs using it will no be compiled ])],
[sex=false
 # and set cmt flags to null
SEX_CFLAGS=""
SEX_LDFLAGS=""
AC_SUBST(SEX_CFLAGS)
AC_SUBST(SEX_LDFLAGS)
],
[ sex=true])
AM_CONDITIONAL(USESEXTRACTOR, test x$sex = xtrue)

# whether to compile the MC code or not
AC_ARG_ENABLE(mc,
[AC_HELP_STRING([--enable-mc],[compile the Monte-Carlo Code])],
[mc=true
 MC_CFLAGS=""
 MC_LDFLAGS=""
 AC_SUBST(MC_CFLAGS)
 AC_SUBST(MC_LDFLAGS)
],)
AM_CONDITIONAL(ENABLE_MC, test x$mc = xtrue)


# Use CMT to link packages, see cmt/configure.in
AC_ARG_ENABLE(cmt,
[AC_HELP_STRING([--disable-cmt],[do not use cmt to find external software (default is yes)])],
[ # else put here the additionnal software to check
if test $enable_cmt = yes; then
  use_cmt=true
else 
  use_cmt=false
fi
],
[
use_cmt=true
])

#may seem complicated, but AC_ARG_ENABLE(only checks if 
if test $use_cmt = true ; then
CHECK_CMT
else
CHECK_CFITSIO(2.400)
if test  x$sex = xtrue ; then
CHECK_SEX
 # and set cmt flags to null
CMT_CFLAGS=""
CMT_LDFLAGS=""
AC_SUBST(CMT_CFLAGS)
AC_SUBST(CMT_LDFLAGS)
fi
fi

# Use imagemagick to output some vignettes in gif
AC_ARG_ENABLE(imagemagick,
[AC_HELP_STRING([--enable-imagemagick],[compile some utils using imagemagick (default is NO)])],
[
imagemagick=true
CHECK_IMAGEMAGICK
],
[imagemagick=false])
AM_CONDITIONAL(IMAGEMAGICK, test x$imagemagick = xtrue)


######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package
AC_ARG_ENABLE(autoinstall,
[AC_HELP_STRING([--disable-autoinstall],[disable automatic installation using symbolic links])],
[
if test $enable_autoinstall = no; then
  autoinstall=false
else
  autoinstall=true
fi
],
[autoinstall=true])
AM_CONDITIONAL(AUTOINSTALL, test x$autoinstall != xfalse)

# Check cernlib 
AC_ARG_ENABLE(cernlib,
[AC_HELP_STRING([--disable-cernlib],[Do not use cernlib, programs such as l2tup will no be compiled ])],
[cernlib=false],
[ cernlib=true ; CHECK_CERNLIB ])
AM_CONDITIONAL(USECERNLIB, test x$cernlib = xtrue)

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

CXXFLAGS="-g -Wall"
CPPFLAGS="-g -Wall"
CFLAGS="-g -Wall"
FFLAGS="-g"


# optional features

AC_SUBST(CMTPACKAGE,"`\basename $PWD`")


AC_CONFIG_FILES([Makefile
		 cmt/requirements
		 cmt/version.cmt
                 cern_stuff/Makefile
                 lapack_stuff/Makefile
                 dao/Makefile
                 psf/Makefile
                 simphot/Makefile
                 datacards/Makefile
		 src_base/Makefile
		 src_utils/Makefile
		 flat/Makefile
                 src/Makefile
  		 lc/Makefile
	         mc/Makefile
                 telinst/Makefile
                 utils/Makefile
		 cern_utils/Makefile
		 imagemagick_utils/Makefile
                 doc/Makefile
                 doc/Doxyfile
])


my_top_dir=`dirname $0`
my_abs_top_dir=`cd $my_top_dir && pwd`
prefix=$my_abs_top_dir/installarea

AC_OUTPUT
