# -*- autoconf -*-
# 
# $Id: lapack.m4,v 1.2 2004/02/23 16:32:22 guy Exp $
# 
# autoconf macro to check the lapack installation
# Nicolas Regnault <regnault@in2p3.fr> Feb. 2004.
# 

AC_DEFUN([CHECK_LAPACK],[

 lapack_prefix=""
 lapack_headers=""
 lapack_libs=""

 AC_ARG_WITH(lapack,
  [  --with-lapack=<prefix>         prefix where lapack is installed],
  [  lapack_prefix="$withval"
     lapack_headers="$withval/include"
     lapack_libs="$withval/lib"],
  )

 AC_ARG_WITH(lapack-headers,
  [  --with-lapack-headers=<prefix> prefix where the headers are installed],
  [  lapack_headers=$withval],
  )

 AC_ARG_WITH(lapack-lib,
  [  --with-lapack-lib=<prefix>     prefix where the lapack lib is installed],
  [  lapack_libs=$withval],
  )
 

 LAPACK_CPPFLAGS=""
 LAPACK_LDFLAGS=""
 CPPFLAGS_sav="$CPPFLAGS"
 LDFLAGS_sav="$LDFLAGS"
 LIBS_sav="$LIBS"

 lapack_vreq=[$1]
 lapack_major=`echo $lapack_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
 lapack_minor=`echo $lapack_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`

 if test -n "$lapack_headers" ; then
  LAPACK_CPPFLAGS="-I$lapack_headers"
 elif test -n "$LAPACKHOME" ; then
  LAPACK_CPPFLAGS="-I$LAPACKHOME/include"
 elif test -n "$prefix" && test "$prefix" != "NONE"; then 
  LAPACK_CPPFLAGS="-I$prefix/include"
 fi
 CPPFLAGS="$CPPFLAGS $LAPACK_CPPFLAGS"

 AC_CHECK_HEADER(lafnames.h,,
  [
   echo "*** Header file lafnames.h not found.                                 "
   echo "*** Keep in mind that ./configure looks for lapack in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$LAPACKHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install lapack someplace else, you will have  "
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-lapack=                                                "
   echo "***  --with-lapack-headers=                                        "
   echo "***  --with-lapack-libs=                                           "
   echo "***                                                                 "
   echo "*** Please check your lapack installation and try again.           "
   AC_MSG_ERROR(aborting.)
  ])

 if test -n "$lapack_libs" ; then
  LAPACK_LDFLAGS="-L$lapack_libs -llapack++ $lapack_libs/liblapack.a $lapack_libs/libblas.a  "
 elif test -n "$LAPACKHOME" ; then
  LAPACK_LDFLAGS="-L$LAPACKHOME/lib -llapack++ /home/guy/software/lapack/v1.1a/lib/liblapack.a $LAPACKHOME/lib/libblas.a "
 elif test -n "$prefix" && test "$prefix" != "NONE" ; then
  LAPACK_LDFLAGS="-L$prefix/lib  -llapack++ $prefix/lib/liblapack.a $prefix/lib/libblas.a  "
 fi
 LDFLAGS="$LDFLAGS $LAPACK_LDFLAGS -lm"
 AC_CHECK_LIB(lapack++,main,,
  [echo "*** Unable to link libraries  libblas, liblapack and liblapack++ "
   echo "*** see config.log for details                                    "
   echo "***                                                                "
   echo "*** Keep in mind that ./configure looks for lapack in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$LAPACKHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install lapack someplace else, you will have  "
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-lapack=                                                "
   echo "***  --with-lapack-headers=                                        "
   echo "***  --with-lapack-libs=                                           "
   echo "***                                                                 "
   echo "*** Please check your lapack installation and try again.           "
  AC_MSG_ERROR(aborting.)
  ])

 CPPFLAGS="$CPPFLAGS_sav"
 LDFLAGS="$LDFLAGS_sav"
 LIBS="$LIBS_sav"
 AC_SUBST(LAPACK_CPPFLAGS)
 AC_SUBST(LAPACK_LDFLAGS)
]
)
