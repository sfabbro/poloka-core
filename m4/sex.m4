# -*- autoconf -*-
# 
# $Id: sex.m4,v 1.1 2004/02/23 00:39:38 nrl Exp $
# 
# autoconf macro to check the sextractor installation
# Nicolas Regnault <regnault@in2p3.fr> Feb. 2004.
# 
AC_DEFUN([AM_CHECK_SEX],[

 sex_prefix=""
 sex_headers=""
 sex_libs=""

 AC_ARG_WITH(sex,
  [  --with-sex=<prefix>         prefix where sextractor is installed],
  [  sex_prefix="$withval"
     sex_headers="$withval/include"
     sex_libs="$withval/lib"],
  )

 AC_ARG_WITH(sex-headers,
  [  --with-sex-headers=<prefix> prefix where the sextractor headers are installed],
  [  sex_headers=$withval],
  )

 AC_ARG_WITH(sex-lib,
  [  --with-sex-lib=<prefix>     prefix where the sex lib is installed],
  [  sex_lib=$withval],
  )
 

 SEX_CPPFLAGS=""
 SEX_LDFLAGS=""
 CPPFLAGS_sav="$CPPFLAGS"
 LDFLAGS_sav="$LDFLAGS"

 sex_vreq=[$1]
 sex_major=`echo $sex_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
 sex_minor=`echo $sex_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`

 if test -n "$sex_headers" ; then
  SEX_CPPFLAGS="-I$sex_headers"
 elif test -n "$FROGSHOME" ; then
  SEX_CPPFLAGS="-I$FROGSHOME/include"
 elif test -n "$prefix" && test "$prefix" != "NONE"; then 
  SEX_CPPFLAGS="-I$prefix/include"
 fi
 CPPFLAGS="$CPPFLAGS $SEX_CPPFLAGS"

 AC_CHECK_HEADER(fitsio.h,,
  [
   echo "*** Header file fitsio.h not found.                                 "
   echo "*** If sextractor is not installed on your system, install it first."
   echo "*** You may retrieve it from:                                       "
   echo "***    ftp://??????????????????????                                 "
   echo "***                                                                 "
   echo "*** Keep in mind that ./configure looks for sextractor in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$FROGSHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install sextractor someplace else, you will have"
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-sex=                                                    "
   echo "***  --with-sex-headers=                                            "
   echo "***  --with-sex-libs=                                               "
   echo "***                                                                 "
   echo "*** Please check your sextractor installation and try again.        "
   AC_MSG_ERROR(aborting.)
  ])

 if test -n "$sex_libs" ; then
  SEX_LDFLAGS="-L$sex_libs -lsex -lwcs"
 elif test -n "$FROGSHOME" ; then
  SEX_LDFLAGS="-L$FROGSHOME/lib -lsex -lwcs"
 elif test -n "$prefix" && test "$prefix" != "NONE" ; then
  SEX_LDFLAGS="-L$prefix/lib -lsex -lwcs"
 fi
 LDFLAGS="$LDFLAGS $SEX_LDFLAGS"
 AC_CHECK_LIB(sex,main,,
  [echo "*** Library libsex not found.                                   "
   echo "*** If sextractor is not installed on your system, install it first."
   echo "*** You may retrieve it from:                                       "
   echo "***    ftp://??????                                                 "
   echo "***                                                                 "
   echo "*** Keep in mind that ./configure looks for sextrator in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$FROGSHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install sextrator someplace else, you will have"
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-sex=                                                    "
   echo "***  --with-sex-headers=                                            "
   echo "***  --with-sex-libs=                                               "
   echo "***                                                                 "
   echo "*** Please check your sextrator installation and try again.           "
  AC_MSG_ERROR(aborting.)
  ])

 # FIXME: check the version


 CPPFLAGS="$CPPFLAGS_sav"
 LDFLAGS="$LDFLAGS_sav"
 AC_SUBST(SEX_CPPFLAGS)
 AC_SUBST(SEX_LDFLAGS)
]
)
