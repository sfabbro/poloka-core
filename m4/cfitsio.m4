# -*- autoconf -*-
# 
# $Id: cfitsio.m4,v 1.5 2004/02/23 12:59:47 nrl Exp $
# 
# autoconf macro to check the cfitsio installation
# Nicolas Regnault <regnault@in2p3.fr> Feb. 2004.
# 
AC_DEFUN([CHECK_CFITSIO],[

 cfitsio_prefix=""
 cfitsio_headers=""
 cfitsio_libs=""

 AC_ARG_WITH(cfitsio,
  [  --with-cfitsio=<prefix>         prefix where cfitsio is installed],
  [  cfitsio_prefix="$withval"
     cfitsio_headers="$withval/include"
     cfitsio_libs="$withval/lib"],
  )

 AC_ARG_WITH(cfitsio-headers,
  [  --with-cfitsio-headers=<prefix> prefix where the headers are installed],
  [  cfitsio_headers=$withval],
  )

 AC_ARG_WITH(cfitsio-lib,
  [  --with-cfitsio-lib=<prefix>     prefix where the cfitsio lib is installed],
  [  cfitsio_lib=$withval],
  )
 

 CFITSIO_CPPFLAGS=""
 CFITSIO_LDFLAGS=""
 CPPFLAGS_sav="$CPPFLAGS"
 LDFLAGS_sav="$LDFLAGS"
 LIBS_sav="$LIBS"

 cfitsio_vreq=[$1]
 cfitsio_major=`echo $cfitsio_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
 cfitsio_minor=`echo $cfitsio_vreq | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`

 if test -n "$cfitsio_headers" ; then
  CFITSIO_CPPFLAGS="-I$cfitsio_headers"
 elif test -n "$CFITSIOHOME" ; then
  CFITSIO_CPPFLAGS="-I$CFITSIOHOME/include"
 elif test -n "$FROGSHOME" ; then
  CFITSIO_CPPFLAGS="-I$FROGSHOME/include"
 elif test -n "$prefix" && test "$prefix" != "NONE"; then 
  CFITSIO_CPPFLAGS="-I$prefix/include"
 fi
 CPPFLAGS="$CPPFLAGS $CFITSIO_CPPFLAGS"

 AC_CHECK_HEADER(fitsio.h,,
  [
   echo "*** Header file fitsio.h not found.                                 "
   echo "*** If cfitsio is not installed on your system, install it first.   "
   echo "*** You may retrieve it from:                                       "
   echo "***    ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/               "
   echo "***                                                                 "
   echo "*** Keep in mind that ./configure looks for cfitsio in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$FROGSHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install cfitsio someplace else, you will have  "
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-cfitsio=                                                "
   echo "***  --with-cfitsio-headers=                                        "
   echo "***  --with-cfitsio-libs=                                           "
   echo "***                                                                 "
   echo "*** Please check your cfitsio installation and try again.           "
   AC_MSG_ERROR(aborting.)
  ])

 if test -n "$cfitsio_libs" ; then
  CFITSIO_LDFLAGS="-L$cfitsio_libs -lcfitsio"
 elif test -n "$CFITSIOHOME" ; then
  CFITSIO_LDFLAGS="-L$CFITSIOHOME/lib -lcfitsio" 
 elif test -n "$FROGSHOME" ; then
  CFITSIO_LDFLAGS="-L$FROGSHOME/lib -lcfitsio"
 elif test -n "$prefix" && test "$prefix" != "NONE" ; then
  CFITSIO_LDFLAGS="-L$prefix/lib -lcfitsio"
 fi
 LDFLAGS="$LDFLAGS $CFITSIO_LDFLAGS -lm"
 AC_CHECK_LIB(cfitsio,main,,
  [echo "*** Library libcfitsio not found.                                   "
   echo "*** If cfitsio is not installed on your system, install it first.   "
   echo "*** You may retrieve it from:                                       "
   echo "***    ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/               "
   echo "***                                                                 "
   echo "*** Keep in mind that ./configure looks for cfitsio in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$FROGSHOME and <install prefix>        "
   echo "***                                                                 "
   echo "*** If you choose to install cfitsio someplace else, you will have  "
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-cfitsio=                                                "
   echo "***  --with-cfitsio-headers=                                        "
   echo "***  --with-cfitsio-libs=                                           "
   echo "***                                                                 "
   echo "*** Please check your cfitsio installation and try again.           "
  AC_MSG_ERROR(aborting.)
  ])


 # check the cfitsio version...
 # We have to run a program to do that, since
 # cfitsio does not export any simple configure script 
 AC_MSG_CHECKING(that cfitsio version >=$cfitsio_major.$cfitsio_minor)
 AC_RUN_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include <fitsio.h>
      #include <stdio.h>
      #include <stdlib.h>]],
    [[float v;
      int major, minor;
      int major_req=atoi("$cfitsio_major"); 
      int minor_req=atoi("$cfitsio_minor");
      v = ffvers(&v);
      major = (int)v;
      v -= major; v*=1000;
      minor = (int)v;
      
      if( (major>=major_req) || 
          (major==major_req && minor>=minor_req) )
	return 0;
      else {
       printf("*** found an old version of libcfitsio (%d.%d)\n", major, minor);
       printf("*** this software requires version (%d.%d) or higher\n", major_req, minor_req);
       printf("*** please download a newer version, for example from:\n");
       printf("***     ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/\n\n");
       printf("*** It may also be that you already have the required version\n");
       printf("*** But that ./configure cannot locate it. Please check your installation\n");
       printf("*** and use the --with-cfitsio-* options\n");
       return 1;
      }]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_ERROR(cfitsio version does not match the requirements)],
 )

 CPPFLAGS="$CPPFLAGS_sav"
 LDFLAGS="$LDFLAGS_sav"
 LIBS="$LIBS_sav"
 AC_SUBST(CFITSIO_CPPFLAGS)
 AC_SUBST(CFITSIO_LDFLAGS)
]
)
