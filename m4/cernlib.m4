# -*- autoconf -*-
# 
# $Id: cernlib.m4,v 1.5 2004/02/23 16:32:48 guy Exp $
# 
# autoconf macro to check the cernlib installation
# Nicolas Regnault <regnault@in2p3.fr> Feb. 2004.
# 
AC_DEFUN([CHECK_CERNLIB],[

 cernlib_prefix=""

 AC_ARG_WITH(cernlib,
  [  --with-cernlib=<prefix>         prefix where cernlib is installed],
  [  cernlib_prefix="$withval"],
  )

 CERNLIB_LDFLAGS=""
 LDFLAGS_sav="$LDFLAGS"
 LIBS_sav="$LIBS"
 
 if test -n "$cernlib_prefix" ; then
  CERNLIB_LDFLAGS="-L$cernlib_prefix -lpacklib -lnsl -lcrypt -ldl"
 elif test -n "$CERN" ; then
  CERNLIB_LDFLAGS="-L$CERN/pro/lib -lpacklib -lnsl -lcrypt -ldl"
 elif test -n "$prefix" && test "$prefix" != "NONE" ; then
  CERNLIB_LDFLAGS="-L$prefix/lib -lpacklib -lnsl -lcrypt -ldl"
 fi
 LDFLAGS="$LDFLAGS $CERNLIB_LDFLAGS"


# FIXME: this code breaks the typechecks. Investigate.

 AC_CHECK_LIB(packlib,main,
  [
   toto=1
  ],
  [echo "*** Library libpacklib not found.                                   "
   echo "*** If the CERNLIB is not installed on your system, install it first."
   echo "*** You may retrieve it from:                                       "
   echo "***    http://wwwasd.web.cern.ch/wwwasd/cernlib/version.html        "
   echo "***                                                                 "
   echo "*** Keep in mind that ./configure looks for the CERNlib in the following"
   echo "*** standard locations:                                             "
   echo "***    /, /usr, /usr/local, \$CERN and <install prefix>             "
   echo "***                                                                 "
   echo "*** If you choose to install cernlib someplace else, you will have  "
   echo "*** to specify its install prefix using one or more of the following"
   echo "*** configure options:                                              "
   echo "***  --with-cernlib=                                                "
   echo "***                                                                 "
   echo "*** Please check your CERNlib installation and try again.           "
  AC_MSG_ERROR(aborting.)
  ])
 
 # TODO: check the cernlib version...

 LDFLAGS="$LDFLAGS_sav"
 LIBS="$LIBS_sav"
 AC_SUBST(CERNLIB_LDFLAGS)
]
)
