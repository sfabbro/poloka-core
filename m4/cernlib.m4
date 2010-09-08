# -*- autoconf -*-
# 
# $Id: cernlib.m4,v 1.9 2010/09/08 22:09:11 seb Exp $
# 
# autoconf macro to check the cernlib installation
# Nicolas Regnault <regnault@in2p3.fr> Feb. 2004.
# 
AC_DEFUN([CHECK_CERNLIB],[
AC_PATH_PROG(WITH_CERNLIB, cernlib, no)
if test "$WITH_CERNLIB" = "no" ; then
    AC_MSG_ERROR([Cannot find program cernlib used to find link options. You may retrieve it from http://wwwasd.web.cern.ch/wwwasd/cernlib/version.html. If you do not want to use it, try configure --disable-cernlib.],)
fi
CERNLIB_COMMAND=$WITH_CERNLIB
AC_SUBST(CERNLIB_COMMAND)


])


