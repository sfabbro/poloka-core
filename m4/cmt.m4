# -*- autoconf -*-
# 
# $Id: cmt.m4,v 1.4 2004/03/10 16:33:50 guy Exp $
# 
# autoconf macro to check  cmt
# 

####################################### 
AC_DEFUN([CHECK_CMT],[

# message
AC_MSG_CHECKING(cmt configuration)
AC_MSG_RESULT( )

# try to find cmt in PATH
AC_PATH_PROG(CMT_PROG, cmt, no)
if test "$CMT_PROG" = "no" ; then
    AC_MSG_ERROR(cannot find cmt)
fi

# check cmt version>=v1r14
CMT_VERSION_MAJOR_REQ="1"
CMT_VERSION_MINOR_REQ="14"

CMT_VERSION_MAJOR=`cmt version | sed 's/v/ /g' | sed 's/r/ /g' | sed 's/p/ /g' | awk '{printf("%d",$[1]);}'`
CMT_VERSION_MINOR=`cmt version | sed 's/v/ /g' | sed 's/r/ /g' | sed 's/p/ /g' | awk '{printf("%d",$[2]);}'`

AC_MSG_CHECKING(cmt version $CMT_VERSION_MAJOR.$CMT_VERSION_MINOR >= requested $CMT_VERSION_MAJOR_REQ.$CMT_VERSION_MINOR_REQ )


if test $CMT_VERSION_MAJOR -lt $CMT_VERSION_MAJOR_REQ ; then
	AC_MSG_ERROR(cmt version does not match the requirements)
elif test $CMT_VERSION_MAJOR = $CMT_VERSION_MAJOR_REQ ; then
	if test $CMT_VERSION_MINOR -lt $CMT_VERSION_MINOR_REQ ; then
		AC_MSG_ERROR(cmt version does not match the requirements)
	fi
fi
AC_MSG_RESULT(yes)



# check for missing packages in cmt
AC_MSG_CHECKING(cmt connected packages)

CMT_PACKAGE_NOT_FOUND=`cd cmt; cmt config | grep package_not_found`
if test -n "$CMT_PACKAGE_NOT_FOUND" ; then
   AC_MSG_ERROR($CMT_PACKAGE_NOT_FOUND)	
fi

# get includes and use_linkopts from cmt
CMT_CFLAGS=`cd cmt; cmt-config.csh --cflags`
CMT_LDFLAGS=`cd cmt; cmt-config.csh --ldflags | sed 's/-Wl,-rpath/-R/g'`

# save these values for usage in Makefile.am
AC_SUBST(CMT_CFLAGS)
AC_SUBST(CMT_LDFLAGS)

AC_MSG_RESULT(ok)
]
)
