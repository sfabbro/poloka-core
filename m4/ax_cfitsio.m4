# -*- autoconf -*-
# 
# use as AX_CFITSIO([MIN VERSION[, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]]])
# ex: AX_CFITSIO(3.24)
#

AC_DEFUN([AX_CFITSIO],[

ax_cfitsio_ok=no
version_req=[$1]
version_check=
if test x$version_req != x; then
   version_check=">= $version_req"
fi

PKG_CHECK_MODULES([CFITSIO], [cfitsio $version_check], ax_cfitsio_ok=yes,[
    if test x$CFITSIO_LIBS = x ; then
       CFITSIO_LIBS="-lcfitsio -lm"
    fi

    AC_ARG_WITH(cfitsio,
	[  --with-cfitsio=DIR        directory where cfitsio source was compiled],
    	[CFITSIO_CFLAGS="-I$withval"
	 CFITSIO_LIBS="-L$withval $CFITSIO_LIBS"],
    )

    AC_ARG_WITH(cfitsio-includedir,
	[  --with-cfitsio-includedir=DIR directory where the headers were installed],
	[CFITSIO_CFLAGS=-I$withval],
    )

    AC_ARG_WITH(cfitsio-libdir,
	[  --with-cfitsio-libdir=DIR directory where the library was installed],
	[CFITSIO_LIBS=-L$withval $CFITSIO_LIBS],
    )

    CPPFLAGS_sav="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $CFITSIO_CFLAGS"

    AC_CHECK_HEADER(fitsio.h, ,
	[AC_MSG_WARN([  *** Header file fitsio.h not found.])])

    LIBS_sav="$LIBS"
    LIBS="$LIBS $CFITSIO_LIBS"

    AC_CHECK_LIB(cfitsio, ffopen, ax_cfitsio_ok=yes,
    	[AC_MSG_WARN([ *** Library libcfitsio not linked properly])])

 # check the cfitsio version...
 # We have to run a program to do that, since
 # cfitsio does not export any simple configure script 
 version_major=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
 version_minor=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`

 AC_MSG_CHECKING(that cfitsio version >=$version_major.$version_minor)
 AC_RUN_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include <fitsio.h>
      #include <stdio.h>
      #include <stdlib.h>]],
    [[float v;
      int major, minor;
      int major_req=atoi("$version_major"); 
      int minor_req=atoi("$version_minor");
      v = ffvers(&v);
      major = (int)v;
      v -= major; v*=1000;
      minor = (int)v;
      
      if( (major>=major_req) || 
          (major==major_req && minor>=minor_req) )
	return 0;
      else {
       printf("*** found an old version of cfitsio (%d.%d)\n", major, minor);
       printf("*** this software requires version (%d.%d) or higher\n", major_req, minor_req);
       return 1;
      }]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_ERROR(cfitsio version does not match the requirements)],
 )
   AC_SUBST(CFITSIO_CFLAGS)
   AC_SUBST(CFITSIO_LIBS)

   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_cfitsio_ok" = xyes; then
        ifelse([$2],,AC_DEFINE(HAVE_CFITSIO,1,[Define if you have CFITSIO library.]),[$2])
        :
else
        ax_cfitsio_ok=no
        $3
fi
])
