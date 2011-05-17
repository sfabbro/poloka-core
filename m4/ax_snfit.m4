# -*- autoconf -*-
# 
# use as AX_SNFIT([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#

AC_DEFUN([AX_SNFIT],[

ax_snfit_ok=no

PKG_CHECK_MODULES([SNFIT], [snfit],[ax_snfit_ok=yes], [ax_snfit_ok=no])

if test x"$ax_snfit_ok" = x"no"; then
   # if neither pkg-config file found, nor SNFIT_LIBS defined
   if test x$SNFIT_LIBS = x ; then
      SNFIT_LIBS="-lsnfit -lm"
      AC_ARG_WITH([snfit-libdir],
	 [  --with-snfit-libdir=DIR directory where the library was installed],
	 [SNFIT_LIBS="-L$withval $SNFIT_LIBS"], )
   fi
   LIBS_sav="$LIBS"
   LIBS="$LIBS $SNFIT_LIBS"
   ax_snfit_lib_ok=no
   ax_snfit_hdr_ok=no
   if test x$SNFIT_CFLAGS = x ; then
      SNFIT_CFLAGS=""
      AC_ARG_WITH([snfit-includedir],
	 [  --with-snfit-includedir=DIR directory where the headers were installed],
	 [SNFIT_CFLAGS="-I$withval"], )
   fi
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $SNFIT_CFLAGS"
   AC_LANG_PUSH([C++])
   AC_MSG_CHECKING(that snfit compiles properly)
   AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM(
	[[#include <param.h> ]],
	[[LcParam param();]])],	
     [AC_MSG_RESULT(yes); ax_snfit_lib_ok=yes],
     [AC_MSG_RESULT(no)])

   # if neither pkg-config file found, nor SNFIT_CFLAGS defined
   AC_CHECK_HEADER([salt2model.h],
      [ax_snfit_hdr_ok=yes],
      [AC_MSG_WARN([  *** snfit headers not found.])])
   AC_LANG_POP([C++])
   if test x$ax_snfit_lib_ok = xyes -a x$ax_snfit_hdr_ok = xyes; then
      ax_snfit_ok=yes
   fi
   AC_SUBST(SNFIT_CFLAGS)
   AC_SUBST(SNFIT_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_snfit_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_SNFIT, [1], [Define if you have SNFIT library.]),[$1])
        :
else
        ax_snfit_ok=no
        $2
fi
])
