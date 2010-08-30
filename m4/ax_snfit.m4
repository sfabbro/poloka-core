# -*- autoconf -*-
# 
# use as AX_SNFIT([MIN VERSION[, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]]])
# all arguments optional
# ex: AX_SNFIT(2.2)
#

AC_DEFUN([AX_SNFIT],[

ax_snfit_ok=no
version_req=[$1]
version_check=
if test x$version_req != x; then
   version_check=">= $version_req"
fi

PKG_CHECK_MODULES([SNFIT], [snfit $version_check], ax_snfit_ok=yes,[
    if test x$SNFIT_LIBS = x ; then
       SNFIT_LIBS="-lsnfit -lm"
    fi

    AC_ARG_WITH(snfit,
	[  --with-snfit=DIR        directory where snfit source was compiled],
    	[SNFIT_CFLAGS="-I$withval"
	 SNFIT_LIBS="-L$withval $SNFIT_LIBS"],
    )

    AC_ARG_WITH(snfit-includedir,
	[  --with-snfit-includedir=DIR directory where the headers were installed],
	[SNFIT_CFLAGS=-I$withval],
    )

    AC_ARG_WITH(snfit-libdir,
	[  --with-snfit-libdir=DIR directory where the library was installed],
	[SNFIT_LIBS=-L$withval $SNFIT_LIBS],
    )

    CPPFLAGS_sav="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $SNFIT_CFLAGS"

    AC_CHECK_HEADER(salt2model.h, ,
	[AC_MSG_WARN([  *** Header file salt2model.h not found.])])

    LIBS_sav="$LIBS"
    LIBS="$LIBS $SNFIT_LIBS"

    AC_CHECK_LIB(snfit, LcParam, ax_snfit_ok=yes,
    	[AC_MSG_WARN([ *** Library libsnfit not linked properly])])

    # FIXME: check the version
    # version_major=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
    # version_minor=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`
    # AC_MSG_CHECKING(that snfit version >=$version_major.$version_minor)

   AC_SUBST(SNFIT_CFLAGS)
   AC_SUBST(SNFIT_LIBS)

   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_snfit_ok" = xyes; then
        ifelse([$2],,AC_DEFINE(HAVE_SNFIT,1,[Define if you have SNFIT library.]),[$2])
        :
else
        ax_snfit_ok=no
        $3
fi
])
