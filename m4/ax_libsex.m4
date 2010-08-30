# -*- autoconf -*-
# 
# use as AX_LIBSEX([MIN VERSION[, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]]])
# ex: AX_LIBSEX(2.2)
#

AC_DEFUN([AX_LIBSEX],[

ax_sex_ok=no

version_req=[$1]
version_check=
if test x$version_req != x; then
   version_check=">= $version_req"
fi

PKG_CHECK_MODULES([SEX], [sex $version_check], ax_sex_ok=yes,[
    if test x$SEX_LIBS = x ; then
       SEX_LIBS="-lsex -lm"
    fi

    AC_ARG_WITH(sex,
	[  --with-sex=DIR        directory where SExtractor source was compiled],
    	[SEX_CFLAGS="-I$withval"
	 SEX_LIBS="-L$withval $SEX_LIBS"],
    )

    AC_ARG_WITH(sex-includedir,
	[  --with-sex-includedir=DIR directory where the headers were installed],
	[SEX_CFLAGS=-I$withval],
    )

    AC_ARG_WITH(sex-libdir,
	[  --with-sex-libdir=DIR directory where the library was installed],
	[SEX_LIBS=-L$withval $SEX_LIBS],
    )

    CPPFLAGS_sav="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $SEX_CFLAGS"

    AC_CHECK_HEADER(define.h, ,
	[AC_MSG_WARN([  *** Header file define.h not found.])])

    LIBS_sav="$LIBS"
    LIBS="$LIBS $SEX_LIBS"

    AC_CHECK_LIB(sex, makeit, ax_sex_ok=yes,
    	[AC_MSG_WARN([ *** Library sex not linked properly])])

    # FIXME: check the version
    # version_major=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
    # version_minor=`echo $version_req | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`
    # AC_MSG_CHECKING(that sex version >=$version_major.$version_minor)
	
   AC_SUBST(SEX_CFLAGS)
   AC_SUBST(SEX_LIBS)

   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_sex_ok" = xyes; then
        ifelse([$2],,AC_DEFINE(HAVE_SEX,1,[Define if you have SExtractor library.]),[$2])
        :
else
        ax_sex_ok=no
        $3
fi
])
