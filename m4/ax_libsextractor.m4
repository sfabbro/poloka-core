# -*- autoconf -*-
# 
# use as AX_LIBSEXTRACTOR([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#

AC_DEFUN([AX_LIBSEXTRACTOR],[

ax_sextractor_ok=no

PKG_CHECK_MODULES([SEXTRACTOR], [sextractor],[ax_sextractor_ok=yes], [ax_sextractor_ok=no])

if test x"$ax_sextractor_ok" = x"no"; then
   # if neither pkg-config file found, nor SEXTRACTOR_LIBS defined
   if test x$SEXTRACTOR_LIBS = x ; then
      SEXTRACTOR_LIBS="-lsextractor -lm"
      AC_ARG_WITH([sextractor-libdir],
	 [  --with-sextractor-libdir=DIR directory where the library was installed],
	 [SEXTRACTOR_LIBS="-L$withval $SEXTRACTOR_LIBS"], )
   fi
   ax_sextractor_lib_ok=no
   LIBS_sav="$LIBS"
   LIBS="$LIBS $SEXTRACTOR_LIBS"
   AC_CHECK_LIB([sextractor], [makeit], [ax_sextractor_lib_ok=yes], [AC_MSG_WARN([ *** SExtractor library not found])])

   # if neither pkg-config file found, nor SEXTRACTOR_CFLAGS defined
   if test x$SEXTRACTOR_CFLAGS = x ; then
      SEXTRACTOR_CFLAGS=""
      AC_ARG_WITH(sextractor-includedir,
	 [  --with-sextractor-includedir=DIR directory where the headers were installed],
	 [SEXTRACTOR_CFLAGS="-I$withval"], )
   fi
   ax_sextractor_hdr_ok=no
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $SEXTRACTOR_CFLAGS"
   AC_CHECK_HEADER([define.h],
      [ax_sextractor_hdr_ok=yes],
      [AC_MSG_WARN([  *** SExtractor library headers not found.])])

   if test x$ax_sextractor_lib_ok = xyes -a x$ax_sextractor_hdr_ok = xyes; then
      ax_sextractor_ok=yes
   fi
   AC_SUBST(SEXTRACTOR_CFLAGS)
   AC_SUBST(SEXTRACTOR_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_sextractor_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_SEXTRACTOR, [1], [Define if you have SExtractor library.]),[$1])
        :
else
        ax_sextractor_ok=no
        $2
fi
])
