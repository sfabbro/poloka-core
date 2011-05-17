# -*- autoconf -*-
# 
# use as AX_LIBSEX([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#

AC_DEFUN([AX_LIBSEX],[

ax_sex_ok=no

PKG_CHECK_MODULES([SEX], [sex],[ax_sex_ok=yes], [ax_sex_ok=no])

if test x"$ax_sex_ok" = x"no"; then
   # if neither pkg-config file found, nor SEX_LIBS defined
   if test x$SEX_LIBS = x ; then
      SEX_LIBS="-lsex -lm"
      AC_ARG_WITH([sex-libdir],
	 [  --with-sex-libdir=DIR directory where the library was installed],
	 [SEX_LIBS="-L$withval $SEX_LIBS"], )
   fi
   ax_sex_lib_ok=no
   LIBS_sav="$LIBS"
   LIBS="$LIBS $SEX_LIBS"
   AC_CHECK_LIB([sex], [makeit], [ax_sex_lib_ok=yes], [AC_MSG_WARN([ *** SExtractor library not found])])

   # if neither pkg-config file found, nor SEX_CFLAGS defined
   if test x$SEX_CFLAGS = x ; then
      SEX_CFLAGS=""
      AC_ARG_WITH(sex-includedir,
	 [  --with-sex-includedir=DIR directory where the headers were installed],
	 [SEX_CFLAGS="-I$withval"], )
   fi
   ax_sex_hdr_ok=no
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $SEX_CFLAGS"
   AC_CHECK_HEADER([sex/define.h],
      [ax_sex_hdr_ok=yes],
      [AC_MSG_WARN([  *** SExtractor library headers not found.])])

   if test x$ax_sex_lib_ok = xyes -a x$ax_sex_hdr_ok = xyes; then
      ax_sex_ok=yes
   fi
   AC_SUBST(SEX_CFLAGS)
   AC_SUBST(SEX_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_sex_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_SEX, [1], [Define if you have SExtractor library.]),[$1])
        :
else
        ax_sex_ok=no
        $2
fi
])
