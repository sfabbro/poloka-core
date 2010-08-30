# -*- autoconf -*-
# 
# use as AX_IMAGEMAGICK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# no version checking

AC_DEFUN([AX_IMAGEMAGICK],[

ax_imagemagick_ok=no

PKG_CHECK_MODULES([IMAGEMAGICK], [ImageMagick++], ax_imagemagick_ok=yes,[
    if test x$IMAGEMAGICK_LIBS = x ; then
       # try to find Magick++-config in PATH
       AC_PATH_PROG(PROG, Magick++-config, no)
       if test "$PROG" = "no" ; then
           AC_MSG_ERROR([cannot find Magick++-config either])
       fi
       AC_MSG_RESULT(yes)
       # get includes and use_linkopts from Magick++
       IMAGEMAGICK_CFLAGS="`Magick++-config --cppflags` `Magick++-config --cflags`"
       IMAGEMAGICK_LIBS="`Magick++-config --libs`"

       # save these values for usage in Makefile.am
       AC_SUBST(IMAGEMAGICK_CFLAGS)
       AC_SUBST(IMAGEMAGICK_LIBS)
    fi
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_imagemagick_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_IMAGEMAGICK,1,[Define if you have IMAGEMAGICK library.]),[$1])
        :
else
        ax_imagemagick_ok=no
        $2
fi
])
