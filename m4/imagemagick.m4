# -*- autoconf -*-
# 
# 
# autoconf macro to check imagemagick
# 

####################################### 
AC_DEFUN([CHECK_IMAGEMAGICK],[

# message
AC_MSG_CHECKING(ImageMagick configuration)

# try to find Magick++-config in PATH
AC_PATH_PROG(PROG, Magick++-config, no)
if test "$PROG" = "no" ; then
    AC_MSG_ERROR(cannot find Magick++-config)
fi

AC_MSG_RESULT(yes)

# get includes and use_linkopts from Magick++
IMAGEMAGICK_CPPFLAGS="`Magick++-config --cppflags`"
IMAGEMAGICK_LIBS="`Magick++-config --libs | sed 's/-lxml2//'`"
IMAGEMAGICK_LDFLAGS="`Magick++-config --ldflags`"

# save these values for usage in Makefile.am
AC_SUBST(IMAGEMAGICK_CPPFLAGS)
AC_SUBST(IMAGEMAGICK_LDFLAGS)
AC_SUBST(IMAGEMAGICK_LIBS)


]
)
