#!/bin/sh 
# 
# $Id: autogen.sh,v 1.1 2004/02/23 10:49:47 nrl Exp $
# 
# autogen.sh 
# 
# run the autotools in the correct order.
#

autoheader
aclocal -I m4
libtoolize --copy --force
autoconf
automake -a --copy
