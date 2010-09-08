#!/bin/sh 
# 
# $Id: autogen.sh,v 1.4 2010/09/08 22:09:11 seb Exp $
# 
# autogen.sh 
# 
# run the autotools in the correct order.
#

test -d config || mkdir config
autoheader
aclocal -I m4
libtoolize --copy --force
autoconf
automake -a --copy
