#!/bin/sh 
# 
# $Id: autogen.sh,v 1.2 2006/12/22 13:35:40 guy Exp $
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
