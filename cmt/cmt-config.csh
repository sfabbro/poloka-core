#!/bin/csh -f

if ( $#argv < 1 || $1 == "--help" ) then
    echo "cmt-config.csh <option>"
    echo "               --ldflags   "
    echo "               --cflags "
    echo "               --help"
    exit
endif

# get includes necessary to compile this package with the help of cmt
# get ldflags

# copy requirements.in in order to configure this package
echo "package poloka" >! requirements
cat requirements.in | grep -v '@' >> requirements
if ( ! -f version.cmt ) echo "v0r0" >  version.cmt

set packages = `cmt show uses | grep -v "#" | grep -v "CMT" | awk '{printf("%s ",$2);}' `


if ( $1 == "--ldflags" ) then
    set LDFLAGS = ""

    foreach package ( $packages ) 
	set package_LDFLAGS = `cmt show macro_value ${package}_linkopts`
	set LDFLAGS = "$LDFLAGS $package_LDFLAGS"
    end 
    echo $LDFLAGS
endif

if ( $1 == "--cflags" ) then
    set CFLAGS = ""

    foreach package ( $packages ) 
	set package_CFLAGS = `cmt show macro_value ${package}_cflags`
	set CFLAGS = "$CFLAGS $package_CFLAGS"
    end 
    echo $CFLAGS
endif
