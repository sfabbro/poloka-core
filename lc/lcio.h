// This may look like C code, but it is really -*- C++ -*-
#ifndef LCUTILS__H
#define LCUTILS__H

#include <reducedimage.h>
#include "refstar.h"

//! read a list of objects and a ReducedImageList from a light curve file 
istream& lc_read(istream& Stream, RefStarList& Objects, ReducedImageList& Images);

//! create a light curve file from a list of objects and a ReducedImageList
ostream& lc_write(ostream& Stream, const RefStarList& Objects, const ReducedImageList& Images);

//!
//! \page lightfile Syntax of the "lightfile"
//! 
//! A light curve file consists of objects and images.
//! Comment are coded with the '#' at the beginning of the line.
//! The expected format is the following:
//! 
//! \code
//! 
//! # X Y [DATE_MIN=DD/MM/YYYY DATE_MAX=DD/MM/YYYY NAME=toto ]
//! OBJECTS
//! 510.20 1043.98 DATE_MIN=10/08/2003 DATE_MAX=12/12/2003 NAME=R4D14-4
//! 
//! # The list of all DbImages in a single filter
//! IMAGES
//! DbImage1_R
//! DbImage2_R
//! DbImage1_I
//! DbImage2_I
//!   .
//!   .
//!   .
//! # Optional photometric reference. The default is the best seeing image.
//! PHOREF
//! PhoRefDbImage
//! 
//! \endcode
//! 


#endif // LCUTILS__H