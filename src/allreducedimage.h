// This may look like C code, but it is really -*- C++ -*-
#ifndef ALLREDUCEDIMAGE__H
#define ALLREDUCEDIMAGE__H

#include <string>
using namespace std;

class ReducedImage;
class ImageSum;
class TransformedImage;

/*! \file 
    \brief utilities around all derived classes ot ReducedImage.
*/

//! 'virtual' constructor. Not totally functionnal unfortunately.
ReducedImage* ReducedImageNew(const string &Name);

//! returns NULL if the argument is not a TransformedImage.
const TransformedImage *IsTransformedImage(const ReducedImage *RImage);
TransformedImage *IsTransformedImage(ReducedImage *RImage);

//! returns NULL if the argument is not an ImageSum.
const ImageSum *IsImageSum(const ReducedImage *RImage);
ImageSum *IsImageSum(ReducedImage *RImage);

#endif /*  ALLREDUCEDIMAGE__H */
