// This may look like C code, but it is really -*- C++ -*-
#ifndef ALLREDUCEDIMAGE__H
#define ALLREDUCEDIMAGE__H

#include <string>
using namespace std;

class ReducedImage;
class ImageSum;
class TransformedImage;
class ImageSubtraction;
class SubImage;
class SwarpStack;

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

//! returns NULL if the argument is not an ImageSubtraction.
const ImageSubtraction *IsImageSubtraction(const ReducedImage *RImage);
ImageSubtraction *IsImageSubtraction(ReducedImage *RImage);

//! returns NULL if the argument is not a SubImage.
const SubImage *IsSubImage(const ReducedImage *RImage);
SubImage *IsSubImage(ReducedImage *RImage);

//! returns NULL if the argument is not a SwarpStack.
const SwarpStack *IsSwarpStack(const ReducedImage *RImage);
SwarpStack *IsSwarpStack(ReducedImage *RImage);

#endif /*  ALLREDUCEDIMAGE__H */
