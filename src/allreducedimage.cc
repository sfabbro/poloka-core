#include <iostream>

#include "allreducedimage.h"

/* this file handles reloading previously written ReducedImage's
and inheriters (TransformedImage and ImageSum).
The basic thing that happens is to read the actual type
construct the corresponding object, and return it
as a pointer on the base type (ReducedImage) here.
OO databases (at least Root) do that automagically. 
If we use root for all I/O's, this file just disappears. */

#include "reducedimage.h"
#include "imagesum.h"
#include "transformedimage.h"
#include "imagesubtraction.h"
#include "swarpstack.h"
#include "subimage.h"

ReducedImage* ReducedImageNew(const string &Name)
{
ReducedImage redImage(Name);
if (!redImage.IsValid()) return NULL;
 string typeName;
 if (FileExists(redImage.Dir()+"type.name"))
   typeName = redImage.StoredTypeName();
 else
   typeName = redImage.TypeName();

if (typeName == "ReducedImage")
  {
    return new ReducedImage(Name);
  }
else if (typeName == "TransformedImage")
  {
    return new TransformedImage(Name);
  }
else if (typeName == "ImageSum")
  {
    return new ImageSum(Name);
  }
else if (typeName == "ImageSubtraction")
  {
    return new ReducedImage(Name);// cause ImageSubtraction read is buggy 
  }
else if (typeName == "SubImage")
  {
    return new SubImage(Name);
  }
else if (typeName == "SwarpStack")
  {
    return new SwarpStack(Name);
  }
else if (redImage.HasCatalog()&& redImage.HasImage() && redImage.HasWeight()) // reducedimage that does not know it...
  {
    return new ReducedImage(Name);
  }
else
 cerr << " trying to reload " << Name << " as a ReducedImage, which it does not seem to be " << typeName << endl;
return NULL;
}

const TransformedImage *IsTransformedImage(const ReducedImage *RImage)
{
  return dynamic_cast<const TransformedImage*>(RImage);
}

const ImageSum *IsImageSum(const ReducedImage *RImage)
{
  return dynamic_cast<const ImageSum*>(RImage);
}

TransformedImage *IsTransformedImage(ReducedImage *RImage)
{
  return dynamic_cast<TransformedImage*>(RImage);
}

ImageSum *IsImageSum(ReducedImage *RImage)
{
  return dynamic_cast<ImageSum*>(RImage);
}

ImageSubtraction *IsImageSubtraction(ReducedImage *RImage)
{
  return dynamic_cast<ImageSubtraction*>(RImage);
}

const ImageSubtraction *IsImageSubtraction(const ReducedImage *RImage)
{
  return dynamic_cast<const ImageSubtraction*>(RImage);
}

SubImage *IsSubImage(ReducedImage *RImage)
{
  return dynamic_cast<SubImage*>(RImage);
}

const SubImage *IsSubImage(const ReducedImage *RImage)
{
  return dynamic_cast<const SubImage*>(RImage);
}

SwarpStack *IsSwarpStack(ReducedImage *RImage)
{
  return dynamic_cast<SwarpStack*>(RImage);
}

const SwarpStack *IsSwarpStack(const ReducedImage *RImage)
{
  return dynamic_cast<const SwarpStack*>(RImage);
}
