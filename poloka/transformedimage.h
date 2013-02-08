// This may look like C code, but it is really -*- C++ -*-
#ifndef TRANSFORMEDIMAGE__H
#define TRANSFORMEDIMAGE__H 

#include <poloka/gtransfo.h>
#include <poloka/sestar.h>
#include <poloka/apersestar.h>
#include <poloka/reducedimage.h>
#include <poloka/frame.h>
#include <poloka/countedref.h>

class Image;
class FitsImage;

//! a (virtual) class for generic image transformation
class ImageTransfo : public RefCount
{
  
public :

  //! checks validity
  virtual bool IsValid() const = 0;
  
  //! applies the transfo to the image.
  virtual void TransformImage(const FitsImage &Source, FitsImage& Transformed,
			      const ReducedImage *RSource, 
			      ReducedImage *Result, double DefaultVal = 0) const = 0;

  //! applies the transfo to the image. accounts for the fact that the image actually conatins weights
  virtual void TransformWeightImage(const FitsImage &Source, 
				    FitsImage& Transformed) const = 0;
  
  //! applies the transfo to the image.
  virtual void TransformBoolImage(const FitsImage &Source, 
				  FitsImage& Transformed) const = 0;

  //! applies the transfo to the list
  virtual void TransformCatalog(const SEStarList &Catalog, 
				SEStarList &Transformed) const = 0;
  
  //! applies the transfo to the list
  virtual void TransformAperCatalog(const AperSEStarList &Catalog, 
				    AperSEStarList &Transformed) const = 0;
    
  //! makes a copy (probably useless with counted ref)
  virtual ImageTransfo* Clone() const = 0;

  //! standard printing
  virtual void dump(ostream &s = cout) const = 0;

  virtual ~ImageTransfo() {}
};

typedef CountedRef<ImageTransfo> ImageTransfoRef;


//! geometric transfo of a reduced image
class ImageGtransfo : public ImageTransfo {

private:

  string geomRefName;
  GtransfoRef transfoFromRef;
  GtransfoRef transfoToRef;
  Frame outputImageSize; // the frame on which we want the input image to be resampled.
  
public:

  //!
  ImageGtransfo(const GtransfoRef TransfoFromRef, const GtransfoRef TransfoToRef, const Frame &OutputImageSize, const string &GeomRefName);

  //! the output image size is the one of the Ref. Finds the transfo(s).
  ImageGtransfo(const ReducedImage &Ref, const ReducedImage& ToAlign);

  //!
  ImageGtransfo() {}

  //!
  ~ImageGtransfo() {}

  //!
  const GtransfoRef TransfoFromRef() const { return transfoFromRef; }

  //!
  const GtransfoRef TransfoToRef() const { return transfoToRef; }

  //!
  const GtransfoRef FromRef() const;

  //! geometric DbImage reference
  string GeomRefName() const { return geomRefName; }

  //! print basic info
  virtual void  dump(ostream &s = cout) const;

  //! the one that transforms the image and update header.
  void TransformImage(const FitsImage &Source, FitsImage& Transformed, 
		      const ReducedImage *RSource, ReducedImage *Result, 
		      double DefaultVal = 0) const ;

  //! transforms a weight image.
  void TransformWeightImage(const FitsImage &Source, FitsImage& Transformed) const;

  //! transforms a bool Image
  void TransformBoolImage(const FitsImage &Source, FitsImage& Transformed) const ;

  //! transforms sextractor catalog
  void TransformCatalog(const SEStarList &Catalog, SEStarList &Transformed) const;

  //! transforms aperture catalog
  void TransformAperCatalog(const AperSEStarList &Catalog, AperSEStarList &Transformed) const;

  bool IsValid() const;
  double ScaleFactor() const;
  ImageTransfo* Clone() const;

};

typedef CountedRef<ImageGtransfo> ImageGtransfoRef;


//! class that operates the transformation of a ReducedImage (image(s) + list) 
//! As for other descendants of ReducedImage, the actual computations occur
//!  when you request the name of a data file (e.g. FitsName()).
class TransformedImage : public ReducedImage {

private:

  ImageTransfoRef transfo;
  ReducedImageRef source;
  void init(const ReducedImage &Source,const ImageTransfo *Transfo);

public :

  //! to create a new TransformedImage, or locate an existing one.
  /*! If you want to align A on B, the typical constructor call will be:
    TransformedImage(NewName, A, &ImageGtransfo(B,A));
    To align a set of images on the same reference, use ImagesAlign(). */
  TransformedImage(const string &Name, const ReducedImage &Source, const ImageTransfo *Transfo);

  //! assumes that the TransformedImage already exists
  TransformedImage(const string &Name);

  //! empty constructor for persistence
  TransformedImage() {}

  TransformedImage(const TransformedImage &Original);

  ~TransformedImage() {}

  ReducedImageRef Clone() const;

  virtual const string  TypeName() const { return "TransformedImage"; }

  //! Original (untransformed) image name
  string  SourceName() const { return source->Name(); }

  //! involved transformation.
  const ImageTransfoRef Transfo() const { return transfo; }

  //! Gtransfo from reference to transformed image. Assumes that ImageTransfo is an ImageGtransfo.
  const GtransfoRef FromRef() const;

  //! ImageGtransfo from reference to transformed image. Assumes that ImageTransfo is an ImageGtransfo.
  const  ImageGtransfoRef IMAGEGTransfo() const;

  //! standard dump
  virtual void dump(ostream& s = cout) const;

  virtual bool MakeFits() ;
  virtual bool MakeCatalog();
  virtual bool MakeAperCat();
  virtual bool MakeStarCat();
  virtual bool MakeDead();
  virtual bool MakeSatur();
  virtual bool MakeCosmic();
  virtual bool MakeSatellite();
  virtual bool MakeWeight();

  string TransfoFileName() const { return Dir()+"/imagetransfo.dat"; }
};

//! To standardize transformed image names
string TransformedName(const string &ToTransform, const string &Ref);    

// ! Geometrically align images on a given reference. The image is done by default, and eventually more according to ToDo.
/*! ToDo may be constructed using the tags DoFits DoCatalog DoDead DoSatur, e.g. provide
  DoCatalog|DoDead to get catalog and dead frame on top of the default image itself. */
int ImagesAlign(const ReducedImageList &ToAlign, const ReducedImage &Reference, ReducedImageList &Aligned, const int ToDo, bool WcsOnly=false);

#endif /*  TRANSFORMEDIMAGE__H */
