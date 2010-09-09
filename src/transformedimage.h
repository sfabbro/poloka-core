// This may look like C code, but it is really -*- C++ -*-
#ifndef TRANSFORMEDIMAGE__H
#define TRANSFORMEDIMAGE__H 

#include "gtransfo.h"
#include "sestar.h"
#include "apersestar.h"
#include "reducedimage.h"
#include "frame.h"
#include "countedref.h"

class Image;
class FitsImage;

/* a (virtual) class for generic image transformation */

class ImageTransfo : public RefCount
{
 protected:
    string SourceName();
  public :

//!
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

    
  //!  
  virtual ImageTransfo* Clone() const = 0;

  //!
  virtual void dump(ostream &s = cout) const = 0;

  virtual ~ImageTransfo() {} ;

};

typedef CountedRef<ImageTransfo> ImageTransfoRef;


/*! geometric transfo of a reduced image.  */

class ImageGtransfo : public ImageTransfo {
 private:
  string geomRefName; 
  GtransfoRef transfoFromRef;
  GtransfoRef transfoToRef;
  double scaleFactor; // 1d scale factor (== sqrt(jacobian))
  Frame outputImageSize; // the frame on which we want the input image to be resampled.
  
public:

  //!
  ImageGtransfo(const Gtransfo* TransfoFromRef, const Gtransfo* TransfoToRef, const Frame &OutputImageSize, const string &GeomRefName );
  //! the output image size is the one of the Ref. Finds the transfo(s).
  ImageGtransfo(const ReducedImage &Ref, const ReducedImage& ToAlign,float min_match_ratio=0);
  //!
  ImageGtransfo();
  //!
  const GtransfoRef TransfoFromRef() const {return transfoFromRef;}
  //!
  const GtransfoRef TransfoToRef() const {return transfoToRef;}
  //!
  const GtransfoRef FromRef() const;

  //!
  string GeomRefName() const { return geomRefName;}

  //!
  virtual void  dump(ostream &s = cout)const ;

  //! the one that transforms the image and update header.
  void TransformImage(const FitsImage &Source, FitsImage& Transformed, 
		      const ReducedImage *RSource, ReducedImage *Result, 
		      double DefaultVal = 0) const ;

  //! Transforms a weight image.
  void TransformWeightImage(const FitsImage &Source, FitsImage& Transformed) const;
  //! transforms a bool Image
  void TransformBoolImage(const FitsImage &Source, FitsImage& Transformed) const ;
  void TransformCatalog(const SEStarList &Catalog, SEStarList &Transformed) const;
  void TransformAperCatalog(const AperSEStarList &Catalog, AperSEStarList &Transformed) const;



  bool IsValid() const;
  double ScaleFactor() const {return scaleFactor;}
  ImageTransfo* Clone() const;
  // bool Write(const string &Name) const;
  
  ~ImageGtransfo();


};

typedef CountedRef<ImageGtransfo> ImageGtransfoRef;


//! class that operates the transformation of a ReducedImage (image(s) + list) 
/*! As for other descendants of ReducedImage, the actual computations occur
when you request the name of a data file (e.g. FitsName()).
To geometrically align a set of images on the same reference, use
  ImagesAlign(). If you want to sum them, uses ImagesAlignAndSum() */
class TransformedImage : public ReducedImage {
private:
  ImageTransfoRef transfo;
  string sourceName;
  ReducedImageRef source;  //!
  //ReducedImageRef geomRef; //!
  void init(const ReducedImage &Source,const ImageTransfo *Transfo);


  public :
  //! to create a new TransformedImage, or locate an existing one.
  /*! If you want to align A on B, the typical constructor call will be:
    TransformedImage(NewName, A, &ImageGtransfo(B,A));
    To align a set of images on the same reference, use ImagesAlign().
    
  */
  TransformedImage(const string &Name, const ReducedImage &Source, const ImageTransfo *Transfo);

  //! assumes that the TransformedImage already exists
    TransformedImage(const string &Name);

  //! empty constructor for persistence
  TransformedImage(){};

  virtual const string  TypeName() const { return "TransformedImage";}

  //! Original (untransformed) image name
  string  SourceName() const { return sourceName;}

  //! involved transformation.
  const ImageTransfoRef Transfo() const {return transfo;};

  //! Gtransfo from reference to transformed image. Assumes that ImageTransfo is an ImageGtransfo.
  const Gtransfo* FromRef() const;

//! ImageGtransfo from reference to transformed image. Assumes that ImageTransfo is an ImageGtransfo.
  const  ImageGtransfoRef IMAGEGTransfo() const ;


  virtual void dump(ostream& s = cout) const;

  virtual bool MakeFits() ;
  virtual bool MakeCatalog();
  bool MakeAperCat();
  virtual bool MakeDead();
  virtual bool MakeSatur();
  virtual bool MakeCosmic();
  virtual bool MakeSatellite();
  virtual bool MakeWeight();
  ReducedImage* Clone() const;
  TransformedImage(const TransformedImage &Original);
  
  // void SetTranfo(ImageTransfo *transfo);
  // TransformedImage(string &ImageTransfoName, ReducedImage &Source, string Name); 
  ~TransformedImage();

};

#include "imagelist.h"
//! a handy typedef
//#ifdef STORAGE : en storage jusque la parce que
// les routines genre ImageSum manipulent des ReducedImageList
// et non des TransformedImageList.
typedef ImageList<TransformedImage> TransformedImageList;
typedef TransformedImageList::iterator       TransformedImageIterator;
typedef TransformedImageList::const_iterator TransformedImageCIterator;
//#endif

//! To standardize transformed image names
string TransformedName(const string &ToTransform, const string &Ref);    

// ! Geometrically align images on a given reference. The image is done by default, and eventually more according to ToDo.
/*! ToDo may be constructed using the tags DoFits DoCatalog DoDead DoSatur, e.g. provide
   DoCatalog|DoDead to get catalog and dead frame on top of the default image itself. */

int ImagesAlign(const ReducedImageList &ToAlign, const ReducedImage &Reference, ReducedImageList &Aligned, const int ToDo,bool use_wcs=false,float min_match_ratio=0);

void MakeUnionRef(const ReducedImageList& ToAlign, const ReducedImage& Reference, const string& unionName);

#endif /*  TRANSFORMEDIMAGE__H */
