// This may look like C code, but it is really -*- C++ -*-
#ifndef REDUCEDIMAGE__H
#define REDUCEDIMAGE__H


#include <iostream>
#include <string>

using namespace std;

#define DoFits 1
#define DoCatalog 2
#define DoDead 4
#define DoSatur 8
#define DoWeight 16
#define DoCosmic 32
#define DoSatellite 64
#define DoAperCatalog 128

#define PutVar 1
#define PutCosmic 2
#define PutFlat 4
#define PutDead 8
#define PutSattelite 16


#include "dbimage.h" // for inheritance.
#include "frame.h"
#include "fileutils.h" // for FileExists
#include "countedref.h"

class Gtransfo;

class ForSExtractor;

/*! \file reducedimage.cc reducedimage.h */
/*!
    \brief a handle to access data associated to an image: the fits file, the catalog,
    the dead and satur frames, and a set of 'scalars' such as seeing, saturation level &co.

   This class is meant to be derived. One should provide a correct Clone routine in 
   derived classes. It is also important for ReducedImageList code 
   performance that derivations remain small and fast to copy, in particular, 
   they should not contain actual Image's.
   
*/

//class Gtransfo;
class ReducedImage : public DbImage
{


private :
  // take care if you add any pointer here : 
  // the default copy constructor then has to be written
  bool actuallyReduced;//!

public :
  //!
  ReducedImage(): actuallyReduced(false){};
  //!
  explicit ReducedImage(const DbImage&);
  //!
  ReducedImage(const string &Name);


  virtual const string  TypeName() const { return "ReducedImage";}

  //!
  bool HasImage() const { return (FileExists(FitsName()));}
  bool HasBack() const { return (FileExists(FitsBackName()));}
  bool HasMiniBack() const { return (FileExists(FitsMiniBackName()));}
  bool HasCatalog() const {return(FileExists(CatalogName()));}
  bool HasDead() const {return(FileExists(FitsDeadName()));}
  bool HasFlat() const {return(FileExists(FitsFlatName()));}
  bool HasSatur() const {return(FileExists(FitsSaturName()));}
  bool HasCosmic() const {return(FileExists(FitsCosmicName()));}
  bool HasSatellite() const {return(FileExists(FitsSatelliteName()));}
  bool HasWeight() const {return(FileExists(FitsWeightName()));}
  bool HasBad() const {return(FileExists(FitsBadName()));}

//! not const because it may actually compute the image and other things (for derived class)
  string FitsName()  const { return FitsImageName(Calibrated);}
 //! produce fits image
  virtual bool MakeFits();

  //!
  string CatalogName()  const { return ImageCatalogName(SExtractor);}

  /* Pour remplir le carnet d'ordre de SExtractor */
  void  FillSExtractorData(ForSExtractor & data, 
			   bool fond_deja_soustrait, bool sauver_fond,
			   bool use_sigma_header);

  void  FillSExtractorData_2(ReducedImage & rim_det, ForSExtractor & data);

  /* Recover Back eventually from mini back image and 
     re add it if desired */
bool
RecoverBack(bool add_to_im);
  /* As indicated: enables to start from an image with ist background re-added when possible */
bool
ReAddBackground_and_ResetKeys();


/* Runs SExtractor on the image and compute the seeing from
SExtractor Catalog. See belows for details according 
to the options.
argument (argument value in  MakeCatalog()): 
overwrite(false): 
  - if = false: if the catalog already exists, then exit 
immediately. 
  - if = true: any existing catalog is overwritten.

savemasksat (true): 
  - if = true: the SExtractor computed  saturated pixels 
map is saved.

pas_sub_fond(false): 
  - if = true: the SExtractor computed background map 
won't be subtracted from the image. 
  - if = false: the SExtractor computed background 
will be subtracted, providing this hasn't already been 
done (this is tracked with BackSub())
use_sigma_header (false): 
  - if = true:  SExtractor will use the SigmaBackground() 
to set his detection and photometric  levels. 
Usefull in case of artificially smoothed images 
(geometrically transformed images, convolved images).*/

  bool MakeCatalog(bool redo_from_beg, bool overwrite, bool savemasksat,
		   bool pas_sub_fond, 
		   bool use_sigma_header);


  // pour processer 2 images : detection + mesure, en simplifie
  //weight_from_measurement_image : detection et measurement weight map = measurement weight map.
  // catalog_name : Dir()+"/sedble.list" par exemple.
  bool MakeCatalog_2images(ReducedImage & rim_det, bool overwrite, 
			   bool weight_from_measurement_image, 
			   string catalog_name, bool do_segmentation);


//! Produce the Saturated stars pixels mask, subtract the image background, detect with the SExtractor computed sigma. search the cosmics, and update catalog and weight for cosmics. No free coffee.
  virtual bool MakeCatalog();


  //! produces aperture photometry on the objects from the SE catalog
  virtual bool MakeAperCat();


  //! Extracts stars form the AperCat
  virtual bool MakeStarCat();
  
//! MakeCatalog_ImageBizarre() is for sum-images, or convolved images, for which we do not want to subtract the background map, nor compute the saturated pixels map, and for which we provide the value of the sigmabackground. overwrite is set to true.
  bool MakeCatalog_ImageBizarre(){ 
    return MakeCatalog(/*redo_from_beg=*/false, /*overwrite = */ true,
		       /*savemasksat = */ false, 
		       /*pas_sub_fond= */true,
		       /*use_sigma_header = */ true);};
				     

  //! produce satur image
  virtual bool MakeSatur();
  
  //! produce dead image
  virtual bool MakeDead();
  
  //!
  void FlagCosmicsInCatalog(const Image &CosmicImage, 
			    const double dist=2);
  //! produce cosmic image
  virtual bool MakeCosmic();


  //! produce satellite image
  virtual bool MakeSatellite();

  //! produce weight image
  bool IsToadsWeight();
  virtual bool MakeWeight();

  //! produce bad pixel image
  virtual bool MakeBad();



  //! returns if both fits image and catalog file exist.
  bool ActuallyReduced() const;

  //! shorthand call for Make{Fits,Catalog,Dead,Satur}. ToDo may conveniently be contructed using predefined tags DoFits DoCatalog DoDead DoSatur.
  bool Execute(const int ToDo);

  //! dumps basic info.
  virtual void dump(ostream &s = cout ) const;
  friend ostream& operator << (ostream &stream, const ReducedImage &red)
    {red.dump(stream); return stream; }
    
  //! size in x-axis
  int XSize() const;
  //! size in y-axis
  int YSize() const;

 

 //! basic seeing
  double Seeing() const;
  void RemoveSeeing();
  bool SetSeeing(const double &Value, const string Comment = "");

 //! bseeing from gaussian fits to the objects+ star clump finding
  double GFSeeing() const;
  void RemoveGFSeeing();
  bool SetGFSeeing(const double &Value, const string Comment = "");



  //! basic PSF shape parameters
  bool GetPsfShapeParams(double &SigmaX, double &SigmaY, double &ThetaXY ) const;
  bool SetPsfShapeParams(const double &SigmaX, const double &SigmaY, 
			 const double &ThetaXY, const string Comment = "");

  //! basic coordinates
  bool GetRaDecEpoch(double &Ra, double &Dec, double &Epoch) const;
  bool SetRaDecEpoch(const double &Ra, const double &Dec, 
			 const double &Epoch, const string Comment = "");

  //! the (average) sky level as it should appear on the image.
  void RemoveBackLevel();

  //! the current back level
  double BackLevel() const;

  //! the back level before it was subtracted
  double BackLevelNoSub() const;

  double SESky() const;
  bool SetBackLevel(const double &Value, const string Comment= "");
  bool SetSESky(const double &Value,  const string Comment = "");
  void RemoveSESky();

  //! the (average) sky level at it would appear if not subtracted
  double OriginalSkyLevel() const;
  bool SetOriginalSkyLevel(const double &Value, const string Comment = "");
 
  //! check if the sky background has been subtracted
  bool IsSkySub() const;

  //! r.m.s of background
  double SigmaBack() const;

  bool SetSigmaBack(const double &Value, const string Comment = "");


  bool SetSESigma(const double &Value, const string Comment="");
  void RemoveSESigma();
  //! actual noise in the image
  double NoisePow() const;
  bool SetNoisePow(const double &Value, const string Comment = "");

  //! wether background was subtracted or not
  bool BackSub() const;
  bool SetBackSub(const bool &Value, const string Comment="");

  //! current saturation level
  double Saturation() const;
  bool SetSaturation(const double &Value, const string Comment = "");

  //! orginal saturation level
  double OriginalSaturation() const;
  bool SetOriginalSaturation(const double &Value, const string Comment = "");

  //! exposure time
  double Exposure() const;
  bool SetExposure(const double &Value, const string Comment = "");
  
 //! A lot of zero points .....................
  
  //! zero point as measured with USNO Catalog
  double ZeroPoint() const;
  bool SetZeroPoint(const double &Value, const string Comment = "");
  
//! Zero Point as computed from instrument specifications
  double Zerop() const;
  bool SetZerop(const double &Value, const string Comment = "");
 
 //! zero point from ZP0 key, a Zero Point that was once thought to be good ......
  double ZP0() const;
  bool SetZP0(const double &Value, const string Comment = "");
  bool HasZP0() const;

//! zero point from ZP key, supposed to be better than ZP0
  double ZP() const;
  bool SetZP(const double &Value, const string Comment = "");
  bool HasZP() const;

//! Zero Point set by toads : read for 1 image (photometric ref)  with AnyZeroPoint routine and then propagated according to the relationship to this image. 
//! The key is ZPTOADS, NOT to be set by hand or by an exterior job.
  //! thus is present only if it was set priorily in a TOADS program (for example on the sub in sub.cc).
  //! Otherwise use AnyZeroPoint()
  double ZZZeroP() const;
  bool HasZZZeroP() const;
  bool SetZZZeroP(const double &Value, const string Comment = "");
  void RemoveZZZeroP();

//! gives Zero Point by order of prefrence: ZZZeroP, then ZP, then ZP0, then Zerop.
  double AnyZeroPoint() const ;


  
  //! date of observation
  string Date() const;
  bool SetDate(const string &Value, const string Comment = "");

  //! time of observation
  string TimeObs() const;
  bool SetTimeObs(const string &Value, const string Comment = "");

  //! reduced julian date of observation
  double JulianDate() const;
  bool SetJulianDate(const double &Value, const string Comment = "");
  double ModifiedModifiedJulianDate() const;
  double ModifiedJulianDate() const;
  
  //! signal to noise at magnitude 23
  double SignalToNoise23() const;
  bool SetSignalToNoise23(const double &Value, const string Comment = "");
  
  //! instrument used 
  string Instrument() const;

  //! telescope used 
  string Telescope() const;

  //! chip used 
  int Chip() const;
  bool SetChip(const int &Value, const string Comment = "");


  //! target aimed at 
  string Target() const;
  bool SetTarget(const string &Value, const string Comment = "");


  //! filter band 
  string Filter() const;
  bool SetFilter(const string &Value, const string Comment = "");


  //! filter band 
  string Band() const;
  bool SetBand(const string &Value, const string Comment = "");
  
  //! flatfielding noise
  double FlatFieldNoise() const;
  bool SetFlatFieldNoise(const double &Value, const string Comment = "");

  //! PSF error noise
  double ProfileError() const;
  bool SetProfileError(const double &Value, const string Comment = "");

  //! read out noise
  double ReadoutNoise() const;
  bool SetReadoutNoise(const double &Value, const string Comment = "");
 
  //! gain
  double Gain() const;
  bool SetGain(const double &Value, const string Comment = "");
  bool SetOldGain(const double &Value, const string Comment = "");
 
  //! pixel size in arcsec
  double PixelSize() const;
  bool SetPixelSize(const double &Value, const string Comment = "");
  
  //! right ascension in degree (J2000)
  double RaDeg2000() const;
  bool SetRaDeg2000(const double &Value, const string Comment = "");
  
  //! declination in degree (J2000)
  double DecDeg2000() const;
  bool SetDecDeg2000(const double &Value, const string Comment = "");

  //! epoch
  double Epoch() const;
  bool SetEpoch(const double &Value, const string Comment = "");
  
  double Airmass() const;
  bool SetAirmass(const double &Value, const string Comment = "");

  Gtransfo *RaDecToPixels() const;
  Gtransfo *PixelsToRaDec() const;
 
  //! photometric reference (i.e. image that should have the same flux)
  string PhotomReference() const;
  bool SetPhotomReference(const string &Value, const string Comment = "");



  //! usable part defined by a frame keyword in the header 
  Frame UsablePart() const;

  bool SetUsablePart(const Frame &NewFrame);

  //! the actual physical size (pixels) of the image
  Frame PhysicalSize() const;

  bool IsGoodImage() const;
  //! returns the overlapping area in arcmin^2 with another reducedimage
  double OverlapArcmin2(const ReducedImage& Other) const;

  //! fitsheader equivalent
  bool SameChipFilterInst(const ReducedImage &Another,const bool Warn = true) const;

  //! returns true if image is same filter, chip, night, instrument
  bool SameChipFilterInstNight(const ReducedImage &Another,const bool Warn = true) const;

  //! check that (Fits) images have the same sizes
  bool SamePhysicalSize(const ReducedImage &OtherImage) const;

  //! Multiply images by gain when not =1 (after stacking for example)
  double MultiplyGain();

  //!
  virtual ReducedImage* Clone() const;

#ifdef USE_ROOT
#ifndef SWIG
  ClassDef(ReducedImage,1);
#endif /*SWIG */
#endif

  virtual ~ReducedImage();

 private:
  bool has_key(const char *KeyName, const string &RoutineName)const ;
  void remove_key(const char *KeyName, const string &RoutineName);
  int    read_int_key   (const char *KeyName, const string &RoutineName) const;
  double read_double_key(const char *KeyName, const string &RoutineName) const; 

  bool set_double_key(const double &Value, const char *KeyName, const string &RoutineName, 
                      const string Comment);
  bool set_int_key(const int &Value, const char *KeyName, const string &RoutineName, 
                      const string Comment);
  bool set_bool_key(const bool &Value, const char *KeyName, 
		   const string &RoutineName, 
		   const string Comment);
  string read_string_key(const char *KeyName, const string &RoutineName) const; 
  bool set_string_key(const string &Value, const char *KeyName, const string &RoutineName, 
                      const string Comment);
  void init();

  void operator = (const ReducedImage&); // same comment
  
};

typedef CountedRef<ReducedImage> ReducedImageRef;


//! to Reload an already existing ReducedImage
//ReducedImageRef ReducedImageRead(const char *Name);

//! to Reload an already existing ReducedImage
//ReducedImageRef ReducedImageRead(const string &Name);

//! allows to sort a list in increasing seeing order
bool IncreasingSeeing(const ReducedImage* one, const ReducedImage* two); 

//! allows to sort a list in increasing pixel scale order
bool IncreasingResolution(const ReducedImage* one, const ReducedImage* two); 

//! allows to sort a list in increasing julian date order
bool IncreasingJulianDate(const ReducedImage* one, const ReducedImage* two);

//! allows to sort a list in decreasing area
bool DecreasingArea(const ReducedImage* one, const ReducedImage* two);


/***************** ReducedImageList *************************/

#include "imagelist.h"
#include "stringlist.h"

class ReducedImageList : public ImageList<ReducedImage> {
public :
  ReducedImageList(const StringList &Names)
  {
    for (StringCIterator i=Names.begin(); i != Names.end(); ++i)
      push_back(new ReducedImage(*i));
  }
  ReducedImageList(const bool ShouldDelete = true) 
    : ImageList<ReducedImage>(ShouldDelete) {};

  const_iterator Find(const string &Name) const;

};

typedef ReducedImageList::iterator       ReducedImageIterator;
typedef ReducedImageList::const_iterator ReducedImageCIterator;


//! the common pixel frame to all images of the list
Frame CommonFrame(ReducedImageList &RList);

#ifndef SWIG 
// SWIG does not like pointer to member functions
//! Performs the boolean OR of boolean images (Dead, Satur). call the "Make" routine (3rd argument) if the file (2nd argument) of one component is missing 
bool BoolImageOr(ReducedImageList &List,
                 string (ReducedImage::*InFitsFileName)() const,
		 bool (ReducedImage::*MakeInFits)() ,
		 const string OutFitsName);

//! Performs the boolean AND of boolean images (Dead, Satur). call the "Make" routine (3rd argument) if the file (2nd argument) of one component is missing is missing

bool BoolImageAnd(ReducedImageList &List,
		  string (ReducedImage::*InFitsFileName)() const,
		 bool (ReducedImage::*MakeInFits)() ,
		 const string OutFitsName);

#endif


#endif /* REDUCEDIMAGE__H */
