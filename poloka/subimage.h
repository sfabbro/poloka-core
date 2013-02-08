#ifndef SUBIMAGE__H
#define SUBIMAGE__H

#include <string>

#include <poloka/reducedimage.h>

//! a class that allows to cut a subimage from a ReducedImage
class SubImage : public ReducedImage {

 private:
  Frame subFrame; // the coordinates in the large one
  Frame largeFrame; // to store the size of the large image
  string largeImageName;

  bool cut_in_fitsimage(string (ReducedImage::*GetFitsFileName)() const,
			bool (ReducedImage::*MakeInFits)());
 public:
  SubImage(); // hopefully never used

  //! 
  SubImage(const string &Name, 
	   const string &LargeImageName, 
	   const Frame &SubFrame);

  SubImage(const string& Name);
  SubImage(const string &Name, 
	   const string &LargeImageName,
	   const string &SmallImageName,
	   int ExtraMargin=0);

  virtual const string  TypeName() const { return "SubImage";}
  //!
  bool MakeFits();
  //!
  bool MakeWeight();
  //!
  bool MakeCosmic();
  //!
  bool MakeSatellite();
  //! 
  bool MakeSatur();
  //! 
  bool MakeDead();
  //! 
  bool MakeCatalog();
  bool MakeAperCat();
  //! 
  ReducedImageRef Clone() const { return new SubImage(*this);};

};


#endif /* SUBIMAGE__H */
