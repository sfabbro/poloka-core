// This may look like C code, but it is really -*- C++ -*-

#ifndef FITSMODE__H
#define FITSMODE__H
enum FitsFileMode {RO = 0, RW = 1};
#endif


#ifndef FITSIMAGE__H
#define FITSIMAGE__H

#include <string>
#include <cstdio>
#include "image.h"
#include "point.h"

/*! \file
   \brief I/O of images/headers in fits format

The IO of images in the fits format are carried out using
the cfitsio library. This library is written in C, and the most
useful functionnalities have been wrapped here in C++ classes.
*/


/* copied from fitsio.h , that we try not include in client code */
/*
#define RO 0
#define RW 1
*/





class Frame;
class StringList;
typedef struct _FITSFILE fitsfile;// the pointer provided by cfitsio.


//-----------------------------------------------------------------------------------
//****************************Class FitsKey****************************************

/* this class is only there to add the convenience functionnality of reading the
keys throuth statements like :
int nx = image.KeyVal("NAXIS1");
or 
cout << image.KeyVal("NAXIS1") << endl;

setting & deleting keys is done directly in the FitsImage/ FitsHeader classes.
Since there is no point in building such an object directly,
the constructor is private. So before adding any functionnality to
FitsKey, consider adding it in FitsHeader class. */


class GtransfoLin;
class VirtualInstrument;

#define NOVAL string("NOVAL")

//! Auxilary class for accessing fits header keys

class FitsKey { 

  friend class FitsHeader;

  private :
/* this could become a pointer if cfitsio routines use const when applicable (e.g.fits_read_key) */
      char keyName[32];
      fitsfile *fptr;
      double dval;
      string sval;
  bool warn;
  // useless to call directly:
  FitsKey(const string &KeyName, fitsfile *Fptr, const bool Warn) : fptr(Fptr), warn(Warn){ strncpy(keyName,KeyName.c_str(),32);};

public :
  // constructors without a fits file
      FitsKey(const string &KeyName, const double& val) :  
          fptr(0), dval(val) 
         { strncpy(keyName, KeyName.c_str(), 32); char toto[30]; sprintf(toto,"%f",val); sval = toto;};

      FitsKey(const string &KeyName, const int val) :  
          fptr(0), dval(val) 
         { strncpy(keyName, KeyName.c_str(), 32); char toto[30]; sprintf(toto,"%d",val); sval = toto;};

      FitsKey(const string &KeyName, const string &val) :   
            fptr(0), sval(val) 
         { strncpy(keyName, KeyName.c_str(), 32); dval = atof(sval.c_str());};
      
  // conversion operators
  // the 5 next ones enable {int,double,string} nx = image.KeyVal("NAXIS1");
  operator int() ; /* those operators are not const because fits_read_key  does not pass const arguments */
  operator float() ;
  operator double() ;
  operator string () ;
  operator bool() ;
  string KeyName() const { return string(keyName);};
  friend ostream& operator <<(ostream& stream, const FitsKey &Key);

};

/***********  FitsKeyArray ************/
/**** used to retrieve key batches at once (FitsHeader::KeyMatch)*/

#include <vector>
class FitsKeyArray : public vector<FitsKey>
{
};
  




//-------------------------------------------------------------------------------
//*************************Class FitsHeader*************************************

//! Fits files and header keys.
/*! The FitsHeader class does not deal with the data part of
fits files. The actual engine is cfitsio. */

class FitsHeader {

  friend class FitsImageArray;

  public :

  FitsHeader();
  //! opens the file, by default in readonly mode. 
  FitsHeader(const string &FileName, const FitsFileMode Mode = RO);
  
  //! opens a new file and copies a old header into it. 
  FitsHeader(const FitsHeader &a_header, const string & NewFileName);
 
  //! returns if the file could be opened. 
  bool IsValid() const { return (fptr!=0);};

  //! access routine to the fits file name. 
  string FileName() const { return fileName;}

  
  //! return the file mode (RO or RW)
  FitsFileMode FileMode() const;

  //! read a key. warns on request 
  /*! use (e.g.) \verbatim double ra = a_header.KeyVal("RA"); \endverbatim
     adequate converters are applied for int, double, float, string and bool
     variables.  */
  FitsKey KeyVal(const string &KeyName, const bool Warn = false) const;

  //!modifies an existing key value and comment. 
  /*!  If comment is empty the comment is unchanged.*/
  int ModKey(const string &KeyName, const int Value, const string Comment = "") const ;
  //! -
  int ModKey(const string &KeyName, const double Value, const string Comment = "") const ;
  //! -
  int ModKey(const string &KeyName, const char *Value, const string Comment = "") const ;
  //! -
  int ModKey(const string &KeyName, const bool Value, const string Comment = "") const;

  //! adds (at the end of the header) a new key. 
  /*! Checks before
     that this key does not exist yet. No comment if Comment is NULL or absent.
     Value can be int, double, char*, or bool. */
  int AddKey(const string &KeyName, const char* KeyVal, const string Comment = "");
  int AddKey(const string &KeyName, const double KeyVal, const string Comment = "");
  int AddKey(const string &KeyName, const int KeyVal, const string Comment = "");
  int AddKey(const string &KeyName, const bool KeyVal, const string Comment = "");

  //!  modifies an existing key or adds it if it does not exist yet.
  int AddOrModKey(const string &KeyName, const char  *Value, const string Comment = "");
  //!  modifies an existing key or adds it if it does not exist yet.
  int AddOrModKey(const string &KeyName, const string &Value, 
		  const string Comment = "")
  {return AddOrModKey(KeyName, Value.c_str(), Comment);}

  //!  modifies an existing key or adds it if it does not exist yet.
  int AddOrModKey(const string &KeyName, const double Value, const string Comment = "");
  //!  modifies an existing key or adds it if it does not exist yet.
  int AddOrModKey(const string &KeyName, const int    Value, const string Comment = "");
  //!  modifies an existing key or adds it if it does not exist yet.
  int AddOrModKey(const string &KeyName, const bool   Value, const string Comment = "");

  //! To read keys following a pattern
  /*!  * matches anything, ? a single character, # successive decimal digits.
    Array[i]. size() returns the number of matched keys.
    Array[i].KeyName   and {double,string,int}(Array[i]) enable to accees 
    key names and values */
  int KeyMatch(const string &KeyPattern, FitsKeyArray &Array) const;


  //! changes the Key itself. 
  int ModKeyName(const string &OldKeyName, const string &NewKeyName) const ;
  
  //! changes the comment.
  int ModKeyComment(const string &KeyName, const string &NewComment) const ;

  //! enables to check the presence of a key. warns on request
  bool HasKey(const string &KeyName, const bool Warn=false) const ;
  
  //! has a genuine fits key.
  bool HasActualKey(const string &KeyName, const bool Warn=false) const;

  //! deletes a key 
  int RmKey(const string &KeyName) const ;

  //! return the number of keys (without counting the END key). 
  int NKeys() const;

  //! copy verbatim a key (name, value, comment)
  bool CopyKey(const std::string &KeyName, FitsHeader &To) const;

  //! add a COMMENT keyword. 
  /*! It will be split over multiple 
      COMMENT lines if longer than 70 characters. */
  int AddCommentLine(const string &AVeryUsefulComment);

  //! add a HISTORY keyword. 
  /*! It will be split over multiple 
      HISTORY lines if longer than 70 characters. */
  int AddHistoryLine(const string &HistoryStuff);


  //! flush file buffers.
  int Flush();
  bool ReadCard(const std::string &KeyName, string &Card) const;

  //! add a whole card, name+value+comment (or modifies an existing one)
  int AddOrModCard(const string &KeyName, const string &Card);


  //! dumps the keys contained in WhichKeys to an ascii file.
  bool AsciiDump(const std::string &AsciiFileName, const StringList &WhichKeys) const;



  //!
  void ImageSizes(int &Xsize, int &YSize) const;


  //!
  bool SameImageSizes(const FitsHeader &Other) const;

  //! The geometric center of the image.
  Point ImageCenter() const;


  //! checks that both images refer to the same chip, filter and instrument
  bool SameChipFilterInst(const FitsHeader &Other, const bool Warn = true) const;

  //! checks that both images refer to the same chip and filter.
  bool SameChipFilter(const FitsHeader &Other, const bool Warn = true) const;

  //!
  bool SameChipFilter(const string &OtherFitsName, const bool Warn = true) const;

  //! checks that both images refer to the same filter.
  bool SameFilter(const FitsHeader &Other, const bool Warn = true) const;

  //!
  bool SameFilter(const string &OtherFitsName, const bool Warn = true) const;

  //! checks that both images refer to the same chip.
  bool SameChip(const FitsHeader &Other, const bool Warn = true) const;

  //! same as above, using file name.
  bool SameChip(const string &OtherFitsName, const bool Warn = true) const;


   //! enables : \verbatim cout << FitsHeader("my_image.fits"); \endverbatim
friend ostream& operator << (ostream &stream, const FitsHeader &Header);
  
  VirtualInstrument *TelInst() const;

  //! for a FitsImage opened RW, enables to forbid writing (default behaviour) when destructor is called.

  /*! Example of use : the flatfielding opens (RW) the flatfielded
    FitsImage, and only allows writing on successful flatfielding
    completion. */
  void EnableWrite(const bool YesOrNo);
  
  friend class FitsKey; friend class FitsImage; 
  friend class FitsSlice; friend class FitsOutSlice;

  /* routines used to split multiple image fits files (used in split_fits)*/
  void Append_LowPriority(const FitsHeader& ToAppend);
  int MoveHDU(int HowMany = 1);
  int CopyCHDUTo(FitsHeader &OutHeader);
  int CopyDataTo(FitsHeader &OutHeader);
  int NHDU() const;

  //!
  void DeleteFile();

  // closes file
  ~FitsHeader();


protected:


private:


  int create_file();

  //! reads the image data. 0,0,nx,ny reads the whole image.
  int read_image(const int xmin, const int ymin, const int xmax, const int ymax,
		 float *data);

  string fileName;
  fitsfile *fptr; /* internal pointer for cfitsio */
  bool writeEnabled;
  bool compressedImg;
  VirtualInstrument *telInst; 
  
  FitsFileMode fileModeAtOpen;

 int mod_key(int type, const string &KeyName, void *Value, const string &Comment) const;
 int write_key(const string &KeyName, void* KeyVal, const string &Comment, const int type);  /* the generic one */
 void minimum_header();
  
/* disable for the moment the copy, because it provides a way to create a new header without providing a filename*/
#ifndef SWIG
  FitsHeader(const FitsHeader &);
  FitsHeader& operator=(const FitsHeader &Right);
#endif

};


//---------------------------------------------------------------------------------------
//****************************Class FitsImage********************************************

//!This class enables basic manipulation of images stored in fits files.
/*! For The image itself,
see the Image class. The files opened in mode RW are actually written to disk
by the destructor. FitsHeader::EnableWrite() allows to alter this default behavior.
*/ 
class FitsImage : public FitsHeader, public Image {

  public :

  //!  opens in (ReadOnly mode by default, use RW to modifiy or create) and loads the pixel data.
FitsImage(const string &FileName, const FitsFileMode Mode = RO);


  //! constructor for a new FitsImage.
  /*! To be used to assemble a FitsImage from a header and an image
        obtained separately. The associated mode is RW. The actual size
        of the saved image is the one of the image, not the one in the header.*/ 
FitsImage(const string &NewFileName, const FitsHeader &a_fits_header, const Image & an_image);

  //! constructor of a new FitsImage using an existing header, a new name, and a zeroed Image. RW mode.
FitsImage(const string &NewFileName, const FitsHeader &AHeader);

//! constructor for a new FitsImage with a minimal header from an existing image. 
FitsImage(const string &FileName, const Image& an_image);

//! constructor for a new FitsImage with a minimal header from scratch. 
  FitsImage(const string &FileName, const int Nx, const int Ny);

//! Extract and replaces the trimed Image into a FitsImage 
void Trim(const Frame &Region);

  //! This routine is to be called if one wants to preserve the value of 0 on I/O. Assumes contents >=0, as it is Used for weight maps.
  void PreserveZeros();


~FitsImage();

/* This image will be saved as floating point data when written. */
bool SetWriteAsFloat();


  int Write(bool force_bscale = false);
  int Write(const double &Bscale, const double &Bzero);


 private :
   int written ;

};


#endif /* FITSIMAGE__H */

