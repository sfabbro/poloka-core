// This may look like C code, but it is really -*- C++ -*-
#ifndef DAOPHOT__H
#define DAOPHOT__H

#include <string>

#include "fitssubs.h"
#include "daophotio.h"
#include "daophotoption.h"

class ReducedImage;

//! A wrapper to most DAOPHOT routines. More DAOPHOT documentation can be found 
//! in the DAOPHOT II manual, and in the original papers

class Daophot {
  
  int     open;                       // to check whether fits image is open
  float  *data;                       // a pointer to the fits image data
  float   lowbad, threshold, ap1;     // some numbers to write the daophot files
  float   global_sky;                 // background value of the whole image
  string  rootname;                   // the name of the file without extension
  
  Daophot(const Daophot&);            // assignment not properly handled 
  Daophot& operator=(const Daophot&); // copy not properly handled 

public:

  //! empty constructor make sure pointers are pointing to 0
  Daophot() : open(0), data(0) {}

  //! initialize with a ReducedImage, the recommanded option
  Daophot(const ReducedImage &Rim);

  //! free allocated memory and close opened files
  ~Daophot();

  //! Options to set for daophot, see DaophotOptions for use.
  DaophotOptions opt;

  //! open (but does not read) the FITS image -  equivalent to DAOPHOT/ATTACH
  void Attach(const string& FitsFileName);
  void Attach() { Attach(rootname+".fits"); }

  //! compute and print the sky and r.m.s.    -  equivalent to DAOPHOT/SKY
  void Sky();

  //! detect stars with a gaussian filter     -  equivalent to DAOPHOT/FIND
  void Find(const string& CooFileName) const;
  void Find() const { Find(rootname+".coo"); }

  //! perform aperture photometry             -  equivalent to DAOPHOT/PHO
  void Photometry(const string& CooFileName, const string& MagFileName, const string &OptFileName) const;
  void Photometry() const { Photometry(rootname+".coo", rootname+".ap", "photo.opt"); }

  //! select isolated stars for PSF           -  equivalent to DAOPHOT/PICK 
  void Pick(const string& MagFileName, const string& LstFileName, const int Nstar, const float& MagLim) const;
  void Pick(const int Nstar=40, const float& MagLim=16.) const { Pick(rootname+".ap", rootname+".lst", Nstar, MagLim); }
  
  //! fit PSF parameters                      -  equivalent to DAOPHOT/PSF
  void Psf(const string& ApFileName, const string& LstFileName, const string& PsfFileName, const string& NeiFileName);
  void Psf() { Psf(rootname+".ap", rootname+".lst", rootname+".psf", rootname+".nei"); }

  //! fit the PSF to a list of stars          -  equivalent to DAOPHOT/PEAK
  void Peak(const string& MagFileName, const string& PsfFileName, const string& PkFileName) const;
  void Peak() const { Peak(rootname+".ap", rootname+".psf", rootname+".pk"); }

  //! group stars for simultaneous fitting    -  equivalent to DAOPHOT/GROUP
  void Group(const string& MagFileName, const string& PsfFileName, const string& GrpFileName, const float& CritOverlap) const;
  void Group(const float& CritOverlap=0.1) const { Group(rootname+".ap", rootname+".psf", rootname+".grp", CritOverlap); }

  //! simultaneous fit of grouped stars       -  equivalent to DAOPHOT/NSTAR
  void Nstar(const string& GrpFileName, const string& PsfFileName, const string& NstFileName) const;
  void Nstar() const { Nstar(rootname+".grp", rootname+".psf", rootname+".nst"); }

  //! subtract stars from an image            -  equivalent to DAOPHOT/SUB*
  void Substar(const string& PsfFileName, const string& NstFileName, const string& LstFileName, const string& SubPicFileName) const;
  void Substar() const { Substar(rootname+".psf", rootname+".nst", rootname+".lst", rootname+"_sub.fits"); }

  //! add fake stars and noise from a PSF     -  equivalent to DAOPHOT/ADD*
  void Addstar(const string& PsfFileName, const string& AddFileName, const string& AddPicName, const int InSeed);
  void Addstar() { Addstar(rootname+".psf", rootname+".add", rootname+"_add.fits", 10); } 

  //! fit and iteratively subtract stars      -  equivalent to ALLSTAR
  void Allstar(const string& ApFileName, const string& PsfFileName, const string& AlsFileName, const string& SubPicFileName) const;
  void Allstar() const { Allstar(rootname+".ap", rootname+".psf", rootname+".als", rootname+"_sub.fits"); }

  //! dump pixel intensities on screen        -  equivalent to DAOPHOT/DUMP
  void Dump(const Point& Pt, const float& Size) const;

  //! get value of sky with 3 estimators
  void GetSky(float &SkyMean, float &SkyMedian, float &SkyMode, float &SkySigma) const;
  
  //! write a DAOPHOT style star file
  template<DbImageCatalogKind filetype>
  void WriteSEStarList(const SEStarList& Stars) const
  {    
    const string filename = rootname + "." + DaoFileExtension(filetype);
    // cout << " Daophot::WriteSEStarList() : Writing " << filename << endl;
    write_dao<filetype>(filename, SIZE.ncol, SIZE.nrow, lowbad, opt[AduHighDatum].Value(), 
			threshold, ap1, opt[Gain].Value(), 
			opt[ReadNoise].Value(), opt[FitRadius].Value(), Stars);
  }

  //! read a DAOPHOT style star file into a SEStarList
  template<DbImageCatalogKind filetype>
  void ReadSEStarList(SEStarList& Stars)
  {    
    const string filename = rootname + "." + DaoFileExtension(filetype);
    ifstream daostream(filename.c_str());
    /*
      float gain,readnoise,fitrad;
      int ncol,nrow;
      read_dao_header(daostream, ncol, nrow, lowbad, threshold, ap1, gain, readnoise, fitrad);
    */
    string dum;
    read_dao_header(daostream, dum);
    read_dao_starlist<filetype>(daostream, Stars);
    // cout << " Daophot::ReadSEStarList() : Read " << Stars.size() << " stars in " << filename << endl;
  }
};


#endif // DAOPHOT__H
