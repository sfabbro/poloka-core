// This may look like C code, but it is really -*- C++ -*-
#ifndef WCSUTILS__H
#define WCSUTILS__H

#include <string>

class FitsHeader;
class GtransfoLin;
class TanPix2RaDec;
class Gtransfo;
class Frame;
class StringList;

/*! \file
    \brief World Coordinate Transfo  read and write

    Only supports linear (RA--TAN) WCS now (until somebody understands any WCS documentation).
*/

//! returns a transformation that combines availbale information in the header
bool WCSFromHeader(const FitsHeader &Head, Gtransfo* &Pix2RaDec);

//! same as above
bool WCSFromHeader(const std::string &FitsName, Gtransfo* &Pix2RaDec);



//! Fills the fits file header with provided WCS transfo 
int WCSTransfo2Header(const std::string &FitsImageName, const GtransfoLin &Pix2RaDec);

//! 
int WCSTransfo2Header(FitsHeader &header, const GtransfoLin &Pix2RaDec);

//! retrieves the (lin) WCS transfo from a header. 
bool WCSLinTransfoFromHeader(const std::string &FitsImageName, GtransfoLin &Pix2RaDec);

//! 
bool WCSLinTransfoFromHeader(const FitsHeader& Header, GtransfoLin &Pix2RaDec);


/*********** tangent point WCS ***********/

//! return the tangent point WCS transfo.
bool TanLinWCSFromHeader(const FitsHeader &Head, TanPix2RaDec &TanWcs, 
			 const bool Warn= false);

//!read it
bool TanWCSFromHeader(const FitsHeader &Head, TanPix2RaDec &TanWcs, 
		      const bool Warn = false);

//! write it
int TanWCS2Header(const std::string &FitsImageName, const TanPix2RaDec &TanWcs);

//! write it
void TanWCS2Header(FitsHeader &Head, const TanPix2RaDec &TanWcs);

/************/


//! retrieves the (lin) transfo from header1 to header2, using WCS
bool WCSTransfoBetweenHeader(const FitsHeader &Header1, const FitsHeader &Header2, GtransfoLin &Transfo1to2);

//! check the presence of World Coordinate System data 
bool HasLinWCS(const FitsHeader &Header);

//! remove entirely a WCS transfo from a header
// void RemoveWCS(FitsHeader &Header);

//! returns key names that describe the WCS in this header
void WCSKeyList(const FitsHeader &Head, StringList &KeyNames);

//! copy entirely a WCS transfo from a header to another header
bool CopyWCS(const FitsHeader &FromHeader, FitsHeader &ToHeader);

//! fetch pixel sizes in arcsecond from a header using WCS info
void GetPixelSize(const FitsHeader& Head, double &SizeX, double &SizeY);

//! Get Ra and DEC from WCS on center of image
void RaDecFromWCS(const FitsHeader &Header, double &Ra, double &Dec);

//! so far returns average pixel size in arcsecond 
double PixelSize(const FitsHeader& Head);

//!returns the area in arcmin^2 of a frame in a FitsHeader using WCS info
double Arcmin2Area(const Frame &aFrame,const FitsHeader &Header);

//!returns the area in arcmin^2 of the default Frame in a header using WCS info
double Arcmin2Area(const FitsHeader &Header);

//! returns roughly the overlap area between two FitsHeader in arcmin^2 using WCS info
double Arcmin2Overlap(const FitsHeader& Head1,const FitsHeader& Head2);


//! update RA, DEC and EPOCH based on WCS transfo and Frame in header
bool UpdateRaDec(FitsHeader &Header);

bool GuessLinWCS(const FitsHeader &Header, TanPix2RaDec &Guess);

Frame SkyRegion(const FitsHeader &Header);

// check that the WCS is more accurate than the default one.
bool  CheckWCSIsAccurate( const FitsHeader & header );


#endif /* WCSUTILS__H */

