#ifndef FITSTOAD__H
#define FITSTOAD__H


#include <string>
#include <alltelinst.h>

class TanPix2RaDec;
class FitsHeader;
class VirtualInstrument;

//! The FITS key translator for TOADkeys.
FitsKey ToadsKeyVal(const FitsHeader &Head, const string &KeyName, const bool Warn=false);

//! A name for the telescope
string TelescopeName(const FitsHeader &Head);


//! A Name for the Instrument
string InstrumentName(const FitsHeader &Head);

//! A Name for the Tel/Instrument
string TelInstName(const FitsHeader &Head);


//! test if Head is of kind "Inst" (where Inst is one of the virtual instruments)
template<class Inst> bool IsOfKind(const FitsHeader &Head);


//! returns the overscan region associated with the given amplifier
Frame OverscanRegion(const FitsHeader &Head, const int iAmp);

//! returns the total illuminated region of image
Frame TotalIlluRegion(const FitsHeader &Head);

//! returns the illuminated region associated with the given amplifier
Frame IlluRegion(const FitsHeader &Head, const int Iamp);

//! returns the illuminated region associated with the given amplifier once the image is trimed
Frame AmpRegion(const FitsHeader &Head, const int Iamp);

  //! returns the gain corresponding to given amplifier
double AmpGain(const FitsHeader &Head, const int Iamp);

//! the routine called by FitsHeader to identify which Tel/Inst the header belongs to
VirtualInstrument  *SniffTelInst(const FitsHeader &Head);

// necessary because the class is not in this header file.
void VirtualInstrumentDestructor(VirtualInstrument *);

// define WCS test function 
int SetTestWCS(void (*toto)(const FitsHeader &Head));

#endif /* FITSTOAD__H */
