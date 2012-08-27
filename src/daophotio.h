// This may look like C code, but it is really -*- C++ -*-
#ifndef DAOPHOTIO__H
#define DAOPHOTIO__H

#include <iomanip>
#include <fstream>
#include <map>

#include "sestar.h"
#include "daostar.h"
#include "reducedimage.h"

//! enum for DAO catalogs are created as <Routine><Extension>
enum DaoCatalogEnum {
  FindCoo = 1,
  AllstarAls,
  PhotAp,
  PickLst,
  PsfNei,
  DaoUndef};

//! the option enum to index the vector of options
typedef enum {ReadNoise=0,
	      Gain,
	      SigmaLowDatum,
	      AduHighDatum,
	      Fwhm,
	      SigmaThreshold,
	      LowSharpness,
	      HighSharpness,
	      LowRoundness,
	      HighRoundness,
	      WatchProgress,
	      FitRadius,
	      PsfRadius,
	      VariablePsf,
	      SkyEstimator,
	      AnalyticPsf,
	      ExtraCleanPsf,
	      UseSaturated,
	      PercError,
	      ProfError,
	      ClipExponent,
	      Recentroid,
	      ClipRange,
	      MaxGroup,
	      InnerSky,
	      OutterSky,
	      AperRad01,
	      AperRad02,
	      AperRad03,
	      AperRad04,
	      AperRad05,
	      AperRad06,
	      AperRad07,
	      AperRad08,
	      AperRad09,
	      AperRad10,
	      AperRad11,
	      AperRad12} DaoOptionEnum;


//! a small class to hold options name, limits and value
class DaoOption { 

  string name;     // option name, daophot style
  float  omin;     // min value
  float  omax;     // max value 
  float  odefault; // default value
  float  value;    // current value

public:

  DaoOption() {}

  //! fill up the option
  DaoOption(const string& Name, const float& Min, 
	    const float&  Max,  const float& Default);
 
  //! use this function to set the option value, it checks against limits
  bool SetValue(const float& Value);

  //! the current value of the option
  float& Value() { return value; }

  //! the current value of the option
  const float& Value() const { return value; }

  //! the first two letter of the option name
  string ShortName() const;

  //! check 2 letters against option short name 
  bool operator == (const string& OpName) const { return ShortName() == OpName;}

  //! enable cout << DaoOption << endl;
  friend ostream& operator << (ostream &stream, const DaoOption& opt);

  //! enable cin >> DaoOption
  friend istream& operator >> (istream &stream, DaoOption& opt);
};

typedef map<DaoOptionEnum,DaoOption> DaoOptions;

//! enable dumping of all options, ala daophot 
ostream& operator << (ostream &stream, const DaoOptions& Options);

//! input
istream& operator >> (istream &stream, DaoOptions& Options);

//! stdout DAOPHOT like of verbose definitions and values of options
void verbose_dao_options(const DaoOptions& options);

//! make a decent daophot option file from header
void WriteDaophotOptions(const FitsHeader& Head, const string& FileName="daophot.opt");
void WriteDaophotOptions(const ReducedImage& Im, const string& FileName="daophot.opt");

//! make a decent allstar option file from header
void WriteAllstarOptions(const FitsHeader& Head, const string& FileName="allstar.opt");
void WriteAllstarOptions(const ReducedImage& Im, const string& FileName="allstar.opt");

//! make a decent daophot aperture option file from header
void WriteDaoAperOptions(const FitsHeader& Head, const string& FileName="photo.opt");
void WriteDaoAperOptions(const ReducedImage& Im, const string& FileName="photo.opt");

//! write all the above with default filenames in directory DirName
void WriteDaoOptions(const FitsHeader& Head, const string& DirName=".");
//! write all the above with default filenames in dbimage directory if not specified
void WriteDaoOptions(const ReducedImage& Im, const string& DirName="");

//! the type of star list in the daophot private style files
int DaoFileNumber(const DaoCatalogEnum FileType);

//! the extension to associate with a daophot style file
string DaoFileExtension(const DaoCatalogEnum filetype);

//! the extension to associate with a daophot style file
DaoCatalogEnum DaoFileType(const string& FileName);

//! returns a dao file name Preference order: alf, als, ap, coo.
string FindDaoList(const string& FileBase);

//! read a two line header in daophot star list files
void read_dao_header(istream &daostream, int& Ncol, int& Nrow, float& LowBad, float& Threshold, 
		     float& Ap1, float& Gain, float& ReadNoise, float& FitRad);

//! read and keep a daophot 2 lines header
void read_dao_header(istream &daostream, string& header);

//! reads a DaoStar from a stream, according to the filetype
template<DaoCatalogEnum filetype> 
void read_dao_star(istream &daostream, DaoStar &star);

//! reads and fills in the star list from a stream given a filetype
template<DaoCatalogEnum filetype> 
void read_dao_starlist(istream &daostream, DaoStarList &Stars) {
  char c;
  while (daostream >> c)  { // end-of-file test
    daostream.unget();
    DaoStar *star = new DaoStar;
    double mag;
    read_dao_star<filetype>(daostream, *star); 
    Stars.push_back(star);
  }
}

//! reads and fills in the star list and the header from a stream given a filetype and the name of the file
template<DaoCatalogEnum filetype> 
void read_dao(const string &FileName, DaoStarList &Stars) {
  ifstream daostream(FileName.c_str());
  if (!daostream) {
    cerr << " read_dao : unable to open file '" 
	 << FileName << "'\n";
    return;
  }
  string dummy;
  getline(daostream, dummy);
  getline(daostream, dummy);
  read_dao_starlist<filetype>(daostream, Stars); 
}

//! reads a dao file and guess the type by its suffix
void ReadDaoList(const string& FileName, DaoStarList& Stars);

//! writes a DaoStar to a stream, according to the filetype
template<DaoCatalogEnum filetype> 
void write_dao_star(ostream &daostream, const DaoStar &star);

//! writes a properly formatted daophot style header
template<DaoCatalogEnum filetype> 
void write_dao_header(ostream &daostream, const int Nx, const int Ny, 
		      const float& LowBad, const float& HighBad, 
		      const float& Threshold, const float& Ap1, 
		      const float& Gain, const float&  ReadNoise, const float& FitRad) {
  int nl = DaoFileNumber(filetype);
  daostream << " NL   NX   NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n"; 
  int ndigits=0;
  if (HighBad > 99999.9) ndigits=0;
  daostream << setiosflags(ios::fixed);
  daostream << " " 
	    << setw(2) << nl 
	    << setw(5) << Nx 
	    << setw(5) << Ny 
	    << setw(8) << setprecision(1) << LowBad 
	    << setw(8) << setprecision(ndigits) << HighBad 
	    << setw(8) << setprecision(2) << Threshold 
	    << setw(8) << Ap1 << setw(8) << Gain << setw(8) << ReadNoise << setw(8) << FitRad;
  daostream << endl << endl;
  daostream << resetiosflags(ios::fixed);
}

//! writes a properly formatted daophot style header from a ReducedImage
template<DaoCatalogEnum filetype> 
void write_dao_header(ostream &daostream, const ReducedImage &Rim) {
  write_dao_header<filetype>(daostream, Rim.XSize(), Rim.YSize(),
			     Rim.BackLevel() - 7* Rim.SigmaBack(),
			     Rim.Saturation(), 4*Rim.SigmaBack(), 
			     Rim.Seeing()*2.3548, Rim.Gain(), Rim.ReadoutNoise(), 
			     Rim.Seeing()*2.3548);
}

//! writes a DaoStarList into a stream given a Daophot filetype
template<DaoCatalogEnum filetype> 
void write_dao_starlist(ostream &daostream, const DaoStarList &Stars) {
  daostream << setiosflags(ios::fixed);
  for (DaoStarCIterator it=Stars.begin(); it!=Stars.end(); ++it) {
    write_dao_star<filetype>(daostream, **it);
    daostream << endl;
  }
  daostream << resetiosflags(ios::fixed);
}

void SEStar2DaoStarList(const SEStarList& SEStars, DaoStarList& DaoStars);
void DaoStar2SEStarList(const DaoStarList& DaoStars, SEStarList& SEStars);

//! converts se.list into a daophot star list
template<DaoCatalogEnum filetype> 
void write_dao_starlist(ostream &daostream, const SEStarList &Stars) {
  DaoStarList daostars;
  SEStar2DaoStarList(Stars, daostars);
  write_dao_starlist<filetype>(daostream, daostars);
}

//! writes a daophot starlist and header
template<DaoCatalogEnum filetype> 
void write_dao(const string FileName, const int ncol, const int nrow, const float& lowbad, 
	       const float& highbad, const float& threshold, const float& ap1, const float& gain, 
	       const float& rdnoise, const float& fitrad, const DaoStarList &Stars) {
  ofstream daostream(FileName.c_str());
  write_dao_header  <filetype>(daostream, ncol, nrow, lowbad, highbad, threshold, ap1, gain, rdnoise, fitrad);
  write_dao_starlist<filetype>(daostream, Stars); 
}

//! writes a daophot starlist and creates a non-sense bogus header from the list
template<DaoCatalogEnum filetype> 
void write_dao(const string FileName, const DaoStarList &Stars) {
  double minx=1e30, maxx=-1e30, miny=1e30, maxy=-1e30, mins=1e30;
  for (DaoStarCIterator it = Stars.begin(); it != Stars.end(); ++it) {
    const DaoStar* s;
    if (s->x <= minx) minx = s->x;
    if (s->y <= miny) miny = s->y;
    if (s->x >= maxx) maxx = s->x;
    if (s->y >= maxy) maxy = s->y;
    if (s->sky <= mins) mins = s->sky;
  }
  ofstream daostream(FileName.c_str());
  write_dao_header  <filetype>(daostream, int(maxx-minx), int(maxy-miny),
			       mins-7*sqrt(mins), 65536, 4, 4, 1, 0.1, 4);
  write_dao_starlist<filetype>(daostream, Stars); 
}

//! writes a DaoStarList for a ReducedImage with proper file type
template<DaoCatalogEnum filetype> 
void write_dao(const ReducedImage &Rim, const string& FileName, const DaoStarList &Stars) {
  ofstream daostream(FileName.c_str());
  write_dao_header  <filetype>(daostream, Rim);
  write_dao_starlist<filetype>(daostream, Stars); 
}

//! writes a DaoStarList for a ReducedImage with proper file type
template<DaoCatalogEnum filetype> 
void write_dao(const ReducedImage &Rim, const DaoStarList &Stars) {
  write_dao<filetype>(Rim.Dir() + Rim.Name() + DaoFileExtension(filetype), Rim, Stars);
}

void WriteDaoList(const ReducedImage& Im, const string& FileName, const DaoStarList& Stars);
void WriteDaoList(const string& FileName, const DaoStarList& Stars);

//! merge a SExtractor and DAOPHOT into the SExtractor catalog
void MergeSexDao(SEStarList& Stars, const DaoStarList& DaoStars);
void MergeSexDao(const string& SexName, const string& DaoName);

//! converts lists
void Dao2Sex(const string& DaoName, const string& SexName);
void Sex2Dao(const string& SexName, const string& DaoName);

//! prepare everything to run daophot within a DbImage or in a separate directory
void DaoSetup(ReducedImage& Im, const string& Dir="");

#endif // DAOPHOTIO__H
