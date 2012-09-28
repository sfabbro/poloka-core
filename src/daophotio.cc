#include <iterator>
#include <algorithm>

#include "fitsimage.h"
#include "fastfinder.h"
#include "sestar.h"
#include "daophotio.h"
#include "reducedutils.h"
#include "photoratio.h"
#include "gtransfo.h"
#include "sextractor_box.h"

int DaoFileNumber(const DaoCatalogEnum FileType) {
  switch (FileType) {      
  case FindCoo:
  case AllstarAls: return 1;
  case PhotAp:  return 2;
  case PickLst: 
  case PsfNei: return 3;
  default: cerr << " DaoFileNumber(" << FileType << ") Error:  file type not implemented \n";
  }
  return 0;
}

string DaoFileExtension(const DaoCatalogEnum FileType) {
  switch (FileType) {
  case AllstarAls: return "als";
  case PsfNei: return "nei";
  case PickLst: return "lst";
  case PhotAp : return  "ap";
  case FindCoo: return "coo";
  default: cerr << " DaoExtension(" << FileType << ") : Error : file type not implemented\n";
  }
  return "";
}

string FindDaoList(const string& FileBase) {
  string filename = FileBase + ".alf";
  if (FileExists(filename))
    return filename;
  filename = FileBase + ".als";
  if (FileExists(filename))
    return filename;
  filename = FileBase + ".ap";
  if (FileExists(filename))
    return filename;
  filename = FileBase + ".coo";
  if (FileExists(filename))
    return filename;
  return "";
}

DaoCatalogEnum DaoFileType(const string& FileName) {
  string ext = FileExtension(FileName);
  if (ext == "als" || ext == "alf" || ext == "pk" || ext == "nst") 
    return AllstarAls;
  if (ext == "coo") return FindCoo;
  if (ext == "nei") return PsfNei;
  if (ext == "lst") return PickLst;
  if (ext ==  "ap") return PhotAp;
  cerr << " DaoFileType: error: can't guess daophot file type for " << FileName << endl;
  return DaoUndef;
}

ostream& operator << (ostream &stream, const DaoOption& opt) {
  size_t oldp = stream.precision();
  stream << setiosflags(ios::fixed);
  stream << string(27-opt.name.length(),' ')  << opt.name
	 << " =" << setw(9) << setprecision(2) << opt.value;
  stream << resetiosflags(ios::fixed) << setprecision(oldp);
  return stream;
}

ostream& operator << (ostream &stream, const pair<DaoOptionEnum,DaoOption>& opt) {
  stream << opt.second.ShortName() << "=" << opt.second.Value();
  return stream;
}

istream& operator >> (istream &stream, DaoOption& opt) {
  string line;
  if (!getline(stream, line)) return stream;
  vector<string> s;
  DecomposeString(s, line, "=");
  if (s.size() != 2) {
    cerr << " DaoOption : error reading line '" << line << endl;
    return stream;
  }
  opt.name = s[0];
  opt.value = atof(s[1].c_str());
  return stream;
}

DaoOption::DaoOption(const string& Name, const float& Min, 
		     const float& Max, const float& Default)
  : name(Name), omin(Min), omax(Max), odefault(Default), value(Default) {}

bool DaoOption::SetValue(const float& Value) {
  if (Value >= omin && Value <= omax) { value=Value; return true;}
  cout << " DaoOption::SetValue() : Replacing " << ShortName() << "=" << Value
       << " by default value = " << odefault << endl;
  value = odefault;
  return false;
}

string DaoOption::ShortName() const  {
  return name.substr(0,2);
}

static void init_daophot_opt(DaoOptions& opt) {
  opt[ReadNoise]      = DaoOption("READ NOISE (ADU; 1 frame)" , 1e-20, 1e20, 2.   );
  opt[Gain]           = DaoOption("GAIN (e-/ADU; 1 frame)"    , 1e-20, 1e20, 1.   );
  opt[SigmaLowDatum]  = DaoOption("LOW GOOD DATUM (in sigmas)", 0.   , 1e20, 7.   );
  opt[AduHighDatum]   = DaoOption("HIGH GOOD DATUM (in ADU)"  , 0.   , 1e20, 1.5e5);
  opt[Fwhm]           = DaoOption("FWHM OF OBJECT"            , 0.2  , 20. , 2.5  );
  opt[SigmaThreshold] = DaoOption("THRESHOLD (in sigmas)"     , 0.   , 1e20, 4.   );
  opt[LowSharpness]   = DaoOption("LS (LOW SHARPNESS CUTOFF)" , 0.   , 1.  , 0.2  );
  opt[HighSharpness]  = DaoOption("HS (HIGH SHARPNESS CUTOFF)", 0.6  , 2.  , 1.   );
  opt[LowRoundness]   = DaoOption("LR (LOW ROUNDNESS CUTOFF)" , -2.  , 0.  , -1.  );
  opt[HighRoundness]  = DaoOption("HR (HIGH ROUNDNESS CUTOFF)", 0.   , 2.  , 1.   );
  opt[WatchProgress]  = DaoOption("WATCH PROGRESS"            , -2.  , 2.  , 0.   );
  opt[FitRadius]      = DaoOption("FITTING RADIUS"            , 1.   , 30. , 2.5  );
  opt[PsfRadius]      = DaoOption("PSF RADIUS"                , 5.   , 51. , 11.  );
  opt[VariablePsf]    = DaoOption("VARIABLE PSF"              , -1.5 , 3.5 , -1.  );
  opt[SkyEstimator]   = DaoOption("SKY ESTIMATOR"             , -0.5 , 3.5 , 0.   );
  opt[AnalyticPsf]    = DaoOption("ANALYTIC MODEL PSF"        , -6.5 , 6.5 , 3.   );
  opt[ExtraCleanPsf]  = DaoOption("EXTRA PSF CLEANING PASSES" , 0.   , 9.5 , 5.   );
  opt[UseSaturated]   = DaoOption("USE SATURATED PSF STARS"   , 0.   , 1.  , 0.   );
  opt[PercError]      = DaoOption("PERCENT ERROR (in %)"      , 0.   , 100., 0.75 );
  opt[ProfError]      = DaoOption("PROFILE ERROR (in %)"      , 0.   , 100., 3.   );
}

static void init_allstar_opt(DaoOptions& opt) {
  opt[FitRadius]      = DaoOption("FITTING RADIUS"            , 1.   , 30. , 2.5  );
  opt[ClipExponent]   = DaoOption("CE (CLIPPING EXPONENT)"    , 0.   , 8.  , 6.   );
  opt[Recentroid]     = DaoOption("REDETERMINE CENTROIDS"     , 0.   , 1.  , 1.   );
  opt[ClipRange]      = DaoOption("CR (CLIPPING RANGE)"       , 0.   , 10. , 2.5  );
  opt[WatchProgress]  = DaoOption("WATCH PROGRESS"            , -2.  , 2.  , 0.   );
  opt[MaxGroup]       = DaoOption("MAXIMUM GROUP SIZE"        , 1.   , 100., 70.  );
  opt[PercError]      = DaoOption("PERCENT ERROR (in %)"      , 0.   , 100., 0.75 );
  opt[ProfError]      = DaoOption("PROFILE ERROR (in %)"      , 0.   , 100., 3.   );
  opt[InnerSky]       = DaoOption("IS (INNER SKY RADIUS)"     , 0.   , 35. , 0.   );
  opt[OutterSky]      = DaoOption("OS (OUTER SKY RADIUS)"     , 0.   , 100., 0.   );
}

static void init_aper_opt(DaoOptions& opt) {
  opt[AperRad01]      = DaoOption("A1  RADIUS OF APERTURE  1" , 1e-30, 1e30, 0.   );
  opt[AperRad02]      = DaoOption("A2  RADIUS OF APERTURE  2" , 0.   , 1e30, 0.   );
  opt[AperRad03]      = DaoOption("A3  RADIUS OF APERTURE  3" , 0.   , 1e30, 0.   );
  opt[AperRad04]      = DaoOption("A4  RADIUS OF APERTURE  4" , 0.   , 1e30, 0.   );
  opt[AperRad05]      = DaoOption("A5  RADIUS OF APERTURE  5" , 0.   , 1e30, 0.   );
  opt[AperRad06]      = DaoOption("A6  RADIUS OF APERTURE  6" , 0.   , 1e30, 0.   );
  opt[AperRad07]      = DaoOption("A7  RADIUS OF APERTURE  7" , 0.   , 1e30, 0.   );
  opt[AperRad08]      = DaoOption("A8  RADIUS OF APERTURE  8" , 0.   , 1e30, 0.   );
  opt[AperRad09]      = DaoOption("A9  RADIUS OF APERTURE  9" , 0.   , 1e30, 0.   );
  opt[AperRad10]      = DaoOption("AA  RADIUS OF APERTURE 10" , 0.   , 1e30, 0.   );
  opt[AperRad11]      = DaoOption("AB  RADIUS OF APERTURE 11" , 0.   , 1e30, 0.   );
  opt[AperRad12]      = DaoOption("AC  RADIUS OF APERTURE 12" , 0.   , 1e30, 0.   );
  opt[InnerSky]       = DaoOption("IS       INNER SKY RADIUS" , 0.   , 1e30, 0.   );
  opt[OutterSky]      = DaoOption("OS       OUTER SKY RADIUS" , 0.   , 2e30, 1.   );
}

ostream& operator << (ostream &stream, const DaoOptions& options) {
  copy(options.begin(), options.end(), ostream_iterator< pair<DaoOptionEnum,DaoOption> >(stream, "\n"));  
  return stream;
}

void write_dao_options(const DaoOptions& options, const string& FileName) {
  ofstream optstream(FileName.c_str());
  optstream << options;
  optstream.close();
}

void verbose_dao_options(const DaoOptions& options) {
  for (map<DaoOptionEnum,DaoOption>::const_iterator it = options.begin(); it != options.end(); ) {
    cout << it->second << "    " ; ++it;
    if (it != options.end())
      cout << it->second << endl; ++it;
  }
}

void WriteDaophotOptions(const FitsHeader& Head, const string& FileName) {
  DaoOptions opt;
  init_daophot_opt(opt);
  opt[ReadNoise].SetValue(Head.KeyVal("TOADRDON"));
  opt[Gain].SetValue(Head.KeyVal("TOADGAIN"));
  if (Head.HasKey("SATURLEV"))
      opt[AduHighDatum].SetValue(Head.KeyVal("SATURLEV"));
  if (Head.HasKey("SESEEING")) {
    float fwhm = Head.KeyVal("SESEEING");
    fwhm *=  2.3548;
    opt[Fwhm].SetValue(fwhm);
    opt[FitRadius].SetValue(fwhm);
    opt[PsfRadius].SetValue(4.0 * fwhm);
  }
  opt[WatchProgress].SetValue(-2);
  write_dao_options(opt, FileName);
}

void WriteDaophotOptions(const ReducedImage& Im, const string& FileName) {
  FitsHeader head(Im.FitsName());
  WriteDaophotOptions(head, FileName);
}

void WriteAllstarOptions(const FitsHeader& Head, const string& FileName) {
  DaoOptions opt;
  init_allstar_opt(opt);
  if (Head.HasKey("SESEEING")) {
    float fwhm = Head.KeyVal("SESEEING");
    fwhm *=  2.3548;
    opt[FitRadius].SetValue(fwhm);
    opt[InnerSky].SetValue(1.5 * fwhm);
    opt[OutterSky].SetValue(7.0 * fwhm);
  }
  opt[WatchProgress].SetValue(0);
  write_dao_options(opt, FileName);
}

void WriteAllstarOptions(const ReducedImage& Im, const string& FileName) {
  FitsHeader head(Im.FitsName());
  WriteAllstarOptions(head, FileName);
}

void WriteDaoAperOptions(const FitsHeader& Head, const string& FileName) {
  DaoOptions opt;
  init_aper_opt(opt);
  if (Head.HasKey("SESEEING")) {
    float fwhm = Head.KeyVal("SESEEING");
    fwhm *=  2.3548;
    opt[AperRad01].SetValue(fwhm);
    opt[AperRad02].SetValue(fwhm*1.2);
    opt[InnerSky].SetValue(1.5 * fwhm);
    opt[OutterSky].SetValue(7.0 * fwhm);
  }
  write_dao_options(opt, FileName);
}

void WriteDaoAperOptions(const ReducedImage& Im, const string& FileName) {
  FitsHeader head(Im.FitsName());
  WriteDaoAperOptions(head, FileName);
}

void WriteDaoOptions(const FitsHeader& Head, const string& DirName) {
  WriteDaophotOptions(Head, DirName+"/daophot.opt");
  WriteDaoAperOptions(Head, DirName+"/photo.opt");
  WriteAllstarOptions(Head, DirName+"/allstar.opt");
}

void WriteDaoOptions(const ReducedImage& Im, const string& DirName) {
  FitsHeader head(Im.FitsName());
  if (DirName.empty())
    WriteDaoOptions(head, Im.Dir());
  else
    WriteDaoOptions(head, DirName);
}

//! the shift in pixels between DAOPHOT coordinate system and TOADS one.
const double DAOPHOT_TOADS_SHIFT = -1.;

//! the zeropoint used by DAOPHOT when computing instrumental magnitudes
const double DAOPHOT_ZP = 25.;

static void read_dao_basestar(istream& daostream, DaoStar& star) {
  static double mag;
  daostream >> star.num >> star.x >> star.y >> mag;
  // DAOPHOT to TOADS conversion x,y and flux
  star.x += DAOPHOT_TOADS_SHIFT; 
  star.y += DAOPHOT_TOADS_SHIFT;
  star.flux = pow(10., (0.4*(DAOPHOT_ZP - mag)));
}

// partial specialization for daophot find file
template<> void read_dao_star<FindCoo>(istream &daostream, DaoStar &star) {
  read_dao_basestar(daostream, star);
  double dy;
  daostream >> star.sharp >> star.round >> dy;
}

// partial specialization for daophot aperture file
template<> void read_dao_star<PhotAp>(istream &daostream, DaoStar &star) {
  read_dao_basestar(daostream, star);
  double skyerror, skyskew, emag;
  daostream >> star.sky >> star.esky >> star.skyskew >> emag;
  star.eflux = 0.921034 * emag * star.flux;
}

// partial specialization for daophot pick file
template<> void read_dao_star<PickLst>(istream &daostream, DaoStar &star) {
  read_dao_basestar(daostream, star);
  daostream >> star.sky;
}

// partial specialization for daophot neighbor file
template<> void read_dao_star<PsfNei>(istream &daostream, DaoStar &star) {
  read_dao_basestar(daostream, star);
  daostream >> star.sky >> star.chi;
  char flag;
  daostream >> flag;
  if (isdigit(flag)) daostream.unget();
  else star.FlagAsBlended(); 
}

// partial specialization for daophot allstar file
template<> void read_dao_star<AllstarAls>(istream &daostream, DaoStar &star) {
  read_dao_basestar(daostream, star);
  double emag;
  daostream >> emag;
  star.eflux = 0.921034 * emag * fabs(star.flux);
  daostream >> star.sky;
  double iter;
  daostream >> iter;
  star.iter = int(iter);
  daostream >> star.chi >> star.sharp;
}

static void write_dao_basestar(ostream& daostream, const DaoStar& star) {
  double mag = (star.flux > 0)? -2.5*log10(star.flux) + DAOPHOT_ZP : DAOPHOT_ZP;
  size_t oldp = daostream.precision();
  daostream << setiosflags(ios::fixed)
	    << setw(6) << star.num
	    << setprecision(3) 
	    << setw(9) << star.x - DAOPHOT_TOADS_SHIFT 
	    << setw(9) << star.y - DAOPHOT_TOADS_SHIFT
	    << setw(9) << mag
	    << setprecision(oldp);
}
// partial specialization for daophot find file
template<> void write_dao_star<FindCoo>(ostream &daostream, const DaoStar &star) {
  write_dao_basestar(daostream, star);
  daostream << setprecision(3) 
	    << setw(9) << star.sharp 
	    << setw(9) << star.round
	    << setw(9) << 0; //that last one is a dummy to be compatible
}

// partial specialization for daophot aperture file
template<> void write_dao_star<PhotAp>(ostream &daostream, const DaoStar &star) {
  write_dao_basestar(daostream, star);
  double eflux = (star.eflux > 0) ? star.eflux : 10000;
  double emag = min(1.0857 * eflux / fabs(star.flux), 9.9999);
  daostream << endl << "    " << setw(9) << star.sky 
	    << setprecision(2) << setw(6) << 0.5 << setw(6) << 0.5  // bogus errors and skewness on sky
	    << setprecision(4) << setw(9) << emag << endl;
}

// partial specialization for daophot peak file
template<> void write_dao_star<PickLst>(ostream &daostream, const DaoStar &star) {
  write_dao_basestar(daostream, star);
  int ndigits=3;
  if (star.sky > 9999.999) ndigits=2;
  daostream << setw(9) << setprecision(ndigits) << star.sky;
}

// partial specialization for daophot neighbor file
template<> void write_dao_star<PsfNei>(ostream &daostream, const DaoStar &star) {
  write_dao_basestar(daostream, star);
  char flag = ' ';
  if (star.HasNeighbours()) flag = '*';
  daostream << setw(9) << star.sky
	    << setw(9) << setprecision(3) << star.chi
	    << "   " << flag;
}

// partial specialization for daophot allstar file
template<> void write_dao_star<AllstarAls>(ostream &daostream, const DaoStar &star) {
  write_dao_basestar(daostream, star);
  int ndigits;
  if (star.sky > 0) ndigits = min(3, max(0, 6-int(log10(star.sky))));
  else ndigits = min(3, max(0, 5-int(log10(-star.sky+0.001))));
  double eflux = (star.eflux > 0) ? star.eflux : 10000;
  double emag = min(1.0857 * eflux / fabs(star.flux), 9.9999);
  daostream << setw(9) << setprecision(4) << emag
	    << setw(9) << setprecision(ndigits) << star.sky
	    << setw(9) << setprecision(0) << star.iter
	    << setw(9) << setprecision(3) << star.chi
	    << setw(9) << star.sharp;
}

void read_dao_header(istream &daostream, int& Ncol, int& Nrow, float& LowBad, float& Threshold, 
		     float& Ap1, float& Gain, float& ReadNoise, float& FitRad) {
  string dummy;
  getline(daostream, dummy);
  // have not seen an utility of nl yet, so do not pass it as argument
  int nl;
  daostream >> nl >> Ncol >> Nrow >> LowBad >> Threshold >> Ap1 >> Gain >> ReadNoise >> FitRad;
}

void read_dao_header(istream &daostream, string& header) {
  getline(daostream, header);
  string dummy;
  getline(daostream, dummy);
  header += "\n" + dummy;
}

void ReadDaoList(const string& FileName, DaoStarList& Stars) {
  DaoCatalogEnum filetype = DaoFileType(FileName);
  switch (filetype) {
  case AllstarAls: return read_dao<AllstarAls>(FileName, Stars);
  case FindCoo: return read_dao<FindCoo>(FileName, Stars);
  case PsfNei: return read_dao<PsfNei>(FileName, Stars);
  case PickLst: return read_dao<PickLst>(FileName, Stars);
  case PhotAp: return read_dao<PhotAp>(FileName, Stars);
  default:
    cerr << " ReadDaoList: no file type found for " << FileName << endl;
  }
}

void WriteDaoList(const string& FileName, const DaoStarList& Stars) {
  DaoCatalogEnum filetype = DaoFileType(FileName);
  switch (filetype) {
  case AllstarAls: return write_dao<AllstarAls>(FileName, Stars);
  case FindCoo: return write_dao<FindCoo>(FileName, Stars);
  case PsfNei: return write_dao<PsfNei>(FileName, Stars);
  case PickLst: return write_dao<PickLst>(FileName, Stars);
  case PhotAp: return write_dao<PhotAp>(FileName, Stars);
  default:
    cerr << " WriteDaoList: no file type found for " << FileName << endl;
  }
}

void WriteDaoList(const ReducedImage& Im, const string& FileName, const DaoStarList& Stars) {
  DaoCatalogEnum filetype = DaoFileType(FileName);
  switch (filetype) {
  case AllstarAls: return write_dao<AllstarAls>(Im, Im.Dir()+"/"+FileName, Stars);
  case FindCoo: return write_dao<FindCoo>(Im, Im.Dir()+"/"+FileName, Stars);
  case PsfNei: return write_dao<PsfNei>(Im, Im.Dir()+"/"+FileName, Stars);
  case PickLst: return write_dao<PickLst>(Im, Im.Dir()+"/"+FileName, Stars);
  case PhotAp: return write_dao<PhotAp>(Im, Im.Dir()+"/"+FileName, Stars);
  default:
    cerr << " WriteDaoList: no file type found for " << FileName << endl;
  }
}

void dao2sex(const DaoStar& daostar, SEStar& sexstar) {
  sexstar.x = daostar.x;
  sexstar.y = daostar.y;
  sexstar.flux = daostar.flux;
  sexstar.EFlux() = daostar.eflux;
  sexstar.Fond() = daostar.sky;
  sexstar.Iter() = daostar.iter;
  sexstar.Chi() = daostar.chi;
  sexstar.Sharp() = daostar.sharp;
  if (daostar.IsSaturated()) sexstar.FlagAsSaturated();
}

void sex2dao(const SEStar& sexstar, DaoStar& daostar) {
  daostar.x = sexstar.x;
  daostar.y = sexstar.y;
  daostar.flux = sexstar.flux;
  daostar.eflux = sexstar.EFlux();
  daostar.num = sexstar.N();
  daostar.sky = sexstar.Fond();
  daostar.iter = sexstar.Iter();
  daostar.chi = sexstar.Chi();
  daostar.sharp = sexstar.Sharp();
  if (sexstar.IsSaturated()) daostar.FlagAsSaturated();
}

void DaoStar2SEStarList(const DaoStarList& DaoStars, SEStarList& SEStars) {
  for (DaoStarCIterator it = DaoStars.begin(); it != DaoStars.end(); ++it) {
    SEStar* sestar = new SEStar();
    dao2sex(**it, *sestar);
    SEStars.push_back(sestar);
  }
}


void SEStar2DaoStarList(const SEStarList& SEStars, DaoStarList& DaoStars) {
  for (SEStarCIterator it = SEStars.begin(); it != SEStars.end(); ++it) {
    DaoStar* daostar = new DaoStar();
    sex2dao(**it, *daostar);
    DaoStars.push_back(daostar);
  }
}

void MergeSexDao(SEStarList& Stars, const DaoStarList& DaoStars) {

  FastFinder finder(*Dao2Base(&DaoStars));
  DaoStarList leftOver;
  DaoStars.CopyTo(leftOver);
  for (SEStarIterator it = Stars.begin(); it != Stars.end(); ++it) {
    SEStar *star = *it;
    const DaoStar* daostar = (const DaoStar *) finder.FindClosest(*star, 0.3);
    if (daostar)  { 
      dao2sex(*daostar, *star);
      leftOver.FindAndRemoveClosest(daostar->x, daostar->y);
    }
  }

  size_t counter = Stars.size();
  for (DaoStarIterator it = leftOver.begin(); it != leftOver.end(); ++it) {
    SEStar* star = new SEStar();
    star->N() = ++counter;
    dao2sex(**it, *star);
    Stars.push_back(star);
  }
}

void MergeSexDao(const string& SexName, const string& DaoName) {
  SEStarList sexList(SexName);
  DaoStarList daoList;
  ReadDaoList(DaoName, daoList);
  MergeSexDao(sexList, daoList);
  sexList.write(SexName);
}

void Dao2Sex(const string& DaoName, const string& SexName) {
  SEStarList sexstars;
  DaoStarList daostars;
  ReadDaoList(DaoName, daostars);
  DaoStar2SEStarList(daostars, sexstars);
  sexstars.write(SexName);
}

void Sex2Dao(const string& SexName, const string& DaoName) {
  SEStarList sexstars(SexName);
  DaoStarList daostars;
  SEStar2DaoStarList(sexstars, daostars);    
  WriteDaoList(DaoName, daostars);
}

void DaoSetup(ReducedImage& Im, const string& Dir) {
  
  string daoname = Dir.empty() ? "daoimage" : Im.Name();

  // replace . with _, daophot is not robust with it
  for(int i = 0; i < daoname.length(); i++)
    if (daoname[i] == '.') daoname[i] = '_';

  string daodir  = Dir.empty() ? Im.Dir() : Dir;
  if (!IsDirectory(daodir)) MKDir(daodir.c_str());
  string fullname = daodir + "/" + daoname;
  cout << " DaoSetup: set-up " << Im.Name() << " in " << daodir << endl;
  float sky = 0;
  bool backsub = Im.BackSub();
  if (backsub) sky = Im.OriginalSkyLevel();

  FitsImage im(Im.FitsName());

  if (!FileExists(fullname + ".fits")) {
    // daophot does not read compressed files
    // and has problem decoding some cfht 16 bits image
    FitsImage daoim(fullname + ".fits", im, im);
    daoim.SetWriteAsFloat();
    // add sky to fits if subtracted
    if (backsub) {
      if (Im.HasMiniBack()) {
	cout << " DaoSetup: re-adding subtracted sky to " << Im.Name() << endl;
	AddMiniBack(daoim, Im.FitsMiniBackName());
      } else { // imagesum and others
	float sky,rms;
	if (daoim.HasKey("SKYSIGEX"))
	  rms = daoim.KeyVal("SKYSIGEX");
	else if  (daoim.HasKey("SEXSIGMA"))
	  rms = daoim.KeyVal("SEXSIGMA");
	else
	  daoim.SkyLevel(&sky, &rms);
	daoim += rms*rms;
      }
      daoim.AddOrModKey("SATURLEV", Im.Saturation() + sky,
			"Saturation level corrected from sky subtraction");
      if (daoim.HasKey("BACK_SUB")) daoim.RmKey("BACK_SUB");
      if (daoim.HasKey("BACKLEV")) daoim.RmKey("BACKLEV");
      if (daoim.HasKey("SEXSKY")) daoim.RmKey("SEXSKY");
      if (daoim.HasKey("SEXSIGMA")) daoim.RmKey("SEXSIGMA");
    }
  }
  
  if (!FileExists(fullname + ".ap") && Im.HasCatalog()) {
    SEStarList stars(Im.CatalogName());
    // add sky to stars as well
    if (backsub)
      for (SEStarIterator it = stars.begin(); it != stars.end(); ++it)
	(*it)->Fond() += sky;
    ofstream daostream((fullname + ".ap").c_str());
    write_dao_header<PhotAp>(daostream, Im);
    write_dao_starlist<PhotAp>(daostream, stars); 
  }

  if (!FileExists(daodir + "/daophot.opt"))
    WriteDaophotOptions(im, daodir+"/daophot.opt");

  if (!FileExists(daodir + "/photo.opt"))
    WriteDaoAperOptions(im, daodir+"/photo.opt");

  if (!FileExists(daodir + "/allstar.opt"))
    WriteAllstarOptions(im, daodir+"/allstar.opt");
}

void TransformDaoList(const ReducedImage&Ref, const ReducedImage& Im) {

  GtransfoRef refToIm = FindTransfo(Ref, Im);
  double eratio;
  double ratio = PhotoRatio(Ref, Im, eratio, refToIm);
  cout << " Photometric ratio: "
       << ratio << " +/- " << eratio << endl;
  Frame frame(Im.UsablePart());
  string filename = FindDaoList(Ref.Dir() + "calibrated");
  if (filename.empty()) {
    cerr << "TransformDaoList: no daophot like star list found\n";
    return;
  }
  DaoStarList stars;
  ReadDaoList(filename, stars);

  for (DaoStarIterator it = stars.begin(); it != stars.end(); ) {
    DaoStar *star = *it;
    refToIm->TransformStar(*star);
    if (frame.InFrame(*star)) {
      star->flux *= ratio;
      star->eflux *= ratio + eratio * star->flux;
      ++it;
    } else 
      it = stars.erase(it);
  }
}
