#include <iostream>
#include <fstream>

#include <sestar.h>
#include <fitsimage.h>
#include <reducedimage.h>
#include <fileutils.h>

#include "daophot.h"
#include "cdaophot.h"

#define CONVERTSTRING(A) char A##_str[256]; strcpy(A##_str,A.c_str());

static void not_open_error(const string &RoutineName)
{
  cerr << " Daophot::" << RoutineName << "() : fits file not yet open. \n";
}

//***************************** Daophot ****************************
Daophot::Daophot(const ReducedImage &Rim) 
  : open(0), data(0), opt(Rim)
{

  if (!Rim.ActuallyReduced()) 
    {
      cerr << " Daophot::Daophot() : " << Rim.Name() << " was not reduced \n";
      return;
    }

  cout << " Daophot::Daophot(" << Rim.Name() <<")\n";
  
  rootname = CutExtension(Rim.FitsName());
  if (Rim.IsSkySub()) global_sky = Rim.OriginalSkyLevel();

  float rms = Rim.SigmaBack();
  float sky = Rim.BackLevel();
  ap1       = opt[Fwhm].Value();
  lowbad    = sky-opt[SigmaLowDatum].Value()*rms;
  threshold = sky+opt[SigmaThreshold].Value()*rms;
}

Daophot::~Daophot()
{
  if (data) delete [] data;
  CLPIC("DATA");
}

void Daophot::Attach(const string& FitsFileName)
{
  if (!FileExists(FitsFileName))
    {
      cerr << "Daophot::Attach() : Error: " << FitsFileName << " does not exist \n";
      return;
    }
  
  if (open) CLPIC("DATA");

  cout << " Daophot::Attach(" << FitsFileName <<")\n";

  ATTACH(FitsFileName.c_str(), &open);

  if (!open) { not_open_error("Attach"); return; }
  if (!data) data = new float[SIZE.ncol*SIZE.nrow];
}

void Daophot::Find(const string& CooFileName) const
{
  if (!open) {not_open_error("Find"); return;}

  cout << " Daophot::Find() : detect stars" << endl;

  const int maxp   = 3000000;  // enough
  const int maxbox = 13;       // box for convolution
  float *options   = opt.GetArray(GETNOPT());
  float *g         = new float[maxbox*maxbox];
  short int *skip  = new short int[maxbox*maxbox];
  int  *jcyln      = new int[maxp];
  CONVERTSTRING(CooFileName);
  FIND(data,(SIZE.ncol),(SIZE.nrow), data+maxp, jcyln, g, skip, options,CooFileName_str, global_sky);
  
  delete [] jcyln;
  delete [] g;
  delete [] skip;
  delete [] options;
}

void Daophot::Photometry(const string& CooFileName, const string& MagFileName, const string &OptFileName) const
{
  if (!open) {not_open_error("Photometry"); return;}

  cout << " Daophot::Photometry() : aperture photometry" << endl;

  float *sky = new float[GETMAXSKY()];
  int *index = new int[GETMAXSKY()];
  
  CONVERTSTRING(CooFileName);
  CONVERTSTRING(MagFileName);
  CONVERTSTRING(OptFileName);
  
  PHOTSB(data, sky, index, SIZE.ncol, SIZE.nrow, opt[AduHighDatum].Value(), 
	 opt[WatchProgress].Value(), opt[SkyEstimator].Value(), CooFileName_str, 
	 MagFileName_str, OptFileName_str, global_sky);

  delete [] sky;
  delete [] index;
}

void Daophot::Pick(const string& MagFileName, const string& LstFileName, const int Nstar, const float& MagLim) const
{
  const int maxp = 1000; // max number of stars to look
  int *id        = new int[maxp];
  int *index     = new int[maxp];
  CONVERTSTRING(MagFileName);
  CONVERTSTRING(LstFileName);
  PCKPSF(id, data+maxp, data+2*maxp, data+3*maxp, data+4*maxp, 
	 index, maxp, opt[FitRadius].Value(), opt[PsfRadius].Value(), opt[VariablePsf].Value(),
	 MagFileName_str, LstFileName_str, Nstar, MagLim);

  delete [] id;
  delete [] index;
}

void Daophot::Psf(const string& ApFileName, const string& LstFileName, const string& PsfFileName, const string& NeiFileName)
{
  if (!open) {not_open_error("Psf");return;}
  cout << " Daophot::Psf() : Profile : " << opt[AnalyticPsf].Value() 
       << " Variability: " << opt[VariablePsf].Value() << endl;
  float *options = opt.GetArray(GETNOPT());

  /*
    GETPSF(data, &SIZE.ncol, &SIZE.nrow, options,
	 ApFileName.c_str(), LstFileName.c_str(), PsfFileName.c_str(), NeiFileName.c_str(), &global_sky);
  */
  CONVERTSTRING(ApFileName);
  CONVERTSTRING(LstFileName);
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(NeiFileName);
  GETPSF(data,SIZE.ncol,SIZE.nrow, options,ApFileName_str,LstFileName_str,PsfFileName_str,NeiFileName_str,global_sky);
  delete [] options;
}

void Daophot::Peak(const string& MagFileName, const string& PsfFileName, const string& PkFileName) const
{
  if (!open) {not_open_error("Peak");return;}
  cout << " Daophot::Peak() : Fitting a PSF to all stars" << endl;

  CONVERTSTRING(MagFileName);
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(PkFileName);
  DAOPK(data,opt[WatchProgress].Value(),opt[FitRadius].Value(), 
	opt[PercError].Value(),opt[ProfError].Value(), 
	MagFileName_str, PsfFileName_str, PkFileName_str, global_sky);
}

void Daophot::Group(const string& MagFileName, const string& PsfFileName, const string& GrpFileName, const float& CriticalOverlap) const
{

  cout << " Daophot::Group() : clustering stars to fit simultaneously" << endl;
  CONVERTSTRING(MagFileName);
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(GrpFileName);
  GROUP(opt[FitRadius].Value(), opt[PsfRadius].Value(), 
	MagFileName_str, PsfFileName_str, GrpFileName_str, 
	CriticalOverlap, global_sky);
}

void Daophot::Nstar(const string& GrpFileName, const string& PsfFileName, const string& NstFileName) const
{
  if (!open) {not_open_error("Nstar");return;}
  cout << " Daophot::Nstar() : Fitting a PSF to all star groups" << endl;

  CONVERTSTRING(GrpFileName);
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(NstFileName);
  NSTAR(data, SIZE.ncol, SIZE.nrow,
	opt[WatchProgress].Value(), opt[FitRadius].Value(), opt[PercError].Value(), opt[ProfError].Value(), 
	GrpFileName_str, PsfFileName_str, NstFileName_str, global_sky);
}

void Daophot::Substar(const string& PsfFileName, const string& NstFileName, const string& LstFileName, const string& SubPicFileName) const
{
  if (!open) {not_open_error("Substar");return;}
  cout << " Daophot::Substar() : subtract PSF stars " << endl;
  
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(NstFileName);
  CONVERTSTRING(LstFileName);
  CONVERTSTRING(SubPicFileName);
  SUBSTR(data,SIZE.ncol,SIZE.nrow,opt[WatchProgress].Value(), 
	 PsfFileName_str,NstFileName_str, LstFileName_str, SubPicFileName_str);
}

void Daophot::Addstar(const string& PsfFileName, const string& AddFileName, const string& AddPicName, const int InSeed)
{

  Attach(AddPicName);
  cout << " Daophot::AddStar() :  Adding fake stars and noise " << endl;
  
  string outstm = CutExtension(AddFileName);
  // unimportant values to pass anyway
  float rmag[]  = {0., 0.};
  int nstar     = 0;
  int nframe    = 1;
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(AddFileName);
  CONVERTSTRING(AddPicName);
  CONVERTSTRING(outstm);
  ADDSTR(data,SIZE.ncol,SIZE.nrow,  
	 opt[WatchProgress].Value(),PsfFileName_str, AddFileName_str, 
	 AddPicName_str, outstm_str,opt[Gain].Value(), 
	 rmag,nstar,nframe,InSeed);
 
}

void Daophot::Allstar(const string& ApFileName, const string& PsfFileName, const string& AlsFileName, const string& SubPicFileName) const
{
  if (!open) {not_open_error("AllStar"); return;}

  cout << " Daophot::AllStar() : Performing simultaneous PSF fits to stars \n" ;

  float pererr = 0.01 * opt[PercError].Value();
  float proerr = 0.01 * opt[ProfError].Value();
  int iexp     = int(opt[ClipExponent].Value());
  int maxgrp   = int(opt[MaxGroup].Value());
  int center   = 0;
  float *subt  = new float[SIZE.ncol*SIZE.nrow];
  float *sigma = new float[SIZE.ncol*SIZE.nrow];
  if (opt[Recentroid].Value() > 0.5) center = 1;
  int status;
  
  CONVERTSTRING(ApFileName);
  CONVERTSTRING(PsfFileName);
  CONVERTSTRING(AlsFileName);
  CONVERTSTRING(SubPicFileName);
  ALLSTR(data,SIZE.ncol,SIZE.nrow, 
	 subt, sigma, opt[PsfRadius].Value(), 
	 opt[WatchProgress].Value(),opt[ClipRange].Value(), iexp, center, maxgrp, 
	 pererr, proerr, opt[InnerSky].Value(), opt[OutterSky].Value(), status, 
	 PsfFileName_str, ApFileName_str, AlsFileName_str, SubPicFileName_str, 
	 global_sky);

  delete [] subt;
  delete [] sigma;
}

void Daophot::Dump(const Point& Pt, const float& Size) const
{
  const float coords[2] = {Pt.x, Pt.y};

  DUMP(data, &SIZE.ncol, &SIZE.nrow, coords, &Size);  
}

void Daophot::GetSky(float& SkyMean, float& SkyMedian, float& SkyMode, float &SkySigma) const
{

  if (!open) {not_open_error("GetSky"); return;}

  int k;
  int maxs   = min(GETMAXSKY(), SIZE.ncol*SIZE.nrow/3);
  int *index = new int[maxs];
  float *s   = new float[maxs];

  GETSKY(data, s, index, &opt[ReadNoise].Value(), &opt[AduHighDatum].Value(), 
	 &SkyMean, &SkyMedian, &SkyMode, &SkySigma, &k);

  delete [] index;
  delete [] s;
}
