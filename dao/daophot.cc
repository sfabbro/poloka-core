#include <iostream>
#include <fstream>

#include <sestar.h>
#include <fitsimage.h>
#include <reducedimage.h>
#include <fileutils.h>

#include "daophot.h"
#include "cdaophot.h"

static void not_open_error(const string &RoutineName)
{
  cerr << " Daophot::" << RoutineName << "() : fits file not yet open. \n";
}

//***************************** Daophot ****************************
// NB: all "hardcoded" numbers are taken to reproduce original daophot programs
// They are usually associated with fortran 77 limitations for dynamical allocation

const int Daophot::MAXSKY = 10000;
const int Daophot::NOPT   = 20;

Daophot::Daophot(const ReducedImage &Rim) 
  : open(0), data(0), psf(0), opt(Rim)
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
  if (psf)  delete psf;

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

  if (!open) {not_open_error("Attach"); return;}
  if (!data) data = new float[SIZE.ncol*SIZE.nrow];
}

void Daophot::Find(const string& CooFileName) const
{
  if (!open) {not_open_error("Find"); return;}

  cout << " Daophot::Find() : detect stars" << endl;

  const int maxp   = 3000000;  // enough
  const int maxbox = 13;       // box for convolution
  const int maxcol = maxp/maxbox;
  float *options   = opt.GetArray(NOPT);
  float *g         = new float[maxbox*maxbox];
  short int *skip  = new short int[maxbox*maxbox];
  int  *jcyln      = new int[maxp];

  FIND(data, data+maxp, jcyln, g, skip, &maxp, 
       &maxbox, &maxcol, &MAXSKY, options, &NOPT, 
       CooFileName.c_str(), &global_sky);

  delete [] jcyln;
  delete [] g;
  delete [] skip;
  delete [] options;
}

void Daophot::Photometry(const string& CooFileName, const string& MagFileName, const string &OptFileName) const
{
  if (!open) {not_open_error("Photometry"); return;}

  cout << " Daophot::Photometry() : aperture photometry" << endl;

  float *sky = new float[MAXSKY];
  int *index = new int[MAXSKY];

  PHOTSB(data, sky, index, &MAXSKY, &SIZE.ncol, &SIZE.nrow, &opt[AduHighDatum].Value(), 
	 &opt[WatchProgress].Value(), &opt[SkyEstimator].Value(), CooFileName.c_str(), 
	 MagFileName.c_str(), OptFileName.c_str(), &global_sky);

  delete [] sky;
  delete [] index;
}

void Daophot::Pick(const string& MagFileName, const string& LstFileName, const int Nstar, const float& MagLim) const
{
  const int maxp = 1000; // max number of stars to look
  int *id        = new int[maxp];
  int *index     = new int[maxp];

  PCKPSF(id, data+maxp, data+2*maxp, data+3*maxp, data+4*maxp, 
	 index, &maxp, &opt[FitRadius].Value(), &opt[PsfRadius].Value(), &opt[VariablePsf].Value(),
	 MagFileName.c_str(), LstFileName.c_str(), &Nstar, &MagLim);

  delete [] id;
  delete [] index;
}

void Daophot::Psf(const string& ApFileName, const string& LstFileName, const string& PsfFileName, const string& NeiFileName)
{
  if (!open) {not_open_error("Psf");return;}
  cout << " Daophot::Psf() : Profile : " << opt[AnalyticPsf].Value() 
       << " Variability: " << opt[VariablePsf].Value() << endl;
  if (!psf) 
    { 
      psf = new DaoPsf(); 
      psf->Allocate(psf->MAXPSF, psf->MAXPAR);
    }

  const float *options = opt.GetArray(NOPT);

  GETPSF(data, &SIZE.ncol, &SIZE.nrow, psf->param, psf->table, options, &NOPT, 
	 ApFileName.c_str(), LstFileName.c_str(), PsfFileName.c_str(), NeiFileName.c_str(), &global_sky);

  delete [] options;
}

void Daophot::Peak(const string& MagFileName, const string& PsfFileName, const string& PkFileName) const
{
  if (!open) {not_open_error("Peak");return;}
  if (!psf)  {cerr << " Daophot::Peak() : Warning PSF not done yet\n"; return; }

  cout << " Daophot::Peak() : Fitting a PSF to all stars" << endl;
  DAOPK(psf->param, &psf->MAXPAR, psf->table, &psf->MAXPSF, &psf->MAXEXP, data,
	&opt[WatchProgress].Value(), &opt[FitRadius].Value(), &opt[PercError].Value(), &opt[ProfError].Value(), 
	MagFileName.c_str(), PsfFileName.c_str(), PkFileName.c_str(), &global_sky);
}

void Daophot::Group(const string& MagFileName, const string& PsfFileName, const string& GrpFileName, const float& CriticalOverlap) const
{
  if (!psf) {cerr << " Daophot::Group() : Warning PSF not done yet\n"; return;}

  const int maxstr = 10000; // max number of stars to group

  GROUP(psf->param,  &psf->MAXPAR, psf->table, &psf->MAXPSF, &psf->MAXEXP, 
	&maxstr, &opt[FitRadius].Value(), &opt[PsfRadius].Value(), 
	MagFileName.c_str(), PsfFileName.c_str(), GrpFileName.c_str(), 
	&CriticalOverlap, &global_sky);
}

void Daophot::Nstar(const string& GrpFileName, const string& PsfFileName, const string& NstFileName) const
{
  if (!open) {not_open_error("Nstar");return;}
  if (!psf) {cerr << " Daophot::Nstar() : Warning PSF not done yet\n"; return;}
  cout << " Daophot::Nstar() : Fitting a PSF to all star groups" << endl;

  NSTAR(psf->param, &psf->MAXPAR, psf->table, &psf->MAXPSF, &psf->MAXEXP, data, &SIZE.ncol, &SIZE.nrow,
	&opt[WatchProgress].Value(), &opt[FitRadius].Value(), &opt[PercError].Value(), &opt[ProfError].Value(), 
	GrpFileName.c_str(), PsfFileName.c_str(), NstFileName.c_str(), &global_sky);
}

void Daophot::Substar(const string& PsfFileName, const string& NstFileName, const string& LstFileName, const string& SubPicFileName) const
{
  if (!open) {not_open_error("Substar");return;}
  if (!psf) {cerr << " Daophot::Substar() : Warning PSF not done yet\n"; return;}
  cout << " Daophot::Substar() : subtract PSF stars " << endl;

  SUBSTR(psf->param, &psf->MAXPAR, psf->table, &psf->MAXPSF, &psf->MAXEXP, data, &SIZE.ncol, &SIZE.nrow, 
	 &opt[WatchProgress].Value(), PsfFileName.c_str(), NstFileName.c_str(), LstFileName.c_str(), SubPicFileName.c_str());
}

/*
void Daophot::Addstar(const string& PsfFileName, const string& AddPicName, const int InSeed)
{

  if (!psf) {cerr << " Daophot::AddStar() : Warning PSF not done yet. Returning. \n"; return;}
  Attach(AddPic);
  cout << " Daophot::AddStar() :  Adding fake stars and noise " << endl;

  ADDSTR(psf->param, &psf->MAXPAR, psf->table, &psf->MAXPSF, &psf->MAXEXP, data, &SIZE.ncol, &SIZE.nrow, 
	 &opt[WatchProgress].Value(), PsfFileName.c_str(), AddPicName.c_str(), ); 
}
*/

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

  ALLSTR(data, &SIZE.ncol, &SIZE.nrow, subt, sigma, &opt[PsfRadius].Value(), 
	 &opt[WatchProgress].Value(),&opt[ClipRange].Value(), &iexp, &center, &maxgrp, 
	 &pererr, &proerr, &opt[InnerSky].Value(), &opt[OutterSky].Value(), &status, 
	 PsfFileName.c_str(), ApFileName.c_str(), AlsFileName.c_str(), SubPicFileName.c_str(), 
	 &global_sky);

  delete [] subt;
  delete [] sigma;
}

void Daophot::Dump(const Point& Pt, const float& Size) const
{
  const float coords[2] = {Pt.x, Pt.y};

  DUMP(data, &SIZE.ncol, &SIZE.nrow, coords, &Size);  
}

void Daophot::PeakFit(SEStar &Star) const
{
  float x      = Star.x - DAOPHOT_TOADS_SHIFT;
  float y      = Star.y - DAOPHOT_TOADS_SHIFT;
  float deltax = (x-1)/psf->xpsf -1 ;
  float deltay = (y-1)/psf->ypsf -1 ;
  float scale  = Star.flux;
  float sky    = Star.Fond();
  float perr   = opt[PercError].Value()*0.01;
  float pkerr  = opt[ProfError].Value()*0.01;
  float radius = min(opt[PsfRadius].Value(), float(((psf->npsf-1.)/2. - 1.)/2.) );
  int lx       = max(1, int(x-radius)+1);
  int ly       = max(1, int(y-radius)+1);
  int nx       = min(SIZE.ncol, int(x-radius)+1) - lx +1;
  int ny       = min(SIZE.nrow, int(y-radius)+1) - ly +1;
  x += -lx+1;
  y += -ly+1;
 
  const int maxbox = 69; // box of the litte vignette around the star
  float *f = new float[maxbox*maxbox];

  int status, niter;
  float errmag, chi, sharp;

  RDARAY("DATA", &lx, &ly, &nx, &ny, &maxbox, f, &status);

  PKFIT(f, &nx, &ny, &maxbox, &x, &y, &scale, &sky, &radius, 
	&lowbad, &opt[AduHighDatum].Value(), &opt[Gain].Value(), &opt[ReadNoise].Value(), &perr, &pkerr, 
	&psf->bright, &psf->type, psf->param, 
	&psf->MAXPAR, &psf->npar, psf->table, 
	&psf->MAXPSF, &psf->MAXEXP, &psf->npsf, &psf->nexp, 
	&psf->nfrac, &deltax, &deltay, &errmag, &chi, 
	&sharp, &niter, &global_sky);

  delete [] f;

  Star.x       = x + DAOPHOT_TOADS_SHIFT;
  Star.y       = y + DAOPHOT_TOADS_SHIFT;
  Star.flux    = scale;
  Star.EFlux() = errmag;
  Star.Chi()   = chi;
  Star.Sharp() = sharp;
  Star.Iter()  = niter;
}

void Daophot::GetSky(float& SkyMean, float& SkyMedian, float& SkyMode, float &SkySigma) const
{

  if (!open) {not_open_error("GetSky"); return;}

  int k;
  int maxs   = min(MAXSKY, SIZE.ncol*SIZE.nrow/3);
  int *index = new int[maxs];
  float *s   = new float[maxs];

  GETSKY(data, s, index, &maxs, &opt[ReadNoise].Value(), &opt[AduHighDatum].Value(), 
	 &SkyMean, &SkyMedian, &SkyMode, &SkySigma, &k);

  delete [] index;
  delete [] s;
}


bool Daophot::LoadPsf(const string& PsfFileName)
{
  if (psf) delete psf;
  psf = new DaoPsf();
  return psf->read(PsfFileName);
}


