#include <iostream>
#include <fstream>

#include "../dao_stuff/cdaophot.h"
#include "daophot.h"
#include "daophotio.h"
#include "sestar.h"
#include "fitsimage.h"
#include "reducedimage.h"
#include "fileutils.h"

//***************************** AperOpt  ****************************

AperOpt::AperOpt() 
  : Naper(12), InnerSky(0), OutterSky(0)
{
  Radius = new float[Naper];
}

AperOpt::AperOpt(const double &Fwhm)
  : Naper(1), InnerSky(Fwhm*2), OutterSky(Fwhm*7)
{
  Radius = new float[Naper];
  Radius[0] = Fwhm;
}

AperOpt::~AperOpt()
{
  if (Radius) delete [] Radius;
}

void AperOpt::Write(const string ApFile) const 
{
  ofstream wtaper(ApFile.c_str());
  string aperstr ="123456789ABC";

  for (int i=0; i<Naper; ++i) 
    wtaper << "A" << aperstr[i] << "=" << Radius[i] << endl;
  
  wtaper << "IS=" << InnerSky << endl   
	 << "OS=" << OutterSky << endl ;
}

//***************************** Daophot ****************************

static void not_open_error(const string &RoutineName)
{
  cerr << " Daophot::" << RoutineName << "() : fits file not yet open. \n";
}

int DaophotPsfVariability(const int nstars)
{
  if (nstars > NMINVARPSF*2) return 2;
  else if (nstars > NMINVARPSF) return 1;
  return 0;
}

Daophot::Daophot()
  : open(0), opt_name(0), opt_min(0), opt_max(0), opt_default(0), opt(0), data(0), psf(0)
{
  init_opt();
  // WriteOptions();  
  // Option();
}

Daophot::Daophot(const string &FitsName) 
  : open(0), opt_name(0), opt_min(0), opt_max(0), opt_default(0), opt(0), data(0), psf(0)
{
  init_opt();
  WriteOptions();  
  Option();
  Attach(FitsName);
}

Daophot::Daophot(const ReducedImage &Rim) 
  : open(0), opt_name(0), opt_min(0), opt_max(0), opt_default(0), opt(0), data(0), psf(0)
{
  if (!Rim.ActuallyReduced()) 
    {
      cerr << " Daophot::Daophot : " << Rim.Name() << " was not reduced \n";
      return;
    }

  cout << endl << " Daophot : Starting reduction for " << Rim.Name() << endl;

  init_opt();
  SetOptions(Rim);
  rootname = CutExtension(Rim.FitsName());
  if (Rim.IsSkySub()) global_sky = Rim.OriginalSkyLevel();
  Attach(Rim.FitsName());
  if (FileExists(Rim.ImagePsfName())) 
    {      
      psf = new DaoPsf(Rim);
    }

}

Daophot::~Daophot()
{
  if (data) delete [] data;
  if (psf) delete psf;
  if (opt) delete [] opt;
  if (opt_name) delete [] opt_name;
  if (opt_default) delete [] opt_default;
  if (opt_min) delete [] opt_min;
  if (opt_max) delete [] opt_max;

  CLPIC("DATA");
}

void Daophot::WriteOptions(const string FileName, const int First, const int Last) const
{
  ofstream pr(FileName.c_str());
  for (int i=First; i<=Last; ++i)
    {
      string shortname(opt_name[i]);
      RemovePattern(shortname," ");
      pr << shortname[0] << shortname[1] << "=" << opt[i] << endl;
    }
}

void Daophot::init_opt() 
{
  if (!opt_min) opt_min = new float[NOPT];
  if (!opt_max) opt_max = new float[NOPT];
  if (!opt_default) opt_default = new float[NOPT];
  if (!opt_name) opt_name = new const char*[NOPT];
      
  opt_name[0] = " READ NOISE (ADU; 1 frame)"; opt_min[0] = 1e-20; opt_max[0] = 1e20; opt_default[0] = 2;
  opt_name[1] = "    GAIN (e-/ADU; 1 frame)"; opt_min[1] = 1e-20; opt_max[1] = 1e20; opt_default[1] = 1;
  opt_name[2] = "LOW GOOD DATUM (in sigmas)"; opt_min[2] = 0.   ; opt_max[2] = 1e20; opt_default[2] = 7;
  opt_name[3] = "  HIGH GOOD DATUM (in ADU)"; opt_min[3] = 0.   ; opt_max[3] = 1e20; opt_default[3] = 150000;
  opt_name[4] = "            FWHM OF OBJECT"; opt_min[4] = 0.2  ; opt_max[4] = 20. ; opt_default[4] = 2.5;
  opt_name[5] = "     THRESHOLD (in sigmas)"; opt_min[5] = 0.   ; opt_max[5] = 1e20; opt_default[5] = 4;
  opt_name[6] = " LS (LOW SHARPNESS CUTOFF)"; opt_min[6] = 0.   ; opt_max[6] = 1.  ; opt_default[6] = 0.2;
  opt_name[7] = "HS (HIGH SHARPNESS CUTOFF)"; opt_min[7] = 0.6  ; opt_max[7] = 2.  ; opt_default[7] = 1;
  opt_name[8] = " LR (LOW ROUNDNESS CUTOFF)"; opt_min[8] = -2   ; opt_max[8] = 0.  ; opt_default[8] = -1;
  opt_name[9] = "HR (HIGH ROUNDNESS CUTOFF)"; opt_min[9] = 0.   ; opt_max[9] = 2.  ; opt_default[9] = 1;
  opt_name[10]= "            WATCH PROGRESS"; opt_min[10]= -2   ; opt_max[10]= 2   ; opt_default[10]= -2;
  opt_name[11]= "            FITTING RADIUS"; opt_min[11]= 1    ; opt_max[11]= 30  ; opt_default[11]= 2.5;
  opt_name[12]= "                PSF RADIUS"; opt_min[12]= 5    ; opt_max[12]= 51  ; opt_default[12]= 11;
  opt_name[13]= "              VARIABLE PSF"; opt_min[13]= -1.5 ; opt_max[13]= 3.5 ; opt_default[13]= -1;
  opt_name[14]= "             SKY ESTIMATOR"; opt_min[14]= -0.5 ; opt_max[14]= 3.5 ; opt_default[14]= 0;
  opt_name[15]= "        ANALYTIC MODEL PSF"; opt_min[15]= -6.5 ; opt_max[15]= 6.5 ; opt_default[15]= 3;
  opt_name[16]= " EXTRA PSF CLEANING PASSES"; opt_min[16]= 0.   ; opt_max[16]= 9.5 ; opt_default[16]= 5;
  opt_name[17]= "   USE SATURATED PSF STARS"; opt_min[17]= 0.   ; opt_max[17]= 1.  ; opt_default[17]= 0;
  opt_name[18]= "      PERCENT ERROR (in %)"; opt_min[18]= 0.   ; opt_max[18]= 100 ; opt_default[18]= 0.75;
  opt_name[19]= "      PROFILE ERROR (in %)"; opt_min[19]= 0.   ; opt_max[19]= 100 ; opt_default[19]= 3.0; // 0-19 daophot
  opt_name[20]= "            FITTING RADIUS"; opt_min[20]= 1.   ; opt_max[20]= 30  ; opt_default[20]= 2.5;
  opt_name[21]= "    CE (CLIPPING EXPONENT)"; opt_min[21]= 0.   ; opt_max[21]= 8   ; opt_default[21]= 6.;
  opt_name[22]= "     REDETERMINE CENTROIDS"; opt_min[22]= 0.   ; opt_max[22]= 1   ; opt_default[22]= 1.;
  opt_name[23]= "       CR (CLIPPING RANGE)"; opt_min[23]= 0.   ; opt_max[23]= 10  ; opt_default[23]= 2.5;
  opt_name[24]= "            WATCH PROGRESS"; opt_min[24]= -2   ; opt_max[24]= 2   ; opt_default[24]= -2;
  opt_name[25]= "        MAXIMUM GROUP SIZE"; opt_min[25]= 1.   ; opt_max[25]= 100 ; opt_default[25]= 70;
  opt_name[26]= "      PERCENT ERROR (in %)"; opt_min[26]= 0.   ; opt_max[26]= 100 ; opt_default[26]= 0.75;
  opt_name[27]= "      PROFILE ERROR (in %)"; opt_min[27]= 0.   ; opt_max[27]= 100 ; opt_default[27]= 3;
  opt_name[28]= "     IS (INNER SKY RADIUS)"; opt_min[28]= 0.   ; opt_max[28]= 35  ; opt_default[28]= 0;
  opt_name[29]= "     OS (OUTER SKY RADIUS)"; opt_min[29]= 0.   ; opt_max[29]= 100 ; opt_default[29]= 0;  // 20-29 allstar
  
  if (!opt) opt = new float[NOPT];
  for (int i=0; i<NOPT; ++i) opt[i] = opt_default[i]; 
}

void Daophot::SetOptions(const ReducedImage &Rim)
{
  cout << " Daophot : Setting option parameters from image\n";
  opt[0] = Rim.ReadoutNoise();
  opt[1] = Rim.Gain();
  opt[3] = Rim.Saturation();
  opt[4] = Rim.Seeing() * 2.3548;
  opt[11] = opt[4];   // fit radius
  opt[12] = 4 * opt[4];   // psf radius
  opt[18] = Rim.FlatFieldNoise();
  opt[19] = Rim.ProfileError();
  opt[20] = opt[11];
  // opt[28] = 1.5 * opt[4]; // inner sky radius for allstar
  // opt[29] = 7.0 * opt[4]; // outer sky radius for alltstar
  lowbad = Rim.BackLevel() - opt[2]* Rim.SigmaBack();
  threshold = opt[5]* Rim.SigmaBack();
  check_opt();

}

void Daophot::check_opt() 
{
  // check opt with min and max
  for (int i=0; i<NOPT; ++i)
    {
      if ((opt[i] < opt_min[i]) || opt[i] > opt_max[i]) 
	{
	  cout << " Daophot::opt_check() Replacing \n" << opt_name[i] << "=" << opt[i] 
	       << " by default value = " << opt_default[i] << endl;
	  opt[i] = opt_default[i];
	}
    }

  cout << endl;  

  // print options names and values on scree, daophot style
  cout  << setiosflags(ios::fixed);
  for (int i=0; i<NOPT-1; i+=2) 
    {
      cout << ' ' << opt_name[i] << " =" << setw(9) << setprecision(2) << opt[i]
	   << "    " << opt_name[i+1] << " =" << setw(9) << setprecision(2) << opt[i+1] 
	   << endl;
    }
  cout << resetiosflags(ios::fixed);
  cout << endl;
}

void Daophot::Option(const string OptFile) const
{
  int istat;
  cout << " Daophot : Reading option file " << OptFile << endl;
  const int nchar = 26;
  char *names = new char[NOPT*nchar];

  // inverse array for fortran
  for (int i=0; i<NOPT; ++i) 
    for (int j=0; j<nchar; ++j) 
      {
	names[i*nchar+j] = opt_name[i][j];
      }

  OPTION(OptFile.c_str(), &NOPT, names, opt, opt_min, opt_max, "OPT>", &istat);

  delete [] names;

  // we should not need this stuff very often
  delete [] opt_min;
  delete [] opt_max;
  delete [] opt_default;
  delete [] opt_name;
}

void Daophot::Attach(const string FitsName)
{
  if (!FileExists(FitsName))
    {
      cerr << "Daophot::Daophot(FitsName) : " << FitsName << " does not exist \n";
      return;
    }
  
  if (open) CLPIC("DATA");

  cout << " Daophot : Attaching fits file " << FitsName << endl;

  ATTACH(FitsName.c_str(), &open);

  if (!open) {not_open_error("Attach"); return;}

  if (!data) data = new float[SIZE.ncol*SIZE.nrow];
}


void Daophot::Photometry() const
{
  if (!open) {not_open_error("Photometry"); return;}

  cout << " Daophot : Performing aperture photometry" << endl;

  float *sky = new float[MAXSKY];
  int *index = new int[MAXSKY];

  PHOTSB(data, sky, index, &MAXSKY, &SIZE.ncol, &SIZE.nrow, opt+3, opt+10, opt+14, 
	 &global_sky);

  delete [] sky;
  delete [] index;
}

void Daophot::Pick() const
{
  cerr << " !!! Pick Not yet implemented !!!" << endl;
  //  PCKPSF(id, x, y, m, s, index, maxn, opt+11, opt+12, opt+13);
}

void Daophot::Psf(const int Variability, const bool Manual)
{
  if (!open) {not_open_error("Psf");return;}

  cout << " Daophot : Building a PSF " << endl;

  if (!psf) psf = new DaoPsf;
  psf->Allocate();

  opt[13] = Variability;
  float oldwatch = opt[10];
  if (Manual) opt[10] = 2;

  GETPSF(data, &SIZE.ncol, &SIZE.nrow, psf->param, psf->table, opt, &NOPT, &global_sky);

  opt[10] = oldwatch;

}

void Daophot::Peak() const
{
  if (!open) {not_open_error("Peak");return;}
  if (!psf) {cerr << " Daophot::Peak PSF not done yet\n";return;}

  cout << " Daophot : Fitting a PSF " << endl;

  DAOPK(psf->param, &MAXPAR, psf->table, &MAXPSF, &MAXEXP, data, opt+10, opt+11, opt+18, opt+19, &global_sky);
}

void Daophot::PeakFit(SEStar &Star) const
{

  float x = Star.x - DAOPHOT_TOADS_SHIFT;
  float y = Star.y - DAOPHOT_TOADS_SHIFT;
  float deltax = (x-1)/psf->xpsf -1 ;
  float deltay = (y-1)/psf->ypsf -1 ;
  float scale = Star.flux;
  float sky = Star.Fond();
  float errmag, chi, sharp;
  int niter;
  float perr = opt[18]*0.01;
  float pkerr = opt[19]*0.01;
  float radius = min(opt[11], float(((psf->npsf-1.)/2. - 1.)/2.) );
  const int maxbox = 69;
  int lx = max(1, int(x-radius)+1);
  int ly = max(1, int(y-radius)+1);
  int nx = min(SIZE.ncol, int(x-radius)+1) - lx +1;
  int ny = min(SIZE.nrow, int(y-radius)+1) - ly +1;
  x += -lx+1;
  y += -ly+1;
  int status;
  
  float *f = new float[maxbox*maxbox];

  RDARAY("DATA", &lx, &ly, &nx, &ny, &maxbox, f, &status);
  PKFIT(f, &nx, &ny, &maxbox, &x, &y, &scale, &sky, &radius, 
	&lowbad, opt+3, opt+1, opt, &perr, &pkerr, 
	&psf->bright, &psf->type, psf->param, 
	&MAXPAR, &psf->npar, psf->table, 
	&MAXPSF, &MAXEXP, &psf->npsf, &psf->nexp, 
	&psf->nfrac, &deltax, &deltay, &errmag, &chi, 
	&sharp, &niter, &global_sky);
  delete [] f;

  Star.x = x + DAOPHOT_TOADS_SHIFT;
  Star.y = y + DAOPHOT_TOADS_SHIFT;
  Star.flux = scale;
  Star.EFlux() = errmag;
  Star.Chi() = chi;
  Star.Sharp() = sharp;
  Star.Iter() = niter;
}


void Daophot::IterPsf(const ReducedImage &Rim, const bool Manual)
{
  PsfStars psfStars(Rim);
  PsfConditions cond(Rim);
  int nstars = psfStars.FilterCond(cond);
  if (nstars == 0) 
    {
      cerr << " Daophot::IterPsf  no PSF stars found ! "<< endl;
      return;
    }
  psfStars.CutTail(30);
  nstars = psfStars.size();
  const int maxiter = 10;
  int iter = 0;
  size_t nstarschange = 1;
  float oldfit = opt[11];
  opt[11] = opt[4]*1.5;   // fit radius 1.5*fwhm because pure analytical
  do
    {
      write_dao<DaophotLst>(Rim, psfStars);
      cout << " Daophot::IterPsf  Iteration #" << iter 
	   << " Nstars #" << nstars << endl;
      Psf(-1, Manual); // analytical
      size_t nkept = psfStars.FilterNei();
      psfStars.sort(&DecreasingFlux);
      cout << " Daophot::IterPsf  FilterNei nkept #" << nkept << endl;
      nstarschange = nstars-nkept;
      nstars = nkept;
      iter++;
    }
  while (nstars && iter<maxiter && nstarschange);
  Psf(DaophotPsfVariability(nstars), Manual);
  opt[11] = oldfit;   // fit radius back to normal
  if (!psf) psf = new DaoPsf();
  psf->read(Rim.ImagePsfName());
}

void Daophot::PrecisePsf(const ReducedImage &Rim)
{
  SEStarList psfStars;
  read_dao<DaophotLst>(Rim.ImageCatalogName(DaophotLst), psfStars);
  size_t nstars = psfStars.size();
  if (nstars == 0) 
    {
      cerr << " Daophot::PrecisePsf  no PSF stars found ! "<< endl;
      return;
    }
  opt[15] = -6;
  Psf(DaophotPsfVariability(nstars));

  if (!psf) psf = new DaoPsf();
  psf->read(Rim.ImagePsfName());
}

void Daophot::AllStar(const string& FileWithStarsToFit) const
{
  if (!open) {not_open_error("AllStar"); return;}

  cout << " Daophot : Performing simultaneous PSF fits to stars for " 
       << FileWithStarsToFit << endl;

  strcpy(FILNAM.MAGFIL, FileWithStarsToFit.c_str());
  size_t slen = FileWithStarsToFit.length();
  for(size_t i=slen;i<NCHARFILE;i++) FILNAM.MAGFIL[i]=' ';  

  
  float pererr = 0.01 * opt[26];
  float proerr = 0.01 * opt[27];
  int iexp = int(opt[20]);
  int maxgrp = int(opt[25]);

  int center = 0;
  if (opt[22] > 0.5) center = 1;

  int istat;

  float *subt = new float[SIZE.ncol*SIZE.nrow];
  float *sigma = new float[SIZE.ncol*SIZE.nrow];

  ALLSTR(data, &SIZE.ncol, &SIZE.nrow, subt, sigma, opt+20, 
	 opt+24, opt+23, &iexp, &center, &maxgrp, 
	 &pererr, &proerr, opt+28, opt+29, &istat, &global_sky);

  delete [] subt;
  delete [] sigma;
}


void Daophot::AddStars()
{

  if (!psf) {cerr << " Daophot::AddStars PSF not done yet\n";return;}

  Attach(rootname+"_fake.fits");
  cout << " Daophot : Adding fake stars and noise " << endl;
  ADDSTR(psf->param, &MAXPAR, psf->table, &MAXPSF, &MAXEXP, data, &SIZE.ncol, &SIZE.nrow, opt+10); 
}

void Daophot::Sky(float& SkyMean, float& SkyMedian, float& SkyMode, float &SkySigma) const
{

  if (!open) {not_open_error("Sky"); return;}

  int k, maxs = min(MAXSKY, SIZE.ncol*SIZE.nrow/3);
  int *index = new int[maxs];
  float *s = new float[maxs];

  GETSKY(data, s, index, &maxs, opt, opt+4, 
	 &SkyMean, &SkyMedian, &SkyMode, &SkySigma, &k);

  delete [] index;
  delete [] s;
}
