#include <fstream>

#include <fileutils.h>

#include "daophot.h"
#include "daophotutils.h"
#include "daophotpsf.h"
#include "daophotio.h"


//********* Yet another star selection routine *********************

static bool CouldBeStar(const SEStar* star, 
			const double& Satur,
			const double& MinFwhm, 
			const double& MaxFwhm, 
			const double& MinCstar, 
			const double& MaxChi,
			const double& CutSharp)
{
  //cout << star->Fwhm() << " " << star->Cstar() << " " << star->Chi() << endl;
  
  return (
	  (star->flux      >  1e-10    && star->flux    <  1e38    ) &&
	  (star->Flag()    == 0                                    ) &&
	  (star->FlagBad() == 0                                    ) &&
	  (star->Fluxmax() + star->Fond() <  Satur                 ) &&
	  (star->Fwhm()    >  MinFwhm  && star->Fwhm()  < MaxFwhm  ) &&
	  (star->A()       <  star->B()*1.5                        ) &&
	  (star->Cstar()   >  MinCstar                             ) &&
	  (star->Chi()     <  MaxChi                               ) &&
	  (star->Sharp()   <  CutSharp && star->Sharp() > -CutSharp)
	  );
}

static bool DecreasingFluxmax(const SEStar *S1, const SEStar *S2)
{
  return (S1->Fluxmax() > S2->Fluxmax());
}

static bool DecreasingFwhm(const SEStar *S1, const SEStar *S2)
{
  return (S1->Fwhm() > S2->Fwhm());
}

static bool DecreasingCstar(const SEStar *S1, const SEStar *S2)
{
  return (S1->Cstar() > S2->Cstar());
}

size_t SelectOKStars(SEStarList& Slist, const double& Satur, 
		   const double& Seeing, const double& MinCstar, 
		   const double& MaxChi, const double& CutSharp)
{
  const double fwhm = Seeing*2.3548;
  for (SEStarIterator it=Slist.begin(); it != Slist.end(); )
    {
      if ((CouldBeStar(*it, Satur, fwhm*0.3, fwhm*1.6, 
		       MinCstar, MaxChi, CutSharp))) ++it;
      else it = Slist.erase(it);
    }
  return Slist.size();
}

size_t SelectIsolatedStars(SEStarList& Slist, const Frame& Borders, 
			 const double& Satur, const double& Seeing, 
			 const double& MinCstar, const double& MaxChi, 
			 const double& CutSharp, const double& RadFwhm=2.)
{
  const double fwhm = Seeing*2.3548;
  const double rad  = RadFwhm*fwhm;

  for (SEStarIterator it=Slist.begin(); it != Slist.end(); )
    {
      if ((CouldBeStar(*it, Satur, fwhm*0.3, fwhm*1.6, 
		       MinCstar, MaxChi, CutSharp))   &&
	  (!Slist.HasCloseNeighbor(**it,  rad)) &&
	  (Borders.MinDistToEdges(**it) > rad))
	++it;
      else it = Slist.erase(it);
    }
  return Slist.size();
}

size_t SelectTwoStars(SEStarList& Slist)
{
  Slist.sort(&DecreasingFlux);
  Slist.CutTail(20);
  Slist.sort(&DecreasingFluxmax);
  Slist.CutTail(10);
  Slist.sort(&DecreasingFwhm);
  Slist.CutTail(5);
  Slist.sort(&DecreasingCstar);
  Slist.CutTail(2);
  for (SEStarIterator it=Slist.begin(); it != Slist.end(); ++it)
    cout << "a star=" << (**it) << endl;
  
  
  return Slist.size();
}

static bool increasing_chi(const SEStar *S1, const SEStar *S2)
{
  return (S1->Chi() < S2->Chi());
}

static size_t filter_chi(SEStarList& Stars, const double& maxchi)
{
  
  SEStarIterator si = Stars.begin();
  // remove chi == 0 and chi>0.1 
  while (si != Stars.end()) 
    {
      if ( ((*si)->Chi() > maxchi) || ((*si)->Chi() < 1e-10)) {
	cout << "erasing star cause chi2 = " << (*si)->Chi() << endl;
	si = Stars.erase(si); 
      }
      else ++si;
    }
  if(Stars.size()==0) {
    cout << "ERROR IN filter_chi" << endl;
    return 0;
  }
  Stars.sort(&increasing_chi);
  double cutchi = 2*Stars.front()->Chi();
  si = Stars.begin();
  while (si != Stars.end() && (*si)->Chi() < cutchi) ++si;
  while (si != Stars.end()) si = Stars.erase(si);
  return Stars.size();
}

bool IteratePsf(Daophot& DaoSession, SEStarList& Stars)
{
  size_t nStars = Stars.size();
  if (!nStars) 
    {
      cerr << " IteratePsf() : Error: trying to build a PSF without stars.\n";
      return false;
    }
  size_t nStarsLeft = nStars;
  const unsigned int maxIter = 10;  
  unsigned int curIter = 0;

  DaoSession.opt[WatchProgress].SetValue(-1);
  const double maxchi(1000);
  do
    {
      nStars = nStarsLeft;
      cout << " IteratePsf() : Iteration #" << curIter << " Nstars #" << nStars << endl;
      Stars.FluxSort();
      DaoSession.WriteSEStarList<DaophotLst>(Stars);
      DaoSession.Psf();
      Stars.ClearList();
      DaoSession.ReadSEStarList<DaophotNei>(Stars);
      nStarsLeft = filter_chi(Stars, maxchi);
      cout << " IteratePsf() : After filtering : Left #" << nStarsLeft << endl;
    }
  while ((curIter++ < maxIter) && (nStars-nStarsLeft) && nStarsLeft);

  return true;
}

void MakeDaoPsf(ReducedImage &Rim, const bool Redo)
{
  if (FileExists(Rim.ImagePsfName()) && !Redo)
    {
      cout << " MakeDaoPsf() : PSF already done for " << Rim.Name() << endl;
      return;
    }

  cout << " MakeDaoPsf() : Building PSF" << endl;
  
  Daophot daoSession(Rim);
  SEStarList psfStars(Rim.CatalogName());
  write_dao<DaophotAp>(Rim, psfStars);
  cout << daoSession.opt << endl;
  daoSession.Attach(Rim.FitsName());
  IteratePsf(daoSession, psfStars);

  if (!UpdateSeeingFromDaoPsf(Rim)) 
    cerr << " MakeDaoPsf() : Warning : DbImage PSF parameters were not updated \n";
}

void MakeExperimentalPsf(ReducedImage &Rim)
{

  cout << " MakeExperimentalPsf() : Building PSF" << endl;

  SEStarList psfStars(Rim.CatalogName());
  psfStars.FluxSort();
  float saturratio = 0.90;
  float maxcstar=0.3; // <=> disable
  SelectIsolatedStars(psfStars, Rim.UsablePart(), Rim.Saturation()*saturratio, Rim.Seeing(), maxcstar , 1., 0.5);  
  Daophot daoSession(Rim);
  daoSession.opt.write(Rim.Dir()+"/daophot.opt");
  daoSession.WriteSEStarList<DaophotAp>(psfStars);
  daoSession.Attach(Rim.FitsName());

  if (!FileExists(Rim.ImageCatalogName(DaophotAls)))
    {
#ifdef DEBUG
      cout << "in MakeExperimentalPsf: write 2 psf stars" << endl;
#endif
      SelectTwoStars(psfStars);
      //SelectOKStars(psfStars, Rim.Saturation()*saturratio, Rim.Seeing(), 0.3, 2., 0.3);
      daoSession.WriteSEStarList<DaophotLst>(psfStars);
      daoSession.opt[AnalyticPsf].SetValue(3);  // Moffat
      daoSession.opt[VariablePsf].SetValue(0);
      daoSession.opt[InnerSky].SetValue(0.);
      daoSession.opt[OutterSky].SetValue(0.);
      daoSession.opt[MaxGroup].SetValue(10.);
      cout << daoSession.opt << endl;
#ifdef DEBUG
      cout << "in MakeExperimentalPsf: daoSession.Psf()" << endl;
#endif
      daoSession.Psf();
#ifdef DEBUG
      cout << "in MakeExperimentalPsf: daoSession.Allstar()" << endl;
#endif
      daoSession.Allstar();
    }

  psfStars.ClearList();

  ReadDaoSex<DaophotAls>(Rim, psfStars, true, true, false);  
  SelectOKStars(psfStars, Rim.Saturation()*saturratio, Rim.Seeing(), 0.3, 2., 0.3);
  psfStars.FluxSort();
  daoSession.WriteSEStarList<DaophotAp>(psfStars);
  psfStars.CutTail(200);

  daoSession.WriteSEStarList<DaophotLst>(psfStars);
  daoSession.opt[VariablePsf].SetValue(-1);
  IteratePsf(daoSession, psfStars);

  // Try all profiles
  daoSession.opt[VariablePsf].SetValue(DaoPsfVariability(psfStars.size()));
  daoSession.opt[AnalyticPsf].SetValue(-6);
  daoSession.Psf();
  UpdateSeeingFromDaoPsf(Rim);
}

void MakeDaoPsfAls(ReducedImage &Rim, const bool Merge, const bool Redo)
{
  bool docat = true, dopsf = true;

  if (FileExists(Rim.ImagePsfName()) && !Redo) 
    {
      cout << " MakeDaoPsfAls() : PSF already done " << endl;
      dopsf = false; 
    }
  
  if (FileExists(Rim.ImageCatalogName(DaophotAls)) && !Redo)
    {
      cout << " MakeDaoPsfAls() : ALLSTAR catalog already done " << endl;
      docat = false;
    }

  if ((dopsf) || (docat))
    {
      cout << " MakeDaoPsfAls() : Starting reduction of " << Rim.Name() << endl;
      SEStarList sex(Rim.CatalogName());
      write_dao<DaophotAp>(Rim, sex);
      Daophot daoSession(Rim);
      daoSession.Attach(Rim.FitsName());
      if (dopsf) 
	{
	  SEStarList psfStars(Rim.CatalogName());
	  IteratePsf(daoSession, psfStars);
	  if (!UpdateSeeingFromDaoPsf(Rim)) 
	    {
	      cerr << " MakeDaoPsfAls() : Warning: PSF parameters not updated \n";
	      cerr << " MakeDaoPsfAls() : Error:  stop daophot processing \n";
	      return;
	    }
	}      

      if (docat) daoSession.Allstar();
    }

  if (!UpdateSeeingFromDaoPsf(Rim)) 
    {
      cerr << " MakeDaoPsfAls() : Warning : DbImage PSF parameters were not updated \n";
      return;
    }

  if (Merge && FileExists(Rim.ImageCatalogName(DaophotAls))) 
    {
      cout << " MakeDaoPsfAls() : Merging ALLSTAR catalog into SExtractor catalog " << endl;   
      MergeDaoSex<DaophotAls>(Rim);
    }

  cout << " MakeDaoPsfAls() : complete " << endl;
}

bool UpdateSeeingFromDaoPsf(ReducedImage &Rim)
{
  DaoPsf psf;
  if (!psf.read(Rim.ImagePsfName()))
    {
      cerr << " UpdateSeeingFromDaoPsf() : Error: problem reading PSF file\n";
      return false;
    }
    
  // compute seeing
  double sigmaX = psf.HwhmX() * 0.8493218;
  double sigmaY = psf.HwhmY() * 0.8493218;
  double thetaXY = psf.ThetaXY();

  // thetaXY is not from the gaussian.
  double r = (1-thetaXY*thetaXY);
  double seeing;
  if ((r > 1.e-10) && (sigmaX>1.e-10) && (sigmaY>1.e-10)) 
    seeing = sqrt(sqrt(r)*sigmaX*sigmaY);
  else seeing =  max(sigmaX,sigmaY); // is this ok?

  double oldseeing = Rim.Seeing();
  //Rim.SetSeeing(seeing, "Seeing updated from DAOPHOT PSF file");
  //Rim.SetPsfShapeParams(sigmaX, sigmaY, thetaXY, "Sigmas and theta updated from DAOPHOT PSF file");

  cout << " UpdateFromDaoPsf() : old seeing = " << oldseeing << " new seeing = " << seeing << endl;  

  return true;
}
