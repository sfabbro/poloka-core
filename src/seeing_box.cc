#include <iomanip>

#include "seeing_box.h"
#include "fileutils.h"
#include "sestar.h"
#include "starlist.h"
#include "fitsimage.h"
#include "frame.h"

// Structure pour le calcul du seeing:
// - selection des etoiles
// - calcul de la moyenne de la fwhm, des moments, des sigmas etc.
DatSeeing::DatSeeing() {
  prlevel = 0 ;
  saturation  = 0  ;
  prctage = 0 ;
  Nmax = 0 ;
}



DatSeeing::DatSeeing(const string &DatacardsFileName, const double Saturlev)
{
  if (!FileExists(DatacardsFileName))
    {
      cerr << " Cannot find " << DatacardsFileName << "to read in datacards, using Default Value " << endl;
      this->Default(Saturlev);
      return;
    }
  DataCards dataCards(DatacardsFileName);
  this->LitDataCard(dataCards);
  this->saturation=Saturlev;
}



void
DatSeeing::Default(const double saturlevel)
{
  prlevel = 0 ;
  saturation  = saturlevel  ;
  prctage = 1. ;
  Nmax = 5000 ;
}


DatSeeing::DatSeeing(const double Saturlev)
{
  Default(Saturlev);
}

void
DatSeeing::Print(ostream & s) const 
{
  s  << "*****DatSeeing******" << endl ;
  s << "prlevel: " <<  prlevel  << endl ;
  s << "saturation: " << saturation<< endl ;
  s << "Calcul sur min de : " << Nmax << " et des " << 100.*prctage << "% + brillantes" << endl ;
}



void 
DatSeeing::LitDataCard(DataCards & data)
{
  prlevel = data.IParam("SEEING_PRLEVEL");
  saturation  =  data.DParam("SEEING_SATURATION"); 
  prctage =  data.DParam("SEEING_PRCTAGE") ;
  Nmax =  data.IParam("SEEING_NMAX") ;
}



// structure de sortie du calul du seeing
SortieSeeing::SortieSeeing(){
   Ngardees=0 ;
   seeing=0 ;
}

void
SortieSeeing::Print(ostream & s) const 
{
  ios::fmtflags  old_flags =  s.flags();
  s  << resetiosflags(ios::scientific) ;
  s  << resetiosflags(ios::fixed)  ;
  s  << setiosflags(ios::fixed)  ;
  s  << " *****SortieSeeing******" << endl ;
  s  << " Number of stars kept for histo filling      : " << Ngardees << endl  ;
  s  << " Number of stars kept for seeing computation : " 
     << Ncalcul << endl ; 
  s  << " ** Results ** : " << endl ;
  s  << " seeing: " << seeing << endl ;
  s << endl ;
  s.flags(old_flags);
}


// Modification du calcul du SEEING Julien le 5 oct 2000

#include "histo2d.h"
int
CalculeSeeingSE(const DatSeeing & dat, SortieSeeing & sortie, 
		const SEStarList & stlse, SEStarList & stlg, 
		const double nsigma)
{

 
  // Just keep stars not saturated, not flagged by sextractor
  int Ngardees = KeepIt(dat.saturation, stlse, stlg);
  int nstars_se = stlg.size();
  int ncutse = int(min(dat.prctage*nstars_se,(float) dat.Nmax)); // prctage% brightest in brithness, at most Nmax.
 stlg.sort(&DecFluxMax);

  if (ncutse < nstars_se) stlg.CutTail(ncutse);  
  cerr << "Keep  " << stlg.size() <<" bright stars to compute the seeing."<< endl;
  if ( dat.prlevel > 3)
    {
      const char *dir = getenv("TOADS_WORK");
      if (dir == NULL )
	{
	  dir = "." ;
	  cerr << "TOADS_WORKS  not defined, use local directory  " << endl;
	}
      string sdir = dir ;
      string nnn = sdir+"/se.seeing.cat" ;
      stlg.write(nnn);
    }
 
  // Histogram construction
  // The grid is adapted to the number of stars
  int Bin2d = max(min(100,int(ncutse/10)),1);
  Histo2d histo(Bin2d,0,10,max(Bin2d/2,1),-6,-2); //HC
  double conv =  (sqrt(2.*log(2.))*2.) ;


  // Fill histogram 
  for (SEStarCIterator it= stlg.begin(); it!=stlg.end(); it++)
    {
      const SEStarRef &pstar= *it;
      double flux = pstar->flux;
      double fluxmax = pstar->Fluxmax();
      double fwhm = pstar->Fwhm();
      double shape = -2.5*log10(flux/fluxmax);
      histo.Fill(fwhm ,shape, 1 );
    } 
  

  //  computation of the mean near the histo max
  double X, Y;
  histo.MaxBin(X, Y);
  double BinX, BinY;
  histo.BinWidth(BinX, BinY);
  
  if ( dat.prlevel > 4)
    {
      cout << "Mode en fwhm: " << X << " (Bin: " << BinX 
	   << "), Mode en mag-mu: " << Y <<  " (Bin: " 
	   << BinY <<  ")" << endl ;
      stlg.write("selected.list");
    }


  double moyenne = X + 0.5 * BinX;
  double sigma = 0.5 * BinX;
  double deltamoy=1;
  int iter = 0;
  double precision = 0.01 ;
  int compteur =0;
  
  do
    {
      compteur =0;
      double moy=0, sig=0;
      for (SEStarIterator it= stlg.begin(); it!=stlg.end(); )
	{
	  SEStar *pstar = *it ;
	  double flux = pstar->flux;
	  double fluxmax = pstar->Fluxmax();
	  double fwhm = pstar->Fwhm();
	  double shape = -2.5*log10(flux/fluxmax);
	  if ( (fabs(fwhm-moyenne) <nsigma * sigma) &&  
	       (Y-BinY < shape) && (shape < Y + BinY))
	    {
	      ++compteur;
	      moy += fwhm;
	      sig += fwhm*fwhm; 
	      ++it;
	    }
	  else {it = stlg.erase(it);}
	}
      moy /= compteur;
      deltamoy = fabs(moy - moyenne) / moy; 
      moyenne = moy;
      sigma = sig/compteur - moyenne * moyenne; 
      if (sigma>0) sigma = sqrt(sigma); else sigma=0;
      if (iter>10) break;
      ++iter;
    }
  while((deltamoy > precision ) || (iter<2));


  cout << " FWHM = " << moyenne << " +/- " << sigma 
       << " (" << iter << " iterations)" << endl;
  double seeing = moyenne / conv;
  sortie.seeing=seeing; 
  sortie.Ngardees= Ngardees ;
  sortie.Ncalcul=  compteur ;

return 1;
}


