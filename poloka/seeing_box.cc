#include <iomanip>

#include <poloka/seeing_box.h>
#include <poloka/fileutils.h>
#include <poloka/sestar.h>
#include <poloka/starlist.h>

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
  saturation  =  data.DParam("SATUR_DEFAULT");
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
  s  << "** Results ** : " << endl ;
  s  << setprecision(5) ;
  s  << "seeing: " << seeing << endl ;
  s  << "seeing (histo): " << seeing_histo << endl ;
  s << endl ;
  s.flags(old_flags);
}


// Modification du calcul du SEEING Julien le 5 oct 2000

int
CalculeSeeingSE(const DatSeeing & dat, SortieSeeing & sortie, 
		const SEStarList & stlse, SEStarList & stlg)
{

 
  // Just keep stars not saturated, not flagged by sextractor and not bad
  int Ngardees = KeepOK(dat.saturation, stlse, stlg);
  int nstars_se = stlg.size();
  int ncutse = int(min(dat.prctage*nstars_se,(float) dat.Nmax)); // prctage% brightest in brithness, at most Nmax.
 stlg.sort(&DecFluxMax);

  if (ncutse < nstars_se) stlg.CutTail(ncutse);  
  cerr << "Keep  " << stlg.size() <<" bright stars to compute the seeing."<< endl;
  if ( dat.prlevel > 3)
    {
      const char *dir = getenv("POLOKA_WORK");
      if (dir == NULL )
	{
	  dir = "." ;
	  cerr << "POLOKA_WORK  not defined, use local directory  " << endl;
	}
      string sdir = dir ;
      string nnn = sdir+"/se.seeing.cat" ;
      stlg.write(nnn);
    }
 
  double mfwhm, bin_fwhm, mshape, bin_shape;
  HistoStarFinder(stlg, mfwhm, bin_fwhm, mshape, bin_shape);
 
  double conv = (sqrt(2.*log(2.))*2.) ;
  sortie.seeing_histo = mfwhm/conv ;
  // DEVRAIT ETRE double moyenne = mfwhm; MAIS COMME IL FAUT RESTER HOMOGENE
  // AVEC CE QUI A ETE FAIT AVANT ....
  double moyenne = mfwhm + 0.5 * bin_fwhm ; 
  double sigma = 0.5 * bin_fwhm ;
  double nsigma=3 ;
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
	  double shape = 0.;
	  if ( fabs(fluxmax) > 1.e-10 )
	    shape = -2.5*log10(flux/fluxmax);
	  if ( (fabs(fwhm-moyenne) <nsigma * sigma) && 
	       (fabs(shape - mshape) < bin_shape) )
	    {
	      ++compteur;
	      moy += fwhm;
	      sig += fwhm*fwhm; 
	      ++it;
	    }
	  else {it = stlg.erase(it);}
	}
      if (compteur > 0 )
	moy /= compteur;
      else
	{
	  moyenne = 0.; sigma=0;  break;
	}
	
      if (fabs(moy) > 1.e-10)
	deltamoy = fabs(moy - moyenne) / moy; 
      moyenne = moy;
      sigma = sig/compteur - moyenne * moyenne; 
      if (sigma> 1.e-10) sigma = sqrt(sigma); 
      else {sigma=0;break;}
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


