#include <math.h>
#include <iostream>
#include <string.h>
#include <string>
#include <stdio.h>

#include "image.h"
#include "psf.h"
#include "starlist.h"
#include "dynccd.h"
#include "gtransfo.h"
#include "nbrandom.h"

#include "simulation.h"


// les coordonnees des sn simulees sont en coordonnees XY


DatGen::DatGen() {
  type = 0 ;
  Nstar = 0 ;
  coeffpho = 0 ;
  fond = 0 ;
  frac_surf = 0 ;
  delta_x = 0;
  delta_y = 0;
  sigX = 0 ;
  sigY = 0 ;
  rho = 0 ;
  flux_dyn = 0 ;
  flux_offset = 0 ;
  sig = 0 ;
}

void
DatGen::Print(ostream & pr) const {
  pr << "*****DatGen******" << endl ;
  pr << "type generation (0: aleatoire, 1:damier) " << type << endl ;
  pr << "fond image " << fond << endl ; 
  pr << "coeff phot " << coeffpho<< endl ; 
  pr << "nbre etoile (aleatoire)" << Nstar << endl ;
  pr << "fraction surface image utilisee (aleatoire)" 
     << frac_surf << endl ;
  pr << "inter-distance x pour damier " << delta_x << endl;
  pr << "inter-distance y pour damier " << delta_y << endl;
  pr << "sigma X " << sigX << endl ;
  pr << "sigma Y " <<sigY << endl ;
  pr << "rho " << rho << endl ;
  pr << "flux = random[0;1]* " << flux_dyn  << " + " << flux_offset << endl ;
  pr << "distance sigma from host core: " << sig << endl ;
}


void
Star_Generation(Image & img, PSFStarList & stl, 
		DynCCD const & descriptionCCD, DatGen & gen)
{
  switch (gen.type)
    {
    case GenAleatoire :
      Star_Generation_Aleatoire(img, stl, descriptionCCD, gen);
      break;

    case GenDamier :
      Star_Generation_Damier(img, stl, descriptionCCD, gen);
        break;
    
    case GenUseList :
      Star_Generation_UseList(img, stl, descriptionCCD, gen);
        break;
    
    default:
      cerr << "type de generation non trouve " << gen.type << endl ;
        break;
    }
    
}
void
Star_Generation_UseList_Pave(Image & img,  PSFStarList & stl, 
			     DynCCD const & descriptionCCD, double cphot,
			     double sx, double sy, double rho, Gtransfo *tf )
{

  int Ncollees = 0 ;
  PSFStarIterator it ;
  DynCCD dCCD = descriptionCCD ;
  dCCD.RONoise = 0 ;
  for(it = stl.begin() ; it != stl.end() ; it++)
    {
      PSFStar *star = *it;
      double xsn = star->x;
      double ysn = star->y;
      if (tf != NULL )
	{
	  double xo, yo ;
	  tf->apply(xsn,ysn, xo, yo);
	  xsn = xo ;
	  ysn = yo ;
	}
      double fluxsn = star->flux*cphot ;
      double fluxmax = 0. ;
      int collee = AddSNImage(img, star->GetPSF(), 
			      xsn, ysn, fluxsn, sx, sy, rho, 
			      fluxmax, dCCD);
      // on ne recupere pas le fluxmax qui doit rester celui 
      // sur l'image "totale"
      Ncollees += collee ;
    }

  cout << Ncollees << " supernovae collees. " << endl ;
}





void
Star_Generation_UseList(Image & img,  PSFStarList & stl, 
			DynCCD const & descriptionCCD, DatGen const & gen)
{

  int taille_x = img.Nx();
  int taille_y = img.Ny(); 

  img += gen.fond ;
  double flux_to_fluxmax = gen.sigX*gen.sigY*sqrt(1-gen.rho*gen.rho)*6.28;
  
  for(PSFStarIterator it = stl.begin() ; it != stl.end() ; it++)
    {
      PSFStar *star = *it;
      star->SigX() = gen.sigX ;
      star->SigY() = gen.sigY ;
      star->Rho() =  gen.rho ;
      //cout << "flux1 : " <<  star->flux  ;
      star->flux = star->flux*gen.coeffpho  ;
      star->Fluxmax() = star->flux/flux_to_fluxmax;
      //cout << "  flux2 : " <<  star->flux << endl ;

      star->UpdParmVect();
      if (( star->x < taille_x) && ( star->y < taille_y ))
	{	  
	  star->ImgAdd(img);
	}
    }

  ImgAddNoise(img,descriptionCCD);
}


void
Star_Generation_Aleatoire(PSFStarList & stl,  int taille_x, int taille_y, 
			  DatGen const & gen)
{ 

  for (int i = 0; i< gen.Nstar ; i++)
    {
      PSFStar *pstar = new PSFStar(gen.mypsf) ;

      // remplissage  etoile
      double r = (rand()*1.0)/(RAND_MAX*1.0) ;
      pstar->x = taille_x* (gen.frac_surf*r + (1.- gen.frac_surf)*0.5) ;
      //cout << r << " " <<(1.- gen.frac_surf)*0.5 << " " << gen.frac_surf*r << " " << pstar->x  << endl ;
      r = (rand()*1.0)/(RAND_MAX*1.0) ;
      pstar->y = taille_y* (gen.frac_surf*r + (1.- gen.frac_surf)*0.5) ;
      r = (rand()*1.0)/(RAND_MAX*1.0) ;
      pstar->flux = r*gen.flux_dyn+ gen.flux_offset ;
      pstar->Fond() = 0. ;
      stl.push_back(pstar);
    }


}

void
Star_Generation_inHost(SEStarList & stlse, SEStarList & stlh, 
		       PSFStarList & stl,  int taille_x, int taille_y, 
		       DatGen const & gen)
{ 
  cerr << " Star_Generation_inHost " << endl ;
  double r2 = (rand()*1.0)/(RAND_MAX*1.0) ;
  cerr << " Star_Generation_inHost " << r2 << endl ;



  double margex = 100 ;
  double margey = 100 ;
  int Nmax = 200 ;
  SEStarIterator it ;
  for(it = stlse.begin() ; it != stlse.end() ; it++)
    {
      SEStar *starh = *it;
      SEStar *starhcopy = new SEStar(*starh);
      if ( ( starh->x > margex) && ( starh->x < taille_x - margex) && 
	   ( starh->y > margey) && ( starh->y < taille_y - margey))
	{
   	  PSFStar *pstar = new PSFStar(gen.mypsf) ;
	  //cerr << "remplissage  etoile" << endl ;
	  //cerr << " 1 " << pstar << endl ;
	  // remplissage  etoile
	  double dx =0., dy =0. ;
	  NormGau(&dx, &dy, 0., 0., gen.sig ,gen.sig  , 0.);
	  //cerr << starh << endl ;
	  double u = starh->x + DECALAGE_SE_XY  + dx  ;
	  pstar->x = u ;
	  //cerr << " 2 " << pstar << endl ;
	  //cerr << " 3 " << pstar << endl ;
	  //cerr << " 4 " << pstar << endl ; 
	  double v = starh->y +  DECALAGE_SE_XY  +  dy;
	  pstar->y = v;
	  //cerr << pstar << endl ;
	  //pstar->y = starh->y + DECALAGE_SE_XY  + dy;


	  double r = (rand()*1.0)/(RAND_MAX*1.0) ;
	  pstar->flux = r*gen.flux_dyn+ gen.flux_offset ;
	  //cerr << pstar << endl ;
	  pstar->Fond() = 0. ;
	  //cerr << pstar << endl ;
	  //cerr << "fin remplissage  etoile" << endl ;
	  stl.push_back(pstar);
	  //cerr << "push1" << endl ;
	  stlh.push_back(starhcopy);
	  //cerr << "push2" << endl ;

	  //exit(0);
	  
	}
    }
  cerr << stl.size() << " " << stlh.size() << " hotes selectionnes " 
       << endl ;
  int N = stl.size() ;
  if ( Nmax < N )
    {
      // on delete N-Nmax objets ds la liste 
      int deltaN = N-Nmax ;
      int Ntot = N ;
      for(int i = 0 ;  i < deltaN ; i++ )
	{
	  
	  double r = (rand()*1.0)/(RAND_MAX*1.0) ;
	  int l = (int) (Ntot * r ) ;
	  if ( l < Ntot)
	    { 
	     PSFStarIterator itp ;
	     itp = stl.begin() ;
	     for (int k = 0 ; k < l ; k++)
	       itp++;
	     stl.erase(itp);
	     SEStarIterator ith ;
	     ith = stlh.begin() ;
	     for (int k = 0 ; k < l ; k++)
	       ith++;
	     stlh.erase(ith);
	     Ntot--;
	    }
	}
    }
  cerr << stl.size() << " " << stlh.size() << " hotes selectionnes " 
       << endl ;
}


void
Star_Generation_Aleatoire(Image & img, PSFStarList & stl, 
			  DynCCD const & descriptionCCD, 
			  DatGen const & gen)
{
  int tx = img.Nx();
  int ty = img.Ny();  
  Star_Generation_Aleatoire(stl,tx,ty, gen);
  Star_Generation_UseList(img,  stl, descriptionCCD, gen);
}


void
Star_Generation_Damier(PSFStarList & stl,  int taille_x, 
		       int taille_y, DatGen const & gen)
{

  int nx = (int) (taille_x*1.0/(gen.delta_x*1.0)) - 2 ;
  int ny = (int) (taille_y*1.0/(gen.delta_y*1.0))  - 2;
  if ( (nx < 0 ) || (ny < 0 ))
    {
      cerr << "Pb taille image versus offset entre etoiles ds Star_Generation_Damier " 
	   << taille_x << " " << gen.delta_x << " "
	   << taille_y << " " << gen.delta_y << endl ;
      return;
    }
  float xx = gen.delta_x ;
  float yy = gen.delta_y ;
  for (int i = 0; i< nx ; i++)
    { 
      yy = gen.delta_y ;
      for (int j = 0; j< ny ; j++)
	  {
	    PSFStar *pstar = new PSFStar(gen.mypsf) ;

	    // remplissage  etoile
	    double r = (rand()*1.0)/(RAND_MAX*1.0) ;
	    r = 2.*r -1. ;
	    pstar->x = xx + r ;
	    r = (rand()*1.0)/(RAND_MAX*1.0) ;
	    r = 2.*r -1. ;
	    pstar->y = yy + r ;

	    r = (rand()*1.0)/(RAND_MAX*1.0) ;
	    pstar->flux = r*gen.flux_dyn+ gen.flux_offset ;
	    pstar->Fond() = 0. ;
	    stl.push_back(pstar);
	    yy += gen.delta_y ;
	  }
	xx += gen.delta_x ;
    }

 
}


void
Star_Generation_Damier(Image & img, PSFStarList & stl, 
		       DynCCD const & descriptionCCD, DatGen const & gen)
{

  int taille_x = img.Nx();
  int taille_y = img.Ny(); 
  Star_Generation_Damier(stl, taille_x, taille_y, gen);
  Star_Generation_UseList(img,  stl, descriptionCCD, gen);
 
}



int
AddSNImage(Image & img, PSF *psf, double xsn, double ysn, double fluxsn,
	   double sigx, double sigy, double rho, double & fluxmax, 
	   DynCCD & descriptionCCD)
{
  //cout << "xsn, ysn : " << xsn << " " << ysn << endl ;
  int taille_x = img.Nx();
  int taille_y = img.Ny(); 
  int dtaille = 20;
  int taille = 2*dtaille + 1;
  PSFStar *pstar = new PSFStar(psf) ; 
  pstar->x = xsn - floor(xsn) + dtaille; 
  pstar->y = ysn - floor(ysn) + dtaille;
  //cout << "x, y : " << pstar->x << " " << pstar->y << endl ;
  pstar->flux = fluxsn;
  pstar->SigX() = sigx ;
  pstar->SigY() = sigy;
  pstar->Rho() =  rho ;
  pstar->Fond() = 0. ;
  int xcoin =  (int) (floor(xsn)) - dtaille;
  int ycoin =  (int) (floor(ysn)) - dtaille;
  //cout << "xcoin,ycoin : " << xcoin << " " << ycoin << endl ;
  pstar->UpdParmVect();

  Image pave(taille,taille);
  pstar->ImgAdd(pave);
  delete pstar ;
  ImgAddNoise(pave,descriptionCCD);
  fluxmax = (double) pave.MaxValue() ;

  if ( ( xcoin > 0 ) && 
       ( xcoin < taille_x - taille ) && 
       ( ycoin > 0 ) && 
       ( ycoin + taille  < taille_y ) )
    {	  
      for(int i = 0; i < taille ; i++)
	for(int j = 0; j < taille ; j++)
	  img(xcoin +i, ycoin+j) += pave(i,j);

      return(1);
    }
  else
    return(0);
}



void
GenerationFakeSN(PSFStarList &stl, PSF *mypsf, int taille_x, 
		 int taille_y , int type, SEStarList *pstlse,
		 SEStarList *pstlh )
{
  if (( pstlse == NULL ) || (pstlh == NULL))
    type =  GenDamier ;
 switch (type)
    {
    case GenDamier :
      GenerationFakeSN_Damier(stl, mypsf,taille_x, taille_y );  
      break;
     case GenHost :
      GenerationFakeSN_Host(*pstlse, *pstlh, stl, mypsf,taille_x, taille_y );  
      break;
    default:
      cerr << "type de generation non trouve " << type << endl ;
      break;
    }
 return ;
}




void
GenerationFakeSN_Damier(PSFStarList &stl, PSF *mypsf, int taille_x, 
			int taille_y )
{
  DatGen datgen;
  datgen.mypsf = mypsf ;
  datgen.type = GenDamier ;
  datgen.coeffpho = 1. ;
  datgen.fond = 0. ;
  datgen.flux_dyn = 20000. ;
  datgen.flux_offset = 500 ;
  datgen.frac_surf = 0.8 ;
  datgen.delta_x = 50 ;
  datgen.delta_y = 50 ;
  Star_Generation_Damier(stl, taille_x, taille_y, datgen);
}

void
GenerationFakeSN_Host(SEStarList & stlse,SEStarList & stlh, 
		      PSFStarList &stl, PSF *mypsf, int taille_x, 
		      int taille_y )
{
  DatGen datgen;
  datgen.mypsf = mypsf ;
  datgen.coeffpho = 1. ;
  datgen.fond = 0. ;
  datgen.flux_dyn = 40000. ;
  datgen.flux_offset =  500 ;
  datgen.frac_surf = 0.8 ;
  datgen.sig = 1. ;
  Star_Generation_inHost(stlse,stlh, stl, taille_x, taille_y, datgen);
  return;
}
