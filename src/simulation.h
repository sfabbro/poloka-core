// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMULATION_H
#define SIMULATION_H
#include <iostream>

#include "image.h"
#include "psf.h"
#include "dynccd.h"
#include "gtransfo.h"

enum {GenAleatoire = 0, GenDamier, GenUseList, GenHost};

struct DatGen {
    PSF *mypsf ;

    int type  ;
    int Nstar ;
    double coeffpho ;
    double fond  ;
    double frac_surf  ;
    double delta_x ;
    double delta_y ;
    double sigX  ;
    double sigY  ;
    double rho  ;
    double flux_dyn  ;
    double flux_offset  ; 
  double sig ;

DatGen();
void Print(ostream & pr = cout) const ;

};


void
Star_Generation(Image & img, PSFStarList & stl, 
		DynCCD const & descriptionCCD, DatGen & gen);
void
Star_Generation_Aleatoire(Image & img, PSFStarList & stl, 
		DynCCD const & descriptionCCD, DatGen const & gen);


void
Star_Generation_Damier(Image & img, PSFStarList & stl, 
		       DynCCD const & descriptionCCD, DatGen const & gen);

void
Star_Generation_UseList(Image & img, PSFStarList & stl, 
		DynCCD const & descriptionCCD, DatGen const & gen);
void
Star_Generation_Aleatoire(PSFStarList & stl,  int taille_x, int taille_y, 
			  DatGen const & gen);
void
Star_Generation_Damier(PSFStarList & stl,  int taille_x, 
		       int taille_y, DatGen const & gen);

void
Star_Generation_UseList_Pave(Image & img,  PSFStarList & stl, 
			     DynCCD const & descriptionCCD, double cphot,
			     double sx, double sy, double rho, Gtransfo *tf=NULL);
void
Star_Generation_inHost(SEStarList & stlse, SEStarList & stlh, 
		       PSFStarList & stl,  int taille_x, int taille_y, 
		       DatGen const & gen);


// bruit ajoute au pave puis ajout du pave
int
AddSNImage(Image & img, PSF *psf, double x, double y, double flux,
	   double sigx, double sigy, double rho,  double & fluxmax, 
	   DynCCD & descriptionCCD);


void
GenerationFakeSN_Host(SEStarList & stlse,SEStarList & stlh, 
		      PSFStarList &stl, PSF *mypsf, int taille_x, 
		      int taille_y );
void
GenerationFakeSN_Damier(PSFStarList &stl, PSF *mypsf, int taille_x, 
			int taille_y );
void
GenerationFakeSN(PSFStarList &stl, PSF *mypsf, int taille_x, 
		 int taille_y , int type, SEStarList *pstlse = NULL,
		 SEStarList *pstlh = NULL);






#endif
