#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>
using namespace std;

#include "image.h"
#include "frame.h"
#include "quali_box.h"
#include "gtransfo.h"
#include "sestar.h"
#include "listmatch.h"

#include "starmatch.h"

void
QualiWriteHeader(char * c, ostream & pr)
{
  pr  << "# xc" << c << " : x of center of box" << endl 
      << "# yc" << c << " : y of center of box" << endl 
      << "# x" << c << " : x of pixel" << endl 
      << "# y" << c << " : y of pixel" << endl  
      << "# fd" << c << " : background of image " << endl  
      << "# sigfd" << c << " :sigma of background of image " << endl  
      << "# v" << c << " : pixel value " << endl
      << "# maxval"<< c << " : max value pixel " << endl;
}

void
QualiWriteHeader(int N, ostream & pr)
{
  for(int i = 0 ; i < N ; i++)
    {
      char c[10] ;
      sprintf(c,"%d",i);
      QualiWriteHeader(c,pr);
    } 
  pr << "# end " << endl ;
}






static void
Test_Qual(int N, const Image ** tabimg, double *X, double * Y,  double * Fond,
	  double * SigFond, int dtaille, ostream & pr)
{
  //  int x,y,xc,yc;
  Pixel *maxval = new Pixel[N];
  for(int k = 0; k< N ; k++)
    {
      maxval[k] = -1e30;
      for (int j = -dtaille; j <= dtaille ; j++)
	{
	  for (int i = -dtaille ; i <= dtaille ; i++) 
	    {
	      int x = (int) (X[k]+ 0.5)  + i ;
	      int y = (int) (Y[k]+ 0.5)  + j ;
	      const Image * img = tabimg[k] ;
              Pixel value = (*img)(x,y);
              if (value > maxval[k]) maxval[k] = value;
	    }
	}
    }
  for (int j = -dtaille; j <= dtaille ; j++)
    {
      for (int i = -dtaille ; i <= dtaille ; i++)
 	{
	  for(int k = 0; k< N ; k++)
	    {
	      int x = (int) (X[k]+ 0.5)  + i ;
	      int y = (int) (Y[k]+ 0.5)  + j ;
	      const Image * img = tabimg[k] ;
	      pr << X[k] << " " << Y[k] << " "
		 << x << " " << y << " " 
		 << Fond[k] << " " << SigFond[k] << " " << (*img)(x,y) << " " << maxval[k] << " ";;
	    }
	  pr << endl ;
	}
    }
  delete [] maxval;

}



static 
int GoodForQualTest(StarMatchList & liste, double saturation, 
			  double prctage, int N_max)
{
  // on enleve les saturees, on prend les xx% plus brillantes
  // en brillance de surface.
  int Ngardees = 0 , Nok = 0;
  for (StarMatchIterator it= liste.begin(); it!= liste.end(); )
    {
      StarMatch starm = *it ;
	
	
      SEStar * pstar1 = (SEStar *) starm.s1;
      SEStar * pstar2 = (SEStar *) starm.s2;
	
      if ( (pstar1->IsOK(saturation)) &&
	   (pstar2->IsOK(saturation)) )
	{
	  it++ ; Nok++ ; ;
	}
      else
	{
	    it = liste.erase(it) ;
	}
    }
  liste.sort(&DecFlux);
  Ngardees = min( (int) (prctage*1.0*Nok) ,N_max)  ;
  cout << Ngardees << " plus brillantes pour alignement photom."<< endl ;
  liste.CutTail(Ngardees);  
  return(Ngardees);

} 




// avec meme cut que pour l'alignement photometrique 
void
Test_Qual_N_same(int N, const FitsHeader & header,  
		 SEStarList & stl1, SEStarList & stl2, 
		 Gtransfo *tf, const Image **pimage, double *Fond,double *SigFond,
		 int dtaille, ostream & pr)
{
  double securite = 20 ;
  Frame frame(header);
  //cut sur starlistes
  SEStarList stl11, stl22;
  stl1.CopyTo(stl11);
  stl2.CopyTo(stl22);
  stl11.CutEdges(frame, securite);
  stl22.CutEdges(frame, securite);


  double saturation = header.KeyVal("SATURLEV");
  if (saturation == 0) 
    {
      saturation = 150000;
      cout << " Pas de SATURLEV dans " << header.FileName() << " on prend saturation = " << saturation << endl;
    }
  BaseStarList *l1 = (BaseStarList*)  &stl11;
  BaseStarList *l2 = (BaseStarList*)  &stl22;
  double dist_max = 2. ;
  StarMatchList  *lmatchphot = ListMatchCollect(*l1, *l2, 
					       tf, dist_max) ;

 
  GoodForQualTest(*lmatchphot,saturation,0.3,50 );
  Test_Qual_N_same(N, *lmatchphot, pimage, Fond,SigFond,dtaille,pr);
  delete lmatchphot;
}



// calcule la moyenne de val_sub / val_ref pour
// n*taille_bin <  d < (n+1) * taille_bi n = 0 1 2
// n*taille_bin <  x < (n+1) * taille_bi n = -2 -1 0 1 2
// n*taille_bin <  y < (n+1) * taille_bi n = -2 -1 0 1 2
// pour les pixels >  nsigma * SigFond
void
Test_Qual_Instantane(StarMatchList & liste, Image & imgref, 
		     Image & imgsub,  double Fondref,double Fondsub,
		     double SigFondref, double taille_bin, double n_sigma,
		     int dtaille, ostream & pr)
{
 
  int Nd[3] = {0,0,0} ;
  int Ndx[5] = {0,0,0,0,0} ;
  int Ndy[5] = {0,0,0,0,0} ;
  double Sd[3] = {0.,0.,0.} ;
  double  Sdx[5] = {0.,0.,0.,0.,0.} ;
  double Sdy[5] = {0.,0.,0.,0.,0.} ;
  for (StarMatchIterator it= liste.begin(); it!= liste.end(); it++)
    {
      StarMatch starm = *it ;

      BaseStar * pstar1 = starm.s1;

 
      double X = pstar1->x  ; // coordonnees SExtractor --> image
      double Y = pstar1->y  ;
      for (int i = -dtaille ; i <= dtaille ; i++)
	{
	  for (int j = -dtaille; j <= dtaille ; j++)
 
	    {
	      int x = (int) (X+ 0.5)  + i ;
	      int y = (int) (Y+ 0.5)  + j ;
	      double vref = imgref(x,y);
	      double vsub = imgsub(x,y);
	      if ( vref - Fondref > n_sigma * SigFondref )
		{
		  double u = (vsub - Fondsub) / (vref - Fondref) ;
		  double d = sqrt( (X-x)*(X-x) + (Y-y)*(Y-y) ) ;
		  double dx = x- X ;
		  double dy = y - Y ;
		  for (int ii = 0 ; ii < 3 ; ii++)
		    {
		      if (( d > ii * taille_bin ) &&
			  ( d < (ii+1) * taille_bin ))
			{
			  Sd[ii] += u ; Nd[ii] += 1 ;
			}
		    }

		  for (int ii = -2 ; ii < 3 ; ii++)
		    {
		      if (( dx > ii * taille_bin ) &&
			  ( dx < (ii+1) * taille_bin ))
			{
			  int n = ii + 2 ;
			  Sdx[n] += u ; Ndx[n] += 1 ;
			}
		      if (( dy > ii * taille_bin ) &&
			  ( dy < (ii+1) * taille_bin ))
			{
			  int n = ii + 2 ;
			  Sdy[n] += u ; Ndy[n] += 1 ;
			}
		    }

		} /* fin test */

	    } /* fin bcle j */	  
	    
	}/* fin bcle i */	 
    }/* fin bcle etoiles */
  pr << setiosflags(ios::fixed);
  pr  << setprecision(1) ;
  for (int ii = 0 ; ii < 3 ; ii++)
    {
      if ( Nd[ii] > 0 )
	{
	  Sd[ii] = Sd[ii]/ Nd[ii] ; pr << Sd[ii] << " " ;
	}
      else
	pr << "---- " << endl ;
    }
  
  for (int ii = 0 ; ii < 5 ; ii++)
    {
      if ( Ndx[ii] > 0 )
	{
	  Sdx[ii] = Sdx[ii]/ Ndx[ii] ;  pr << Sdx[ii] << " " ;
	}
      else
	pr << "---- " << endl ;
      if ( Ndy[ii] > 0 )
	{
	  Sdy[ii] = Sdy[ii]/ Ndy[ii] ;  pr << Sdy[ii] << " " ;
	}
      else
	pr << "---- " << endl ;
    }
}




void
Test_Qual_N_same(int N, StarMatchList & liste, const Image **pimage,
		 double *Fond, double *SigFond, int dtaille, ostream & pr)
{
   QualiWriteHeader(N,pr);
   double *X = new double[N];
   double *Y = new double[N];

  for (StarMatchIterator it= liste.begin(); it!= liste.end(); it++)
    {
      StarMatch starm = *it ;


      BaseStar * pstar1 = starm.s1;

 
      for(int i = 0; i < N; i++)
	{
	  X[i] = pstar1->x  ;
	  Y[i] = pstar1->y  ;
	}
      Test_Qual(N, pimage, X, Y, Fond, SigFond, dtaille, pr);
    }
  delete [] X ;
  delete [] Y ;
return;
}

      








