#include "defs.h"
#include <stdlib.h>
#include <stdio.h>

//#include "fmath.h"
#include "nbrandom.h"

#include "image.h"
#include "fitsimage.h"

#include "dynccd.h"

#define   fsqrt(x)    ((float)(sqrt((double)(x))))




// recommendations si fond et sigma fond donnes:
// pixels, fond et sigma_fond doivent etre en meme unites, 
// ie exprimes avec le meme gain, celui donne









//++
// Class	DynCCD
// Lib		Images++ 
// include	dynccd.h
// 
//	Cette classe permet la specification des parametres 
//	definissant la dynamique du CCD, et doit etre utilise
//      pour le calcul des images de bruit.
//	- TypNoise = 1 :
//	Bruit = Sqrt( (RONoise/Gain)**2 + ValPix/Gain )      
//	- TypNoise = 2 :                                     
//	Bruit = Sqrt( RefSFond**2 + (ValPix-RefFond)/Gain )  
//	- TypNoise = 0 
//	Bruit = 1  (Constant pour tous les pixels) 
//	- TypNoise = 3         
//	Bruit = Sqrt(Abs(ValPix)) 
//
//	Les pixels hors dynamique sont marques (Bruit = 0)
//	( < MinADU ou > MaxADU) pour toutes valeurs de TypNoise.
//	Gain est exprime en electron par ADU, RONoise en electron,
//	Bruit, ValPix et RefFond en ADU.
//--


// Titre	Methodes
//--
//++
// DynCCD(int TypNoise=0, float MinADU=-9.e19, float MaxADU=9.e19, float Gain=1., float RONoise=0., float RefFond=0., float RefSFond=0.);
//	Creation d'un objet DynCCD ("typ" determine la methode de calcul du bruit)
//	|Test verbatim 
//
// void Set(int TypNoise=0, float MinADU=-9.e19, float MaxADU=9.e19, float Gain=1., float RONoise=0., float RefFond=0., float RefSFond=0.);
//	Modification des parametres de la dynamique
// void Print()
// float Noise(float pixel) const
//	Renvoie la valeur du bruit pour "pixel" en ADU.
//--

/* ............................................................ */
/* :::::::::::::  methode de la classe  DynCCD :::::::::::::::: */
/* ............................................................ */

/* --Methode-- */
DynCCD::DynCCD(int typ, float min, float max, float g,
               float ron, float rf, float rfs)
{
if ( (typ >= kConstantNoise) && (typ <= kSqrtADUNoise) )  TypNoise = typ;
else  TypNoise = kConstantNoise;
MinADU = min;  MaxADU = max;
Gain = g;  RONoise = ron;
RefFond = rf;  RefSFond = rfs;
}

/* --Methode-- */
void DynCCD::Set(int typ, float min, float max, float g,
               float ron, float rf, float rfs)
{
if ( (typ >= kConstantNoise) && (typ <= kSqrtADUNoise) )  TypNoise = typ;
MinADU = min;  MaxADU = max;
Gain = g;  RONoise = ron;
RefFond = rf;  RefSFond = rfs;
}

/* --Methode-- */
void DynCCD::SetPhotonNoise(FitsHeader const & header)
{
  Gain = header.KeyVal("TOADGAIN") ;
  RONoise = header.KeyVal("TOADRDON");
  MinADU =  -1.e+30 ;
  MaxADU = 1.e+30;
  TypNoise = kPhotonNoise ; 

}

#ifdef IS_IT_USEFUL

void DynCCD::SetSigFondNoise(FitsHeader const & header, 
			     double Fond, double SigmaFond)
{
  Gain = header.KeyVal("TOADGAIN") ;
  RONoise = header.KeyVal("TOADRDON");
  MinADU =  -1.e+30 ;
  MaxADU = 1.e+30;
  RefFond  =  Fond ; 
  RefSFond =  SigmaFond;
  TypNoise = kSigFondNoise ; 

}
// on rajoute le parametre back_sub au cas ou il ne serait
// pas ds le header (le fond vient juste d'etre soustrait par
// exemple).
void DynCCD::SetSigFondNoise(FitsHeader const & header, bool back_sub)
{
  double sky = 0., sigmasky = 0. ;
  header.Fond_SigFond(sky, sigmasky) ;
  if (back_sub)
    sky = 0. ;
  SetSigFondNoise(header,sky, sigmasky) ; 
}
#endif /* IS_IT_USEFUL */

/* --Methode-- */
void DynCCD::Print()
{
printf("DynCCD: Type= %d Min/MaxADU= %g %g Gain/RoN= %g %g\n",
       TypNoise, MinADU, MaxADU, Gain, RONoise);
if (TypNoise == 2)
  printf("... RefFond= %g  RefSFond= %g \n", RefFond, RefSFond);
return;
}

/* --Methode-- */
float DynCCD::Noise(float pixel) const

/* Cette fonction calcule la valeur du bruit pour pixel     */
/*  Si TypNoise = 1 :                                       */
/*     Bruit = Sqrt( (RONoise/Gain)**2 + fabs(ValPix)/Gain )*/
/*  Si TypNoise = 2 :                                       */
/*     Bruit = Sqrt( RefSFond**2 + fabs(ValPix-RefFond)/Gain )*/
/*  Si TypNoise = 0                                         */
/*     Bruit = 1  (Constant pour tous les pixels)           */
/*  Si TypNoise = 3                                         */
/*     Bruit = Sqrt(Abs(PixelADU))                          */
/*  Les pixels hors dynamique sont marques (Bruit = 0)      */
/*  ( < MinADU ou > MaxADU) pour tout valeur de TypNoise    */

{
  float h,s,ronsq;
  float fond;

  if ( (pixel > MaxADU) || (pixel < MinADU) )   return(0.);
  if ( TypNoise == kConstantNoise)  return(1.);
  if ( TypNoise == kSqrtADUNoise )  return(fsqrt(fabsf(pixel)));

  if ( TypNoise == kSigFondNoise)
    { fond = RefFond;
    ronsq = RefSFond * RefSFond; }
  else
    { fond = 0;
    ronsq = RONoise/Gain;  ronsq *= ronsq; }

  //h = (pixel>fond) ? (float)(pixel-fond) : 0.;
  h = fabs((float)(pixel-fond));
  s = ronsq+h/Gain;
  s = fsqrt(s);
  return(s);
}

/* --------------------------------------------------------------  */
/*      Quelques fonctions pour manipuler des images de bruit      */
/* --------------------------------------------------------------  */

//++
// Module	Images de bruit 
// Lib		Images++ 
// include	dynccd.h
//
//	Ces fonctions permettent le calcul d'image de bruit a partir d'une 
//	image  et d'un objet DynCCD 
//--
//++
// Links 	Voir classes
// DynCCD
// RzImage
// Image<T>
//--
//++
// Titre	Les fonctions 
//--

//++
// Function	Image * NoiseImage(Image const * pim, DynCCD const * dynccd)
// Function	ImgAddNoise(Image&, DynCCD const&)
//	Calcule l'image du bruit et le rajoute a l'image originale
//--


// DH: j'aime pas trop les new qui traine: je supprime

/* Nouvelle-Fonction */
void NoiseImage(Image const & img, Image & imgnoise, DynCCD const & dynccd)

/* non-Creation et Calcul d'une image de bruit a partir de l'image */
/* pim et de dynccd. Voir la methode DynCCD::Noise() pour la   */
/* description du calcul                                       */

{
  float h,s,ronsq;
  float fond, min,max;
  int i, npix;
  float minois, offnois;

  if ((img.Nx() != imgnoise.Nx()) || (img.Ny() != imgnoise.Ny()))
    return;

  const Pixel * pix = img.begin();
  npix = img.Nx()*img.Ny();
  Pixel * pno = imgnoise.begin();

  min = dynccd.MinADU;   
  max = dynccd.MaxADU;


  switch (dynccd.TypNoise)
    {
    case kConstantNoise :
      for(i=0; i<npix; i++)
	{
	  if ( (*pix <= max) && (*pix >= min) )  *pno = 1;
	  else *pno = 0;
	  pix++;   pno++;
	}
      break;

    case kSqrtADUNoise :
      for(i=0; i<npix; i++)
	{
	  if ( (*pix <= max) && (*pix >= min) )  *pno = 1;
	  else *pno = (Pixel) fsqrt(fabsf((float)(*pix)));    
	  pix++;   pno++;
	}
      break;

    case kPhotonNoise :
    case kSigFondNoise :
      if ( dynccd.TypNoise == kSigFondNoise)
	{ 
	  fond = dynccd.RefFond;
	  ronsq = dynccd.RefSFond * dynccd.RefSFond; }
      else
	{ 
	  fond = 0;
	  ronsq = dynccd.RONoise/dynccd.Gain;  ronsq *= ronsq; 
	}

// Calcul de minois / offnois pour obtenir un bruit correct malgre 
// les conversions (float) -> (entier)
    /* ICI, Pixel est un float, donc je ne fais rien (cas r_4)
    switch(img.PixelType())
      {
      case kuint_2:
      case kint_2:
      case kint_4:
        minois = 1.001;
        offnois = 0.5;
        break;
      case kr_4:
      case kr_8:
        minois = 1.e-9;
        offnois = 0.;
        break;
      default:
        minois = 1.e-9;
        offnois = 0.;
        break;
      }*/
      /*cout << "NoiseImage RON ^2 " << ronsq << endl ;
      cout << "NoiseImage fond " << fond << endl ;
      cout << "NoiseImage min " << min << endl ;
      cout << "NoiseImage max " << max << endl ;*/
      for(i=0; i<npix; i++)
	{	      
	  if ( (*pix <= max) && (*pix >= min) ) 
	    {
	      //h = (*pix>fond) ? (float)(*pix)-fond : 0.;
	      h = fabs((float)(*pix-fond));
	      s = ronsq + h/dynccd.Gain;
	      s = fsqrt(s)+offnois;
	      *pno = (s > minois) ? (Pixel)s : (Pixel)minois;
	    } 
	  else 
	    *pno = 0;
	  pix++;   pno++;
	}
      break;
    }
return;
}


/* Nouvelle-Fonction */
void ImgAddNoise(Image& img, DynCCD const& dyn)
{
	Pixel* p = img.begin();
	int nPix = img.Nx() * img.Ny();
	
	for (int i=0; i<nPix; i++, p++)
	  {
	   float uu = dyn.Noise(*p) * NorRand();
	    *p += (Pixel) uu;
	  }
}





