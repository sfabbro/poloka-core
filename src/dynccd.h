/*  Caracteristique du CCD      Reza  03/95  -*- C++ -*- */
/* Modif pour toads Delphine 050299 */


/*DOC Structures and Methods for computing or applying noise on image */

/* DOC used when: 1. fitting one image pixels
                  2. adding noise to a simulated image */ 


#ifndef DYNCCD_H_SEEN
#define DYNCCD_H_SEEN

#include "exceptions.h"

#include "fitsimage.h"


enum {kConstantNoise = 0, kPhotonNoise, kSigFondNoise, kSqrtADUNoise};

class DynCCD
{
  public :

  int   TypNoise;        // Type de calcul du bruit 
  float MinADU,MaxADU;   // Valeur min et max des pixels 
  float Gain;            // Gain en electron/ADU 
  float RONoise;         // Bruit de lecture en electron  
  float RefFond;         // Fond de ciel de ref (TypNoise=2) 
  float RefSFond;        // Sigma du fond  (TypNoise=2) */

/*    ===>  Les methodes   */
  DynCCD(int typ=0, float min=-9.e19, float max=9.e19,
         float g=1., float ron=0., float rf=0., float rfs=0.);

  void Set(int typ=0, float min=-9.e19, float max=9.e19,
         float g=1., float ron=0., float rf=0., float rfs=0.);

  void SetPhotonNoise(FitsHeader const & header);

#ifdef IS_IT_USEFUL
  void SetSigFondNoise(FitsHeader const & header, double Fond, double SigmaFond);
  void SetSigFondNoise(FitsHeader const & header, bool fond_sub = false);
#endif /* IS_IT_USEFUL */
  void Print();
  float Noise(float pixel) const;
};

// Quelques fonctions pour manipuler des images de bruit


/* DOCF This function computes the noise on each pixel
   Cette fonction calcule la valeur du bruit pour pixel    
   if TypNoise = 1 :                                       
   Noise = Sqrt( (RONoise/Gain)**2 + abs(ValPix)/Gain )   
   If TypNoise = 2 :                                       
   Noise = Sqrt( RefSFond**2 + abs(ValPix-RefFond)/Gain ) 
   If TypNoise = 0                                         
   Noise = 1  (Constant pour tous les pixels)           
   If TypNoise = 3                                         
   Noise = Sqrt(Abs(PixelADU))                          
   The noise for pixels out of dynamic 
   ( < MinADU ou > MaxADU) is marked out   = 0      
   for every TypNoise value*/ 

void NoiseImage(Image const & img,Image  & imgnoise, 
		  DynCCD const & dynccd);

/* DOCF This function add the noise on each pixel */

void ImgAddNoise(Image&, DynCCD const&);


#endif
