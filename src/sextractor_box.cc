#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>


#include "fileutils.h"
#include "sestar.h"
#include "sextractor_box.h"
#include "fitsimage.h"

/* sextractor header files do not contain any provision for inclusion in C++ sources */

extern "C" {
#include <define.h> /* from sextractor, mandatory for nex one () */
#include <globals.h> /* from sextractor. nice name isn't it ? */
#include <types.h>
	   }


// pour recuperer une liste de SEStar
static SEStarList *Liste_de_SEStar;
static Image *pmask;


/* I added in sextractor a user handle that is called 
for every found object (called user_ana, and set by default to NULL). The data concerning
found objects is stored in 2 different structures
(for the same object), one concerned about the object pixels,
+ a few parameters, and the other with more 
global parameters such as fluxes, errors, astrometry. One has to gather information in both,that is why this user handle has 2 pointer arguments. 
Some description of the variable info can be found in 
types.h, globals.h and param.h in sextractor sources.*/
/* files types.h, analyse.c and scan.c habe been 
   modified in sextractor code. also refine.c, to fixi a bug */


extern "C" {

  // mettre a 1 le pixel (ou les pixels autour de)
  // (x,y) : realisation de masque 
static 
void  Get_Pixel(Image *img, float x, float y)
{
  int d = 1;
  for (int i = -d ; i <=d ; i++)    
    for (int j = -d ; j <=d ; j++)  
      {
	//int xx = rint(x + DECALAGE_SE_IJ +i) ;
	//int yy = rint(y + DECALAGE_SE_IJ +j) ;
	// appele dans le code avant que les coordonnees aient ete decalees
	int xx = (int) (x + i + 0.5 ) ;
	int yy = (int) (y + j+ 0.5 ) ;
	if (( xx > 0 ) && ( yy > 0 ) && 
	    ( xx < img->Nx() )  && ( yy < img->Ny()) )
	  (*img)(xx,yy) = 1;
      }
}


  // convertit les etoiles sextractor en etoile toads
static 
void  Get_SEStar(SEStar *star, objstruct* obj, obj2struct *obj2)
{
  star->x = obj2->sposx + DECALAGE_CAT_SE;
  star->y = obj2->sposy + DECALAGE_CAT_SE;
  // a cause de cet offset, on met DECALAGE_SE_IJ et DECALAGE_IJ_SE a 0.0
  star->flux =   obj2->flux_best;

  star->N()  = obj->number ;
  star->X_Peak() = obj->peakx ;
  star->Y_Peak() = obj->peaky ;
  star->EFlux() = obj2->fluxerr_best  ;
  star->Fluxmax() = obj->peak  ;
  star->Fond() = obj->bkg  ;
  star->Flux_aper() = obj2->flux_auto ;
  star->Eflux_aper() = obj2->fluxerr_auto  ;
  //star->Flux_fixaper() = obj2->flux_aper[0] ;
  //star->Eflux_fixaper() = obj2->fluxerr_aper[0]  ;
  star->Flux_iso() = obj2->flux_iso  ;
  star->Eflux_iso() = obj2->fluxerr_iso  ;
  star->Flux_isocor() = obj2->flux_isocor ;
  star->Eflux_isocor() = obj2->fluxerr_isocor  ;
  star->Fwhm() = obj->fwhm  ;
  star->Kronradius() = obj2->kronfactor ;
  star->Isoarea() = obj->npix  ;
  star->Mxx() = obj->mx2  ;
  star->Myy() = obj->my2  ;
  star->Mxy() = obj->mxy  ;
  star->A() = obj->a ;
  star->B() = obj->b  ;
  star->Gyr_Angle() = obj->theta   ;
  star->Flag() = obj->flag  ;
  star->FlagBad() = obj->imanflag[0] ;
  star->Cstar() = obj2->sprob  ;
  star->Xtrunc() = obj->mxtrunc  ;
  star->Ytrunc() = obj->mytrunc  ;
}



// To recover in SEStar format 1 star of the star list computed by SExtractor, push it in a SEStarList
// 


static 
void  
Get_SEStarList(objstruct* obj, obj2struct *obj2)
{
  // see sestar.h for information on parameters
  SEStar *star = new SEStar; 
  Get_SEStar(star, obj, obj2);
  Liste_de_SEStar->push_back(star); 
}

static 
void  
Get_PixelSat(float x, float y)
{
  Get_Pixel(pmask,   x,y);
}





typedef void (*_SexStarFill)(objstruct* obj, obj2struct *obj2);
  
typedef void (*_SexImgFill)(float x, float y);
  
static 
int 
sex_proc(const AllForSExtractor & data,
	 _SexStarFill StarFill,  _SexImgFill MaskFill, 
          double & Fond, double &SigmaFond)
{

  memset(&prefs,0,sizeof(prefs));
  prefs.pipe_flag = 0;
  prefs.nimage_name = 1; // 1 seule image pour detection & photometrie


  strcpy(prefs.prefs_name, data.SexConfigFileName.c_str());

  char filename[256];
  strcpy(filename,data.FitsFileName.c_str());
  prefs.image_name[0] = filename;

  readprefs(prefs.prefs_name, NULL, NULL, 0);
  
  strcpy(prefs.param_name, data.SexParamName.c_str());
  strcpy(prefs.nnw_name,data.SexNNWName.c_str() );
  strcpy(prefs.filter_name,data.SexFilterName.c_str());

  if (data.back_type_manual)
    {
      prefs.back_type[0]=BACK_ABSOLUTE;
      prefs.nback_type=1;
      prefs.back_val[0] = data.backmean;
      prefs.nback_val=1;
    }
  else
    {
      if (strlen(data.FitsBackName.c_str()) != 0
	  && strlen(data.FitsMiniBackName.c_str()) != 0 )
	{
	  if (prefs.check_name[0] != 0 )
	    free(prefs.check_name[0]);
	  prefs.check_name[0] = (char *)calloc(MAXCHAR,1);      
	  strcpy(prefs.check_name[0],data.FitsBackName.c_str() );
	  prefs.check_type[0] = CHECK_BACKGROUND;
    
	  if (prefs.check_name[1] != 0 )
	    free(prefs.check_name[1]);
	  prefs.check_name[1] = (char *)calloc(MAXCHAR,1);      
	  strcpy(prefs.check_name[1], data.FitsMiniBackName.c_str() );
	  prefs.check_type[1] = CHECK_MINIBACKGROUND;
	  prefs.ncheck_name = 2 ;
	  prefs.ncheck_type = 2 ;
	}
    }

  if (strlen(data.FitsWeightName.c_str()))
    {
      if (prefs.wimage_name[0] != 0 )
	free(prefs.wimage_name[0]);
      prefs.wimage_name[0] = (char *)calloc(MAXCHAR,1);     
      strcpy(prefs.wimage_name[0],data.FitsWeightName.c_str() );
      prefs.nwimage_name = 1 ; 
      prefs.weight_type[0]=WEIGHT_FROMWEIGHTMAP;
      prefs.nweight_type=1;
      prefs.weight_flag = 1;
      prefs.dweight_flag = 1;
      prefs.weightgain_flag=0; // weight map is not a gain map
      prefs.weight_thresh[1] = 0. ;// a mettre a zero
      prefs.nweight_thresh = 1;
      
    }


  if (strlen(data.FitsMaskName.c_str()))
    {
      if (prefs.fimage_name[0] != 0 )
	free(prefs.fimage_name[0]);
      prefs.fimage_name[0] = (char *)calloc(MAXCHAR,1);      
      strcpy(prefs.fimage_name[0],data.FitsMaskName.c_str() );
      prefs.nfimage_name = 1 ;
      prefs.flag_type[0] = FLAG_AND ;
      prefs.nimaisoflag = 1 ;
      prefs.nimaflag= 1 ;
    }
  else
    {
     prefs.nimaflag= 0 ; prefs.nfimage_name = 0;
    }
  
  prefs.user_ana = StarFill; 
  if ( pmask != NULL)
    prefs.user_ana2 = MaskFill;  

  prefs.satur_level = data.saturation ;
  // detection level and photometric level in number of sigmas
  // they wiil be converted by SExtractor in photons,
  // using the sigma it will compute.
  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;
  // to set these levels directly in photons (or adus),
  // using a sigma specified by the user.
   if ( data.sigma_back > 1.e-10 )
     {
       prefs.dthresh[0] = prefs.dthresh[0] * data.sigma_back ;
       prefs.thresh[0] = prefs.thresh[0] * data.sigma_back ;
       prefs.thresh_type[0] = (ThresholdType)(1) ;
     }
  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;

  FitsHeader head(filename, RO);
  prefs.pixel_scale = head.KeyVal("TOADPIXS");
  prefs.gain = head.KeyVal("TOADGAIN");

  makeit();

  Fond = thefield1.backmean;
  SigmaFond= thefield1.backsig ;

  return 1;
} 

} // end of extern "C"


static int file_exists(const char *FileName)
{
if (!FileExists(FileName))
  {
    cerr << "cannot find " << FileName << endl;
    return 0;
  }
return 1;
}


// number the list (mostly for DAOPHOT)
// Two loops to set a number to a star: very optimal!
static void MakeNumbers(SEStarList & stl)
{
  int nmax = 0;
  for (SEStarIterator it=stl.begin(); it!=stl.end(); ++it)
    {
      SEStar *pstar = *it;
      if (pstar->N()> nmax) nmax = pstar->N();
    }

  if (nmax == 0)
    {
      nmax =  1 ;
      for (SEStarIterator it= stl.begin(); it!=stl.end(); ++it)
       {
	 SEStar *pstar = *it;
	 pstar->N() = nmax;
	 nmax++;
       }
    }
}

// make a SEStarList on an image, encapsulage minimal

int
_SEStarListMake(const AllForSExtractor & data,
	       SEStarList &List, double & Fond, 
	       double &SigmaFond,
	       Image * pmask_sat)
{
  Liste_de_SEStar= &List;
  pmask = pmask_sat ;
  if (file_exists(data.FitsFileName.c_str()) &&   file_exists(data.SexConfigFileName.c_str()) &&  file_exists(data.SexParamName.c_str())
       && file_exists(data.SexNNWName.c_str()) && file_exists(data.SexFilterName.c_str()))
    {
     sex_proc(data, Get_SEStarList, Get_PixelSat,
	      Fond, SigmaFond);
     MakeNumbers(List);
     return 1;
    }
  return 0;
}

// encapsulage niveau 1

void
AllForSExtractor::Print()
{
  ForSExtractor::Print();
}

void
ForSExtractor::Print()
{
  cerr << "FitsFileName " << FitsFileName << endl ;
  cerr << "FitsMaskName  " << FitsMaskName << endl ; 
  cerr << "FitsBackName  " << FitsBackName << endl   ;
  cerr << "FitsMiniBackName " << FitsMiniBackName << endl  ;
  cerr << "saturation " << saturation << endl   ;
  cerr << "back_type_manual : " ;
  if (back_type_manual == true)
    cerr << "set " << endl << " a constant value = " << backmean 
	 << " is taken for the background "  << endl   ;
  else
    cerr << "not set " << endl ;
  cerr << "sigma_back : " ;
  if (sigma_back < 0 )
    cerr << "not specified" << endl ;
  else
    cerr << sigma_back  << endl  ;
}

void
AllForSExtractor::FillFromEnvironnement()
{
  const char *sextractor_dir = getenv("TOADSCARDS");
  if (sextractor_dir == NULL )
    {
      sextractor_dir = "." ;
      cerr << "TOADSCARDS not defined, use local directory  " << endl;
    }
  else cerr << " using SExtractor datacards from "  << sextractor_dir << endl;

  string sextr_dir = sextractor_dir ;
  SexConfigFileName = sextr_dir+"/default.sex";
  const char *default1 = getenv("SECARD1");
  if (default1 != NULL )
    {
      SexConfigFileName = sextr_dir+"/default1.sex";
      cout << "Taking datacard: " << SexConfigFileName << endl ;
    }
  SexParamName =sextr_dir+"/default.param" ;
  SexNNWName   = sextr_dir+"/default.nnw";
  SexFilterName   = sextr_dir+"/default.conv" ;
}

  

int
SEStarListMake(const ForSExtractor & shortdata,
	       SEStarList &List, double & Fond, 
	       double &SigmaFond,
	       Image * pmask_sat)
{
  AllForSExtractor data(shortdata);
  data.FillFromEnvironnement();
  int status = _SEStarListMake(data, List, Fond, SigmaFond,
			      pmask_sat);
  return(status);
}

/********************************************/
// Pour pouvoir recuperer le background a partir de la mini carte



extern "C" {
#include <back.h> 
}
/* back_meshx (or width according SExtractor)
is coded in the mini back header as SEXBKGSX,
back_meshy (or height) is coded in the mini back header 
as SEXBKGSY  */
Image *BackFromMiniBack(Image const & minib, int Nx, int Ny, 
			int back_meshx, int back_meshy)
{
  Image *b = new Image(Nx,Ny);
  picstruct f ;
  f.back = minib.begin();
  f.backw =  back_meshx;
  f.backh =back_meshy;
  f.nbackx = minib.Nx();
  f.nbacky =  minib.Ny();
  f.nback = f.nbackx*f.nbacky;
  f.width = Nx;
  f.back_type=BACK_RELATIVE;
  f.dback = makebackspline(&f, f.back);
  float * pback = b->begin();
  float *data = new float[Nx];
  for(int i=0; i < Nx; i++)data[i]=0. ;
  for (int j = 0 ; j < Ny; j++)
    {      
      int y = j;
      f.backline=pback ;
      subbackline(&f, y, data);
      pback += Nx;
    }
  return(b);
}
