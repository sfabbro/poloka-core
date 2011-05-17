#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "fileutils.h"
#include "sestar.h"
#include "sextractor_box.h"
#include "fitsimage.h"

// used when recovering positions from sextractor
#define DECALAGE_CAT_SE -1.0

/* sextractor header files do not contain any provision for inclusion in C++ sources */

extern "C" {
#ifdef OFF_T
#define OFF_T_CFITSIO OFF_T
#undef OFF_T
#endif
#include <sex/define.h> /* from sextractor, mandatory for next one () */
#include <sex/globals.h> /* from sextractor. nice name isn't it ? */
#include <sex/types.h> 
#include <sex/prefs.h> // a ajouter pour la version v244
#include <sex/back.h> 
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
/* files types.h, analyse.c and scan.c have been 
   modified in sextractor code. also refine.c, to fix a bug */


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

  star->vx = obj->poserr_mx2;
  star->vy = obj->poserr_my2;
  star->vxy = obj->poserr_mxy;

  star->N()  = obj->number ;
  star->X_Peak() = obj->peakx + DECALAGE_CAT_SE ;
  star->Y_Peak() = obj->peaky + DECALAGE_CAT_SE ;
  //star->X_Peak() = obj->dbkg ; // pour tester que c'est tjrs =0 !
  star->EFlux() = obj2->fluxerr_best  ;
  star->Fluxmax() = obj->peak  ;
  star->Fond() = obj->bkg  ;
  star->Flux_auto() = obj2->flux_auto ;
  star->Eflux_auto() = obj2->fluxerr_auto  ;

  // doit etre defini dans le fichier de param !
  if (obj2->flux_aper)
    {
      star->Flux_circ_aper() = obj2->flux_aper[0] ;
      star->Eflux_circ_aper() = obj2->fluxerr_aper[0]  ;
    }


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
  // imanflag denotes the number of flagged pixels that this object
  // used. we convert it to a binary flag 0 or 1.
  if (obj->imanflag[0] != 0)  star->FlagBad() = 1;
  star->Cstar() = obj2->sprob  ;
  star->Xtrunc() = obj->mxtrunc  + DECALAGE_CAT_SE;
  star->Ytrunc() = obj->mytrunc  + DECALAGE_CAT_SE;
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

 



#include <vector>

static 
int 
sex_proc(const AllForSExtractor & data,
	 _SexStarFill StarFill,  _SexImgFill MaskFill, 
          double & Fond, double &SigmaFond)
{

  std::vector<void *> toFree;
  memset(&prefs,0,sizeof(prefs));
  prefs.pipe_flag = 0;
  prefs.nimage_name = 1; // 1 seule image pour detection & photometrie


  strcpy(prefs.prefs_name, data.SexConfigFileName.c_str());

  char filename[512];
  strcpy(filename,data.FitsFileName.c_str());
  prefs.image_name[0] = filename;

  // tell SExtractor to read its datacards.
  readprefs(prefs.prefs_name, NULL, NULL, 0);
  

  // and now fill the remainder.
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
      int nchecks = 0;
      /* just a simple comment: SExtractor actually free's the prefs.check_name.
	 so if you free them a second time, a (usualy delayed) crash should happen */
      size_t namesize;
      if ((namesize=strlen(data.FitsBackName.c_str())) != 0)
	{
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks],data.FitsBackName.c_str() );
	  prefs.check_type[nchecks] = CHECK_BACKGROUND;	  
	  nchecks++;
	}
      if ((namesize=strlen(data.FitsMiniBackName.c_str())) != 0 )
	{  
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks], data.FitsMiniBackName.c_str() );
	  prefs.check_type[nchecks] = CHECK_MINIBACKGROUND;
	  nchecks++;
	}
      if ((namesize = strlen(data.FitsSegmentationName.c_str())) != 0 )
	{  
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks], data.FitsSegmentationName.c_str() );
	  cout << "Segmentation image requested " << data.FitsSegmentationName << endl;
	  prefs.check_type[nchecks] = CHECK_SEGMENTATION;
	  nchecks++;
	}
      prefs.ncheck_name = nchecks ;
      prefs.ncheck_type = nchecks ;
    }

  if (strlen(data.FitsWeightName.c_str()))
    {
      prefs.wimage_name[0] = (char *)calloc(MAXCHAR,1);
      toFree.push_back(prefs.wimage_name[0]);
      strcpy(prefs.wimage_name[0],data.FitsWeightName.c_str());
      prefs.nwimage_name = 1 ; 
      prefs.weight_type[0]=WEIGHT_FROMWEIGHTMAP; 
      // if no weight map is to be taken into account
      // prefs.weight_type[0]=WEIGHT_NONE ;
      prefs.nweight_type=1;
      prefs.weight_flag = 1;
      prefs.dweight_flag = 1;
      prefs.weightgain_flag=0; // weight map is not a gain map
      prefs.weight_thresh[1] = 0. ;// a mettre a zero
      prefs.nweight_thresh = 1;      
    }


  if (strlen(data.FitsMaskName.c_str()))
    {
      prefs.fimage_name[0] = (char *)calloc(MAXCHAR,1);
      toFree.push_back(prefs.fimage_name[0]);
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
    {
      prefs.user_ana2 = MaskFill;  
      cout << " We will create satur mask from sextractor pixel lists"<< endl;
    }

  
  // detection level and photometric level in number of sigmas
  // they will be converted by SExtractor in photons,
  // using the  SExtractor computed sigma 
  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;
  // to set these levels directly in photons (or adus),
  // using a sigma specified by the user.
   if ( data.sigma_back > 1.e-10 )
     {
       prefs.dthresh[0] = prefs.dthresh[0] * data.sigma_back ;
       prefs.thresh[0] = prefs.thresh[0] * data.sigma_back ;
       // check if it is needed...
       prefs.thresh_type[0] = THRESH_ABSOLUTE ;
     }
  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;

  prefs.satur_level = data.saturation ;
  FitsHeader head(filename, RO);
  // if you have a problem with next line, try to fix it in
  // the relevant virtual instrument rather than hardcoding a value
  prefs.pixel_scale = head.KeyVal("TOADPIXS");
  prefs.gain = head.KeyVal("TOADGAIN");

  cout << "saturation provided to sextractor : " << prefs.satur_level  << std::endl;
  cerr << "Pixel Scale : " << prefs.pixel_scale << endl ;
  cerr << "Gain : " << prefs.gain << endl ;
  cerr << prefs.param_name << endl ;



  makeit();

  Fond = thefield1.backmean;
  SigmaFond= thefield1.backsig ;
  
  // free the mallocs:
  for (unsigned k=0; k<toFree.size(); ++k) free(toFree[k]);

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
// this routine can only be called with a 
static int
_SEStarListMake(const AllForSExtractor & data,
	       SEStarList &List, double & Fond, 
	       double &SigmaFond,
	       Image * pmask_sat)
{
  Liste_de_SEStar= &List;
  pmask = pmask_sat ;
  if (file_exists(data.FitsFileName.c_str()) )
    {
     sex_proc(data, Get_SEStarList, Get_PixelSat,
	      Fond, SigmaFond);
     MakeNumbers(List);
     return 1;
    }
  return 0;
}

// encapsulage niveau 1



void ForSExtractor::DecompressIfNeeded()
{
  string* names[2] = {&FitsFileName, &FitsWeightName};
  for (unsigned k=0; k<sizeof(names)/sizeof(names[0]); ++k)
    {
      string &name = *(names[k]);
      if (name != "") 
	{
	  string outFileName = AddSlash(TempDir)+UniqueName+"."
	    +CutExtension(BaseName(name))+".fits";
	  name = DecompressImageIfNeeded(name, outFileName, ToRemove);
	}
    }
}


ForSExtractor::~ForSExtractor()
{
  if (ToRemove != "")
    cout << " removing " << ToRemove << endl;
  RemoveFiles(ToRemove);
}


void
AllForSExtractor::Print()
{
  ForSExtractor::Print();
}

void
ForSExtractor::Print()
{
  cout << "FitsFileName " << FitsFileName << endl ;
  cout << "FitsMaskName  " << FitsMaskName << endl ; 
  cout << "FitsBackName  " << FitsBackName << endl   ;
  cout << "FitsMiniBackName " << FitsMiniBackName << endl  ;
  cout << "FitsSegmentationName " << FitsSegmentationName << endl  ;
  cout << "saturation " << saturation << endl   ;
  cout << "back_type_manual : " ;
  if (back_type_manual == true)
    cout << "set " << endl << " a constant value = " << backmean 
	 << " is taken for the background "  << endl   ;
  else
    cout << "not set " << endl ;
  cout << "sigma_back : " ;
  if (sigma_back < 0 )
    cout << "not specified" << endl ;
  else
    cout << sigma_back  << endl  ;
}

bool
AllForSExtractor::FillFromEnvironnement()
{
  string sextractor_dir = ".";
  if (getenv("TOADSCARDS")) {
    sextractor_dir = getenv("TOADSCARDS");
    cout << " using SExtractor datacards from "  << sextractor_dir << endl;
  }
  else
    cout << " TOADSCARDS not defined, use local directory  " << endl;
  SexConfigFileName = sextractor_dir+"/default.sex";
  if (getenv("SECARD1"))
    {
      SexConfigFileName = sextractor_dir+"/default1.sex";
      cout << "Taking datacard: " << SexConfigFileName << endl ;
    }
  SexParamName = sextractor_dir+"/default.param" ;
  SexNNWName   = sextractor_dir+"/default.nnw";
  SexFilterName = sextractor_dir+"/default.conv" ;
  if (!file_exists(SexConfigFileName.c_str()) ||  
      !file_exists(SexParamName.c_str()) ||
      !file_exists(SexNNWName.c_str()) ||
      !file_exists(SexFilterName.c_str()) )
    return(false);
  else
    return(true);
}



int
SEStarListMake(const ForSExtractor & shortdata,
	       SEStarList &List, double & Fond, 
	       double &SigmaFond,
	       Image * pmask_sat)
{
  AllForSExtractor data(shortdata);
  if (!data.FillFromEnvironnement())
    return 0 ;
  data.DecompressIfNeeded();
  int status = _SEStarListMake(data, List, Fond, SigmaFond,
			      pmask_sat);
  return(status);
}

/*************** SExtractor with 2 images ******************/
static
std::string DecompressIfNeeded(const std::string &InFileName,
					      const  std::string &tmpdir,
					      const  std::string & unique,
					      std::string &ToRemove) 
{
  std::string outFileName =  AddSlash(tmpdir)+unique+"."
    +CutExtension(BaseName(InFileName))+".fits";
  return DecompressImageIfNeeded(InFileName, outFileName, ToRemove);
}


static 
int 
sex_proc_2(const AllForSExtractor & data,
	 _SexStarFill StarFill, 
	   double & Fond_0, double &SigmaFond_0, 
	   double & Fond_1, double &SigmaFond_1, bool weight_from_measurement_image)
{
  std::vector<void *> toFree;

  // nombre de check image
  int nchecks=0;

  std::string toRemove;
  memset(&prefs,0,sizeof(prefs));
  prefs.pipe_flag = 0;
  prefs.nimage_name = 2; // 2 images pour detection & photometrie


  strcpy(prefs.prefs_name, data.SexConfigFileName.c_str());

  
  char filename_0[512];
  cerr << "# " << data.FitsFileName_0 << endl ;
  std::string imageName_0 = DecompressIfNeeded(data.FitsFileName_0, data.TempDir_0, data.UniqueName_0, toRemove);
  strcpy(filename_0,imageName_0.c_str());
  prefs.image_name[0] = filename_0;
  cerr << prefs.image_name[0]  << endl ;

  char filename_1[512];
  cerr << "# " << data.FitsFileName_1 << endl ;
  std::string imageName_1 = DecompressIfNeeded(data.FitsFileName_1,  data.TempDir_1, data.UniqueName_1, toRemove); 
  strcpy(filename_1,imageName_1.c_str());
  prefs.image_name[1] = filename_1;
  cerr << prefs.image_name[1]  << endl ;

  cerr << "# Images en entree : " << prefs.image_name[0] << " " << prefs.image_name[1] << endl ;

  readprefs(prefs.prefs_name, NULL, NULL, 0);
  
  strcpy(prefs.param_name, data.SexParamName.c_str());
  strcpy(prefs.nnw_name,data.SexNNWName.c_str() );
  // le filtre est lu dans la datacard (contrairement au traitement dans le cas 1 image, cf sex_proc.
  //strcpy(prefs.filter_name,data.SexFilterName.c_str());

  size_t namesize;
  if (data.back_type_manual)
    {
      prefs.back_type[0]=BACK_ABSOLUTE;
      prefs.back_type[1]=BACK_ABSOLUTE;
      prefs.nback_type=2;
      prefs.back_val[0] = data.backmean;
      prefs.back_val[1] = data.backmean;
      prefs.nback_val=2;
    }
  else
    {
      cerr << "Calcul et soustraction du fond par SE" << endl ;
      prefs.back_type[0]=BACK_RELATIVE;
      prefs.back_type[1]=BACK_RELATIVE;
      prefs.nback_type=2;
      prefs.back_val[0] = data.backmean;
      prefs.back_val[1] = data.backmean;
      prefs.nback_val=2;
      // l acarte de fond sera celle de l'image de mesure
      //(les check image ne concernent que l'image de mesure)
      if ((namesize=strlen(data.FitsMiniBackName.c_str())) != 0 )
	{  
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks], data.FitsMiniBackName.c_str() );
	  prefs.check_type[nchecks] = CHECK_MINIBACKGROUND;
	  nchecks++;
	}
      if ((namesize=strlen(data.FitsBackName.c_str())) != 0 )
	{  
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks], data.FitsBackName.c_str() );
	  prefs.check_type[nchecks] = CHECK_BACKGROUND;
	  nchecks++;
	}
    }
     if ((namesize = strlen(data.FitsSegmentationName.c_str())) != 0 )
	{  
	  prefs.check_name[nchecks] = (char *)calloc(namesize+2,1);
	  strcpy(prefs.check_name[nchecks], data.FitsSegmentationName.c_str() );
	  cerr << "## Segmentation image requested " << data.FitsSegmentationName << endl;
	  prefs.check_type[nchecks] = CHECK_SEGMENTATION;
	  nchecks++;
	}
     
  /*Prefs en sortie quand une seule weight image donnee en dble mode:
    prefs.nwimage_name : 1
    prefs.nweight_type : 1
    prefs.nweight_thresh : 0
    prefs.wimage_name[0] : weight.fits
    prefs.wimage_type[0] : 4
    weight_flag  : 1
    dweight_flag  : 1
    weightgain_flag  : 1
prefs.nweight_thresh  : 0
prefs.weight_thresh[0]: 0
prefs.weight_thresh[1]: 0
*/

  /*Prefs en sortie quand 2 weight maps donnees:
prefs.nwimage_name : 2
prefs.nweight_type : 2
prefs.nweight_thresh : 0
prefs.wimage_name[0] : weight.fits
prefs.wimage_type[0] : 4
prefs.wimage_name[1] : weight2.fits
prefs.wimage_type[1] : 4
weight_flag  : 1
dweight_flag  : 1
weightgain_flag  : 1
prefs.nweight_thresh  : 0
prefs.weight_thresh[0]: 0
prefs.weight_thresh[1]: 0
*/


  prefs.ncheck_name = nchecks ;
  prefs.ncheck_type = nchecks ;

  // ie weight from measurement image used  for measurement
  // and detection
  if ( weight_from_measurement_image   && (strlen(data.FitsWeightName_1.c_str())))
    {
      cerr << "Only Weight Map from measurement image : " << prefs.image_name[1]
	   << " will be used for both detection and  measurement" << endl ;
      prefs.wimage_name[0] = (char *)calloc(MAXCHAR,1);      
      toFree.push_back(prefs.wimage_name[0]);
      std::string weightName_1 = 
	DecompressIfNeeded(data.FitsWeightName_1, data.TempDir_1, data.UniqueName_1,toRemove);
      strcpy(prefs.wimage_name[0],weightName_1.c_str());
      prefs.nwimage_name = 1 ; 
      prefs.weight_type[0]=WEIGHT_FROMWEIGHTMAP;
      prefs.nweight_type=1;
      prefs.weight_flag = 1;
      prefs.dweight_flag = 1;// pas nec, subordonne a weight_type[0]!=NONE
      prefs.weightgain_flag=0;  
    }

  // on utilise des weight pour la detetction et la mesure
   if ( !weight_from_measurement_image  && (strlen(data.FitsWeightName_0.c_str())) && (strlen(data.FitsWeightName_1.c_str())))
    {
      prefs.wimage_name[0] = (char *)calloc(MAXCHAR,1);      
      toFree.push_back(prefs.wimage_name[0]);    
      std::string weightName_0 = 
	DecompressIfNeeded(data.FitsWeightName_0, data.TempDir_0, data.UniqueName_0,toRemove);
      strcpy(prefs.wimage_name[0],weightName_0.c_str());

      prefs.wimage_name[1] = (char *)calloc(MAXCHAR,1);         
      toFree.push_back(prefs.wimage_name[1]);      
      std::string weightName_1 = 
	DecompressIfNeeded(data.FitsWeightName_1, data.TempDir_1, data.UniqueName_1,toRemove);
      strcpy(prefs.wimage_name[1],weightName_1.c_str());
      prefs.nwimage_name = 2 ; 
      prefs.weight_type[0]=WEIGHT_FROMWEIGHTMAP;
      prefs.weight_type[1]=WEIGHT_FROMWEIGHTMAP;
      prefs.nweight_type=2;
      // mis dans prefs.c ?
      prefs.weight_flag = 1;
      prefs.dweight_flag = 1;
      prefs.weightgain_flag=0;
      // normalement la valeur mise dans prefs.c est la bonne pour le thresh,
      // ie : prefs.nweight_thresh = 2; prefs.weight_thresh[1]=prefs.weight_thresh[0]=0. ;      
    }
  
  cerr << "prefs.nwimage_name : " << prefs.nwimage_name << endl  ; 
   cerr << "prefs.nweight_type : " <<  prefs.nweight_type << endl  ; 
   cerr << "prefs.nweight_thresh : " <<  prefs.nweight_thresh << endl  ; 
  cerr << "prefs.wimage_name[0] : " <<  prefs.wimage_name[0] << endl ;
  cerr << "prefs.wimage_type[0] : " <<  prefs.weight_type[0] << endl ; 
  if (prefs.nweight_type >1)
    {
      cerr << "prefs.wimage_name[1] : " <<  prefs.wimage_name[1] << endl ;
      cerr << "prefs.wimage_type[1] : " <<  prefs.weight_type[1] << endl ;
    }  
  cerr << "weight_flag  : " <<  prefs.weight_flag<< endl ; /* do we weight ? */
  cerr << "dweight_flag  : " <<  prefs.dweight_flag	<< endl ; /* detection weight? */
  cerr << "weightgain_flag  : " <<  prefs.weightgain_flag<< endl ; /* weight gain? */
  cerr << "prefs.filter_name : " << prefs.filter_name << endl ;


  // Pas de flag image

  prefs.nimaflag= 0 ; prefs.nfimage_name = 0;
  
  
  prefs.user_ana = StarFill; 

  // no saturation map writen


  // detection level and photometric level in number of sigmas
  // they wiil be converted by SExtractor in photons,
  // using the sigma it will compute.
  // attention de n'en mettre qu'un sinon il le convertit en mag !

  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;
  // to set these levels directly in photons (or adus),
  // using a sigma specified by the user.
   if ( data.sigma_back > 1.e-10 )
     {
       prefs.dthresh[0] = prefs.dthresh[0] * data.sigma_back ;
       prefs.thresh[0] = prefs.thresh[0] * data.sigma_back ;
       prefs.thresh_type[0] = THRESH_ABSOLUTE;
     }
  cout << "DETECTION   THRESH " << prefs.dthresh[0] << endl ;
  cout << "PHOTOMETRIC THRESH " << prefs.thresh[0] << endl ;

  FitsHeader head(filename_1, RO);
  // le pixelscale n'est pas forcement la meme mais on s'en fiche
  // par contre il faut le gain de l'image de photometrie

  // if you have a problem with next line, try to fix it in
  // the relevant virtual instrument rather than hardcoding a value
  prefs.pixel_scale = head.KeyVal("TOADPIXS");
  prefs.gain = head.KeyVal("TOADGAIN");  
  prefs.satur_level = data.saturation ;
  
 

  //if (getenv("SEXTRACTOR_VERBOSE"))
  //prefs.verbose_type = WARN ;

  cout << "saturation provided to sextractor : " << prefs.satur_level  << std::endl;
  cerr << "Pixel Scale : " << prefs.pixel_scale << endl ;
  cerr << "Gain : " << prefs.gain << endl ;
  makeit();

  cerr << "prefs.nwimage_name : " << prefs.nwimage_name << endl  ; 
   cerr << "prefs.nweight_type : " <<  prefs.nweight_type << endl  ; 
   cerr << "prefs.nweight_thresh : " <<  prefs.nweight_thresh << endl  ; 
  cerr << "prefs.wimage_name[0] : " <<  prefs.wimage_name[0] << endl ;
  cerr << "prefs.wimage_type[0] : " <<  prefs.weight_type[0] << endl ; 
  if (prefs.nweight_type >1)
    {
      cerr << "prefs.wimage_name[1] : " <<  prefs.wimage_name[1] << endl ;
      cerr << "prefs.wimage_type[1] : " <<  prefs.weight_type[1] << endl ;
    } 
  cerr << "weight_flag  : " <<  prefs.weight_flag<< endl ; /* do we weight ? */
  cerr << "dweight_flag  : " <<  prefs.dweight_flag	<< endl ; /* detection weight? */
  cerr << "weightgain_flag  : " <<  prefs.weightgain_flag<< endl ; /* weight gain? */
 
  Fond_0 = thefield1.backmean;
  SigmaFond_0= thefield1.backsig ;
  Fond_1 = thefield2.backmean;
  SigmaFond_1= thefield2.backsig ;
  
  // free the mallocs:
  for (unsigned k=0; k<toFree.size(); ++k) free(toFree[k]);

  // cleanup temp files
  RemoveFiles(toRemove);

  return 1;
} 



//makes a SEStarList on 2 images
static int
_SEStarListMake_2(const AllForSExtractor & data,
	       SEStarList &List, double & Fond_0, 
	       double &SigmaFond_0, double & Fond_1, 
	       double &SigmaFond_1, bool weight_from_measurement_image)
{
  Liste_de_SEStar= &List;
  if (file_exists(data.FitsFileName_0.c_str()) && file_exists(data.FitsFileName_1.c_str()) )
    {
     sex_proc_2(data, Get_SEStarList, 
	      Fond_0, SigmaFond_0,Fond_1, SigmaFond_1, weight_from_measurement_image);
     MakeNumbers(List);
     return 1;
    }
  return 0;
}
int
SEStarListMake_2(const ForSExtractor & shortdata,
		 SEStarList &List, double & Fond_0, 
		 double &SigmaFond_0,double & Fond_1, 
		 double &SigmaFond_1,bool use_weight_from_image2)
{
  AllForSExtractor data(shortdata);
  if (!data.FillFromEnvironnement())
    return 0;
  int status = _SEStarListMake_2(data, List, Fond_0, SigmaFond_0, 
				 Fond_1, SigmaFond_1,use_weight_from_image2);
  return(status);
}

/********************************************/
// Pour pouvoir recuperer le background a partir de la mini carte


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

#ifdef OFF_T_CFITSIO
#undef OFF_T
#define OFF_T OFF_T_CFITSIO
#undef OFF_T_CFITSIO
#endif
