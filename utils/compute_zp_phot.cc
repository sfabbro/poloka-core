#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <string>
#include <getopt.h>
#include <string.h>

#include <stdio.h>
#include <math.h>
using namespace std;

#define SATUR_FLAG 4

#include "fitsimage.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include "dicstar.h"
#include "vutils.h"
#include "frame.h"

#include "imageutils.h"
#include "fileutils.h"


using namespace std;

// mag values criteria for star calib selection according to Regnault et al. 2009
int IsOK_MARS09(string band, double mag) ;

int main(int argc, char **argv)
{
 char c;
 string nom, nom_cat,  cat_calib;
  while ((c = getopt(argc, argv, "hi:C:c:b")) != -1) 
    {
      switch (c)
	{
	case 'h' :
	  //usage();
	  break;
	case 'i' :
	  nom = optarg ; 
	  break;
	case 'C' :
	  cat_calib = optarg ; 
	  break;
	case 'c' :
	  nom_cat = optarg ; 
	  break;
	default:
	  //usage();
	  cerr << "bad option " << endl ;
	}
    }
  double zporig = 0. ;

  FitsHeader head(nom,RW);
  string filter = StringToLower(string(head.KeyVal("TOADBAND"))) ;
  string bandname = "m"+filter ;
  cout << "band : " << bandname << endl ;
  DicStarList lcal(cat_calib);    
  DicStarList lcal_cut;    
  DicStarList stl(nom_cat);
  Gtransfo *Pix2RaDec;
  WCSFromHeader(head, Pix2RaDec); 
  cout << *Pix2RaDec << endl ;
  Frame raDecFrame = ApplyTransfo(Frame(head),*Pix2RaDec);
  Frame bigraDecFrame = raDecFrame.Rescale(1.1);
  lcal.ExtractInFrame(lcal_cut,bigraDecFrame );
  // 0.01 is the precision ot invertion (in pixels);
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01,Frame(head) );

  int n = 0 ;
  double *ZPtab_c = new double[lcal_cut.size()] ;
  zporig = head.KeyVal("ZP");

 
  for(DicStarIterator it = lcal_cut.begin(); it!= lcal_cut.end(); ++it)
    {
      DicStar *starc = *it ;

      // position de l'etoile de calib sur l'image (ra et dec dans x et y)
      double xx, yy;
      RaDec2Pix->apply(starc->x, starc->y, xx, yy);	 
      DicStar *mstar  = stl.FindClosest(xx,yy);

      if (mstar != NULL)
	{ 
	  double d = sqrt((mstar->x-xx)*(mstar->x-xx)+(mstar->y-yy)*(mstar->y-yy)) ;
	  double fap6 = mstar->getval("apfl6") ;		       	     
	  double efap6 = mstar->getval("eapfl6");
	  double mag = starc->getval(bandname) ;

	  
	  // test saturation
	  int flag = mstar->getval("flag");
	  if (flag & SATUR_FLAG ) continue ;


	  if (fap6 < 1) continue ;
	  if (d > 2 ) continue ;
	  if (mag < 1) continue ;

	  if ( IsOK_MARS09(filter, mag) < 0 )
	    continue ;

	  double emag = starc->getval("e"+bandname) ;
	  double zpc = mag+ 2.5*log10(fap6) ;
	  double ezpc = (efap6/fap6)*(efap6/fap6)+ emag*emag;

	  if (sqrt(ezpc)>1.)
	    continue ;

	  ZPtab_c[n] =  zpc;
	  n++;
	    
	}
    }
  int n_c=n  ;
  double sigzp_c= 0 ;
  double k = 3 ;
  int niter = 3 ;
  //!!!! ZPtab_c order will be preturbated.
  double zp_c  = -1 ;
  if (n_c>0) zp_c = clipmean(ZPtab_c, n_c , sigzp_c,k, niter);
  cout << "OLD ZP : " << zporig << endl ;
  cout << "ZPphot " << setprecision (6) << zp_c << endl ;
  cout << "sig_ZPphot " << setprecision (6) << sigzp_c << endl ;
  cout << "Nstars_ZPphot  " << n_c << endl ;
  if (n_c<=0) zp_c = -1 ;
  string comment = " w.r.t to phot. calibration catalog " ;
  head.AddOrModKey("ZP_PHOT",zp_c, comment.c_str());
  delete[] ZPtab_c ;
    
  delete Pix2RaDec ;
  delete RaDec2Pix ;

}
int IsOK_MARS09(string band, double mag)
{
  int ok = 1 ;
  if (band == "u" ) 
    {
      if (( mag < 17.5 ) || (mag > 20.5))
	ok = -1 ;
    }
  if (band == "g" ) 
    {
      if (( mag < 17.5 ) || (mag > 22))
	ok = -1 ;
    }
  if (band == "r" ) 
    {
      if (( mag < 17.5 ) || (mag > 21.5))
	ok = -1 ;
    }
     
  if (band == "i" ) 
    {
      if (( mag < 17.5 ) || (mag > 20.8))
	ok = -1 ;
    } 
  if (band == "z" ) 
    {
      if (( mag < 16 ) || (mag > 19))
	ok = -1 ;
    } 
  return(ok);
}
