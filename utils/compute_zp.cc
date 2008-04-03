//#include "standardstar.h"
#include "fitsimage.h"
#include "sestar.h"
#include "reducedimage.h"
#include "dictfile.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
//#define DEBUG 
using namespace std;
 

// need REFCAT variable
static string Get_Ref_Cat_Name(const string field)
{
string name;
//if (!(name = getenv("REFCAT"))) {cerr << " var REFCAT not present" << endl; exit(1);}
name = getenv("REFCAT");
name +=  "/" ;
name += field;
name += "_3.list";


return (name);

 
}


static void usage(const char *pgname)
{
  cerr << pgname << " <image name> " << endl ;
  exit(1);
}
 
 

 
int main(int argc, char **argv)
{
 if (argc != 2) { usage(argv[0]) ; exit(1);}
string ImageName = argv[1];
#ifdef DEBUG
cout << "loading reducedimage : " << ImageName << endl;
#endif  
  
ReducedImage image(ImageName);
        string band;
        band=image.Band();
        if (band.size()>=1)
        {
                if(band[0]> 64 && band[0]< 91) band[0]=band[0]+32;
        }
        else {cerr << "Band error" << endl; exit(1);}
//cout << band << endl;
//if (!image.HasCatalog()) image.MakeCatalog();
FitsImage image_fits(image.FitsName());


Frame current_Frame(image_fits, WholeSizeFrame);
Gtransfo *Pix2RaDec;
if (!WCSFromHeader(/*dynamic_cast <FitsHeader> */image_fits, Pix2RaDec))
        {
                cout << "cannot handle "<<image.FitsName() <<" without a WCS " << endl;
		exit(1);
        }
Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, current_Frame);
 




string field;
 if (image_fits.HasActualKey("OBJECT")) field = string(image_fits.KeyVal("OBJECT"));
	else { cerr << "OBJECT key not present on header " << endl; exit(1);}
 
		#ifdef DEBUG
		cout << "loading catalog : " <<image.CatalogName() << "...  :" ;
		#endif 

SEStarList catalog(image.CatalogName());

		#ifdef DEBUG
		cout <<catalog.size() <<" elements "<< endl;
		cout << "loading standard list... : "<< Get_Ref_Cat_Name(field) <<endl;
		#endif 

  
DictFile ref_stars(Get_Ref_Cat_Name(field) );

int count_ok = 0;
double somme_zp[3]={0.0,0.0,0.0};
double somme_poid[3]={0.0,0.0,0.0};
double somme_err[3]={0.0,0.0,0.0};
double ZP[3][4]={{0.0,0.0,0.0,0.0},
		 {0.0,0.0,0.0,0.0},
		 {0.0,0.0,0.0,0.0}};// zp / err / chi2

for(DictFileCIterator entry=ref_stars.begin();entry!=ref_stars.end();++entry) {

   
    double Ra =entry->Value("x"); // ra (deg)
    double Dec=entry->Value("y"); // dec (deg)
    double x,y;
    RaDec2Pix->apply(Ra,Dec,x,y);
    if (!current_Frame.InFrame(x,y)) continue;
    SEStar *closest = catalog.FindClosest(x,y);
    double dist = (closest->X()-x)*(closest->X()-x)+(closest->Y()-y)*(closest->Y()-y);
    dist = sqrt(dist);
    if (dist>1.5) continue;
    count_ok++;
    string mag_key =  "m" + band;
    string emag_key =  "em" + band;
    double mag=entry->Value(mag_key);
    double emag=entry->Value(emag_key);
    double flux[3];
    double eflux[3];
     flux[0] = closest->flux;
     eflux[0] = closest->EFlux();
     flux[1] = closest->Flux_auto();
     eflux[1] = closest->Eflux_auto();
     flux[2] = closest->Flux_isocor();
     eflux[2] = closest->Eflux_isocor();
    
     for(int i = 0; i < 3;i++){   
    double zp = mag +2.5*log10(flux[i]);
    double zp_err= sqrt(emag*emag + 2.5/log(10.0)*eflux[i]/flux[i]*2.5/log(10.0)*eflux[i]/flux[i]);
    somme_zp[i] = somme_zp[i] + zp/zp_err*zp_err;
    somme_poid[i] = somme_poid[i] +1.0/zp_err*zp_err;
    somme_err[i] = somme_err[i] + zp_err*zp_err;
    }
    #ifdef DEBUG
    cout << zp << "   " << dist << endl;
    #endif

     

  }
      for(int i = 0; i < 3;i++)
      { 
      ZP[i][0]=  somme_zp[i]/somme_poid[i];
      ZP[i][1]= sqrt(somme_err[i]/(double)count_ok);

       } 
  
  
  
for(DictFileCIterator entry=ref_stars.begin();entry!=ref_stars.end();++entry) {

   
    double Ra =entry->Value("x"); // ra (deg)
    double Dec=entry->Value("y"); // dec (deg)
    double x,y;
    RaDec2Pix->apply(Ra,Dec,x,y);
    if (!current_Frame.InFrame(x,y)) continue;
    SEStar *closest = catalog.FindClosest(x,y);
    double dist = (closest->X()-x)*(closest->X()-x)+(closest->Y()-y)*(closest->Y()-y);
    dist = sqrt(dist);
    if (dist>1.5) continue;
    string mag_key =  "m" + band;
    string emag_key =  "em" + band;
    double mag=entry->Value(mag_key);
    double emag=entry->Value(emag_key);
    double flux[3];
    double eflux[3];
     flux[0] = closest->flux;
     eflux[0] = closest->EFlux();
     flux[1] = closest->Flux_auto();
     eflux[1] = closest->Eflux_auto();
     flux[2] = closest->Flux_isocor();
     eflux[2] = closest->Eflux_isocor();
    
     for(int i = 0; i < 3;i++)
     {   
    	double zp = mag +2.5*log10(flux[i]);
    	double zp_err= sqrt(emag*emag + 2.5/log(10.0)*eflux[i]/flux[i]*2.5/log(10.0)*eflux[i]/flux[i]);
	ZP[i][2]=ZP[i][2]+((zp-ZP[i][0])*(zp-ZP[i][0]))/(zp_err*zp_err);
        ZP[i][3]=ZP[i][3]+((zp-ZP[i][0])*(zp-ZP[i][0]));
    }


     

  }  
  
//cout << "no fit!" << endl;  
//cout << count_ok << " stars used"<< endl; 
cout  <<  ZP[0][0] << " +- " << ZP[0][1] << " \tchi2/N : " << ZP[0][2]/(double)count_ok << " variance : " << sqrt(ZP[0][3]/(double)count_ok)  << endl;
//cout << "with aper_flux : " <<  ZP[1][0] << " +- " << ZP[1][1] << " \tchi2/N : " << ZP[1][2]/(double)count_ok << " variance : " << sqrt(ZP[1][3]/(double)count_ok)  <<endl;
//cout << "with isoc_flux : " <<  ZP[2][0] << " +- " << ZP[2][1] << " \tchi2/N : " << ZP[2][2]/(double)count_ok << " variance : " << sqrt(ZP[2][3]/(double)count_ok)  <<endl;

delete RaDec2Pix;
delete Pix2RaDec;
return(0);




}
