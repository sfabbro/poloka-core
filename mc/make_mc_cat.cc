#include "reducedimage.h"
#include "sestar.h"
#include "dictfile.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include "fitsimage.h"
#include "host.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
//#define DEBUG 
#define MAX_DIST 2.0
#define MAX_DELTA_MAG 0.3

static void usage(const char *pgname)
{
  cerr << pgname << " <poloka_cat> <gwyn_cat>  " << endl ;
  exit(1);
}
 
 

 
int main(int argc, char **argv)
{
 if (argc != 4) { usage(argv[0]) ; exit(1);}

string gwyn_cat_name = argv[2];
string ImageName = argv[1];
string out_name= argv[3];
  
ReducedImage image(ImageName);
cout << "loading reducedimage : " << ImageName << endl;
        string band;
        band=image.Band();
        if (band.size()>=1)
        {
                if(band[0]> 64 && band[0]< 91) band[0]=band[0]+32;
        }
        else {cerr << "Band error" << endl; exit(1);}
if(band!="i") { cerr << "wrong image : i filter needed" << endl; exit(1);}
string mag_tag = band+"Mega";
//if (!image.HasCatalog()) image.MakeCatalog();
 
  DictFile gwyn_cat(gwyn_cat_name);
  FitsImage image_fits(image.FitsName());


Frame current_Frame(image_fits, WholeSizeFrame);
Gtransfo *Pix2RaDec;
if (!WCSFromHeader(/*dynamic_cast <FitsHeader> */image_fits, Pix2RaDec))
        {
                cout << "cannot handle "<<image.FitsName() <<" without a WCS " << endl;
		exit(1);
        }
Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, current_Frame);
 



SEStarList catalog(image.CatalogName());

		#ifdef DEBUG
		cout <<catalog.size() <<" elements "<< endl;
		cout << "loading standard list... : "<< Get_Ref_Cat_Name(field) <<endl;
		#endif 
 int in_frame =0;
 int  matched =0;
 int  matched2 =0; 
 int no_matched=0;
 int too_faint=0;
 double mag_lim=25.5;
 double total_delta_mag =0.0;
 int too_faint_no_matched=0;
 int matched_array[5][7]; 
 for ( int i=0 ; i<5;i++) for ( int j = 0 ;j<7 ;j++)matched_array[i][j]=0; 

HostList res_list;
for(DictFileCIterator entry=gwyn_cat.begin();entry!=gwyn_cat.end();++entry) 
{

   
    double Ra =entry->Value("Ra"); // ra (deg)
    double Dec=entry->Value("Dec"); // dec (deg)
    double x,y;
    RaDec2Pix->apply(Ra,Dec,x,y);
    if (!current_Frame.InFrame(x,y)) continue;
    double z = entry->Value("z_phot");
    if(   z>=1.1 || z<=0.3) continue ;    
    
    in_frame++;
    SEStar *closest = catalog.FindClosest(x,y);
    double dist = (closest->X()-x)*(closest->X()-x)+(closest->Y()-y)*(closest->Y()-y);
    dist = sqrt(dist);
    double se_mag = 30.15 - 2.5*log10(closest->flux);
    double gwyn_mag = entry->Value(mag_tag);
    double delta_mag = abs(se_mag-gwyn_mag);
    
   if(   gwyn_mag > 50.0) continue ;
    

    if(gwyn_mag >mag_lim && dist>=MAX_DIST){too_faint_no_matched++;}    
    if (dist>=MAX_DIST) {no_matched++;}
    if(gwyn_mag >mag_lim){too_faint++;}
    
    
    if (dist<5.0)
    {   

    	matched_array[(int)dist][0]++;
       	if (delta_mag < 0.5) matched_array[(int)dist][(int)(10.0*delta_mag)+1]++;
    	else matched_array[(int)dist][6]++;
    }
    
    
    
    if ( dist< MAX_DIST && gwyn_mag < 24.0) { 
    total_delta_mag+=se_mag-gwyn_mag;
    matched2++;
    }
    
    
    if (/*delta_mag < MAX_DELTA_MAG &&*/ dist < MAX_DIST ) 
    { 
    total_delta_mag+=se_mag-gwyn_mag;
    matched++;
       	Host *p = new Host(Ra,Dec);
	p->Z()= z;
	p->Err_Z()=0.01;
	p->Gal_Type()= entry->Value("type");
	p->A()=closest->A();
	p->B()=closest->B();
	p->Theta()=closest->Gyr_Angle();
	p->I_Mag()=entry->Value(mag_tag);
	res_list.push_back(p);
	}
    if (dist > MAX_DIST && gwyn_mag > 25.5) 
    { 
       	Host *p = new Host(Ra,Dec);
	p->Z()= z;
	p->Err_Z()=0.01;
	p->Gal_Type()= entry->Value("type");
	p->A()=1.0;
	p->B()=1.0;
	p->Theta()=0.0;
	p->I_Mag()=entry->Value(mag_tag);
	res_list.push_back(p);
	}





}
res_list.write(out_name);  
cout <<"catalog : " <<  catalog.size() <<" elements "<< endl; 
cout <<"gwyn's cat : " << gwyn_cat.size() <<" elements "<< endl;
cout <<"in frame : " << in_frame << " elements "<< endl;
cout << "no det within "<<MAX_DIST <<" pix : " << no_matched << " elements "<< endl;
cout << "with mag > "<<mag_lim <<" : " << too_faint << " elements "<< endl;
cout << "with mag > "<<mag_lim <<" and no det within "<<MAX_DIST <<" pix : " << too_faint_no_matched << " elements "<< endl;
cout << "delta mag moyen : " << (total_delta_mag/(double)matched2) << endl;

for(int i = 0 ; i<7 ;i++)
cout << matched_array[0][i] << " | "<< matched_array[1][i] << " | "<< matched_array[2][i] << " | "
     << matched_array[3][i] << " | "<< matched_array[4][i] << endl;
cout << matched << " galaxies ok with maxdist = " << MAX_DIST << endl;    
 
 
} 
  
  
  
  
  
  
  
  
  
  
  
  
  
