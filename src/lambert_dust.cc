#include "lambert_dust.h"
#include "datacards.h"
#include "wcscon.h"
#include "fileutils.h"
#include "imageinterpolation.h"
#include "lambert_tag.h"
//#define DEBUG
#ifndef ZEA
#define ZEA 1
#define LAMB 0
#endif
/***************************************** Lambert_transfo ***********************************/





static void RaDec2LonLat(double Ra, double Dec, double &Lon, double &Lat)

{
Lon=Ra;
Lat=Dec;
wcscon (WCS_J2000, WCS_GALACTIC, 0.0, 0.0, &Lon, &Lat, 0.0);
#ifdef DEBUG
cout << " Ra " << Ra << " Dec " << Dec << " -----> Long " << Lon  << " Lat " << Lat << endl; 
#endif

}


static void lambert_lb2xy
  (double    gall,   /* Galactic longitude */
   double    galb,   /* Galactic latitude */
   int       nsgp,   /* +1 for NGP projection, -1 for SGP */
   double    scale,  /* Radius of b=0 to b=90 degrees in pixels */
   double&  pX,     /* X position in pixels from the center */
   double&  pY)     /* Y position in pixels from the center */
{
   double   rho;
   const double dradeg = 180 / 3.1415926534;

   rho = sqrt(1. - nsgp * sin(galb/dradeg));
/* The following two lines were modified by Hans Schwengeler (17-Mar-1999)
   to get this to work on a Tur64 Unix 4.0E (DEC Alpha).  It appears that
   float and double on not the same on this machine.
   *pX = rho * cos(gall/dradeg) * scale;
   *pY = -nsgp * rho * sin(gall/dradeg) * scale;
*/
   pX = rho * cos(((double)gall/dradeg)) * scale;
   pY = -nsgp * rho * sin(((double)gall/dradeg)) * scale;
}





Lambert_transfo::Lambert_transfo(FitsImage* image)
{

  string  ctype1=image->KeyVal(label_ctype1);
  string  ctype2=image->KeyVal(label_ctype2);
  
  
   crval1=image->KeyVal(label_crval1);
   crval2=image->KeyVal(label_crval2);
   crpix1=image->KeyVal(label_crpix1);
   crpix2=image->KeyVal(label_crpix2);
   
   naxis1=image->KeyVal(label_naxis1);
   naxis2=image->KeyVal(label_naxis2); 
   
   
//mode LAMBERT
   if (strcmp(ctype1.c_str(), "LAMBERT--X")  == 0 &&
       strcmp(ctype2.c_str(), "LAMBERT--Y")  == 0)
	{
   		transmode = LAMB;
      		nsgp=image->KeyVal(label_lam_nsgp);
      		scale=image->KeyVal(label_lam_scal);
      
      }
      
//mode ZEA
    else if (strcmp(ctype1.c_str(), "GLON-ZEA")  == 0 && 
             strcmp(ctype2.c_str(), "GLAT-ZEA") == 0)
	{
	
		transmode = ZEA;     
      		if (image->HasActualKey(label_cdelt1) &&
            	image->HasActualKey(label_cdelt2))	   
	   
      		{
          	cd1_1 = image->KeyVal(label_cdelt1 );
          	cd1_2 = 0.0;
          	cd2_1 = 0.0;
          	cd2_2 = image->KeyVal(label_cdelt2 );
      		} 
      		else 
      		{
         	cd1_1=image->KeyVal(label_cd1_1);
         	cd1_2=image->KeyVal(label_cd1_2);
         	cd2_1=image->KeyVal(label_cd2_1);
         	cd2_2=image->KeyVal(label_cd2_2);
      		}
     
     
      		if (image->HasActualKey(label_lonpole))lonpole = image->KeyVal(label_lonpole);
      		else  lonpole = 180.0;
		}
	else { cerr << "unknown projection: ABORT!!" << endl; exit(1);}

}


void Lambert_transfo::LongLat2Pix
  (double    gall,   /* Galactic longitude */
   double    galb,   /* Galactic latitude */
   double&  pX,     /* X position in pixels from the center */
   double&  pY)     /* Y position in pixels from the center */
{
  #ifdef DEBUG
cerr << "LongLat2Pix double"<< endl;
#endif
   double    xr,yr; 
   const double dradeg = 180.0 / 3.1415926534;
switch (transmode)
{
 
case LAMB : 
#ifdef DEBUG
cerr << "LAMB"<< endl;
#endif
      lambert_lb2xy(gall, galb, nsgp, scale, xr, yr);
      pX = xr + crpix1 - crval1 - 1.0;
      pY = yr + crpix2 - crval2 - 1.0;
      break;
      
      
case ZEA : 
#ifdef DEBUG
cerr << "ZEA"<< endl;
#endif
      double theta,phi, Rtheta,denom;
      /* ROTATION */
      /* Equn (4) - degenerate case */
      if (crval2 > 89.9999) {
         theta = galb;
         phi = gall + 180.0 + lonpole - crval1;
      } else if (crval2 < -89.9999) {
         theta = -galb;
         phi = lonpole + crval1 - gall;
      } else {
         cout <<"ERROR: Unsupported projection!!!   Assume it's an NGP projection ..."<< endl;
         theta = galb;
         phi = gall + 180.0 + lonpole - crval1;
      }

      /* Put phi in the range [0,360) degrees */
      phi = phi - 360.0 * floor(phi/360.0);

      /* FORWARD MAP PROJECTION */
      /* Equn (26) */
      Rtheta = 2.0 * dradeg * sin((0.5 / dradeg) * (90.0 - theta));

      /* Equns (10), (11) */
      xr = Rtheta * sin(phi / dradeg);
      yr = - Rtheta * cos(phi / dradeg);

      /* SCALE FROM PHYSICAL UNITS */
      /* Equn (3) after inverting the matrix */
      denom = cd1_1 * cd2_2 - cd1_2 * cd2_1;
      pX = (cd2_2 * xr - cd1_2 * yr) / denom + (crpix1 - 1.0);
      pY = (cd1_1 * yr - cd2_1 * xr) / denom + (crpix2 - 1.0);
	break;
	
default:

      pX = -99.0;
      pY = -99.0;

}
  #ifdef DEBUG
cerr << "LongLat2Pix double"<< endl;
#endif
}




void Lambert_transfo::LongLat2Pix
  (double    gall,   /* Galactic longitude */
   double    galb,   /* Galactic latitude */
   int&  pIX,     /* X position in pixels from the center */
   int&  pIY)     /* Y position in pixels from the center */
{



double    xr;
double    yr;

LongLat2Pix(gall,galb,xr,yr);

pIX = (int)(xr + 0.5);
pIY = (int)(yr + 0.5);

   /* Force bounds to be valid at edge, for ex at l=0,b=0 */

if (pIX >= naxis1) pIX = naxis1 - 1;
if (pIY >= naxis2) pIY = naxis2 - 1;
}




void Lambert_transfo::RaDec2Pix
  (double    Ra,   /* Galactic longitude */
   double    Dec,   /* Galactic latitude */
   double&  pX,     /* X position in pixels from the center */
   double&  pY)     /* Y position in pixels from the center */
{

double a=Ra;
double b=Dec;

wcscon (WCS_J2000, WCS_GALACTIC, 0.0, 0.0, &a, &b, 0.0);
LongLat2Pix(a, b, pX, pY);

cout << " Ra " << Ra << " Dec " << Dec << " -----> Long " << a  << " Lat " << b 
                                       << " -----> x "    << pX << " y "   << pY << endl; 

}



void Lambert_transfo::RaDec2Pix
  (double    Ra,   /* Galactic longitude */
   double    Dec,   /* Galactic latitude */
   int&  pIX,     /* X position in pixels from the center */
   int&  pIY)     /* Y position in pixels from the center */
{
double a=Ra;
double b=Dec;

wcscon (WCS_J2000, WCS_GALACTIC, 0.0, 0.0, &a, &b, 0.0);
LongLat2Pix(a, b, pIX, pIY);

cout << " Ra " << Ra << " Dec " << Dec << " -----> Long " << a  << " Lat " << b 
                                       << " -----> x "    << pIX << " y "   << pIY << endl; 

}









/*********************************************** Lambert_Dust ************************************/

void Lambert_Dust::LoadNorth()
{
#ifdef DEBUG
cerr << "LoadNorth"<< endl;
#endif
if(NorthMap==NULL) NorthMap = new FitsImage(NorthMapName);
if(!transfo_north_ok) {NorthTransfo = new Lambert_transfo(NorthMap);transfo_north_ok = true ;}
#ifdef DEBUG
cerr << "LoadNorth"<< endl;
#endif
}

void Lambert_Dust::LoadSouth()


{
#ifdef DEBUG
cerr << "LoadSouth"<< endl;
#endif
if(SouthMap==NULL) SouthMap = new FitsImage(SouthMapName);
if(!transfo_south_ok) {SouthTransfo = new Lambert_transfo(SouthMap);transfo_south_ok = true ;}
#ifdef DEBUG
cerr << "LoadSouth"<< endl;
#endif
}



Lambert_Dust::Lambert_Dust()
{

#ifdef DEBUG
cerr << "consttructor begin"<< endl;
#endif
 NorthTransfo = NULL;
 SouthTransfo = NULL;
transfo_south_ok =false;
transfo_north_ok =false;
SouthMap = NULL;
NorthMap = NULL;
 
   
   
	string cardname = getenv("DUST_CARD");
	
#ifdef DEBUG	
	cout << "Reading card : " << cardname << endl;
#endif	
	
	if (cardname == "") { cerr << " Need environement variable 'DUST_CARD'  abort!" << endl; exit(1);}
   
  if( FileExists(cardname))
  {
  DataCards card(cardname);
  if (!card.HasKey("NORTH_MAP")|| !card.HasKey("SOUTH_MAP")) { cerr << "missing key!    abort!"<< endl; exit(1);}
  unsigned int memode=card.IParam("MEMORY_MODE");
  Memory_Status = Auto ;
  if (memode == 2) Memory_Status = Full ;
  else if(memode == 1) Memory_Status = Min;

  NorthMapName = card.SParam("NORTH_MAP");
  SouthMapName = card.SParam("SOUTH_MAP");
  if(card.IParam("INTERP_MODE")==0) Interp = false ;
  else Interp = true ;
  Interp_Level = card.IParam("INTERP_LEV");
  
  
#ifdef DEBUG  
  cout << "NorthMapName : "<< NorthMapName <<endl;
  cout << "SouthMapName : "<< SouthMapName <<endl;
  if (Interp) cout << "Interpolation with level : "<< Interp_Level  <<endl;
  else cout << "No interpolation" << endl;  
#endif   


  if(!IsFits (NorthMapName)){ cerr << "error NorthMap"<< endl; exit(1);}
  if(!IsFits (SouthMapName)){ cerr << "error SouthMap"<< endl; exit(1);}
  if(Memory_Status==Full){LoadNorth(); LoadSouth();}
  
  
#ifdef DEBUG
cerr << "consttructor end"<< endl;
#endif
  
  return;



  }
cerr << "Lambert_Dust constructor : No datacard found! abort! " << endl;
exit(1);

}




Lambert_Dust::~Lambert_Dust()
{
   ClearAll();
   if (transfo_north_ok){ delete NorthTransfo; NorthTransfo = NULL;transfo_north_ok=false;}
   if (transfo_south_ok){ delete SouthTransfo; SouthTransfo = NULL;transfo_south_ok=false;}

}

void Lambert_Dust::Set_Memory_Status(Mem_Status status)
{
Memory_Status = status;
switch(status)
{
case Min: ClearAll(); break;
case Full: LoadNorth(); LoadSouth(); break;
case Auto: break;
}
}





double Lambert_Dust::Value_From_Gal(double Lon, double Lat)

{
   #ifdef DEBUG
cerr << "Value_From_Gal"<< endl;
#endif 
 double res;  
FitsImage* map;
Lambert_transfo* transfo;
if(Lat>0.0)  { LoadNorth(); map = NorthMap; transfo= NorthTransfo;}
else         { LoadSouth(); map = SouthMap; transfo= SouthTransfo;}


 

if (Interp == false) 
	       {  /* NEAREST PIXELS */

		  int x,y;		  
		  transfo->LongLat2Pix(Lon, Lat,x,y);
		  #ifdef DEBUG
		  cerr << "x : "<< x << "     y : "<< y   << endl;
		  #endif
		  res= (double)((*map)(x,y));

               } 
else 
	       {  /* INTERPOLATE */
		  double x,y;
		  transfo->LongLat2Pix(Lon, Lat,x,y);
		  #ifdef DEBUG
		  cerr << "x : "<< x << "     y : "<< y   << endl;
		  #endif
		  res=  (double)Interpolate(*map,x,y,Interp_Level);
               }  /* -- END NEAREST PIXEL OR INTERPOLATE -- */


if (Memory_Status==Min) ClearAll();
#ifdef DEBUG
cerr << "Value_From_Gal"<< endl;
#endif 
return res;

}




double Lambert_Dust::Value_From_RaDec(double Ra,double Dec)

{
  #ifdef DEBUG
cerr << "Value_From_RaDec"<< endl;
#endif
double Lon ,Lat;
RaDec2LonLat(Ra,Dec,Lon,Lat);
double res = Value_From_Gal( Lon, Lat);
#ifdef DEBUG
cerr << "Value_From_RaDec"<< endl;
#endif
return res;
}



