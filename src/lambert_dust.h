#ifndef LAMBERT_DUST__H
#define LAMBERT_DUST__H


#include <fstream>
#include <iostream>
#include "lambert_tag.h"
#define WCS_GALACTIC    3 



#include "fitsimage.h"
#ifndef ZEA
#define ZEA 1
#define LAMB 0
#endif

class Lambert_transfo {

public:

Lambert_transfo(FitsImage* image);
virtual ~Lambert_transfo(){};
void LongLat2Pix(double Lon, double Lat, double &x, double &y);
void LongLat2Pix(double Lon, double Lat, int &x, int &y);
void RaDec2Pix(double Lon, double Lat, double &x, double &y);
void RaDec2Pix(double Lon, double Lat, int &x, int &y);

private:

int transmode;
double    crval1;
double    crval2;
double    crpix1;
double    crpix2;
int 	  naxis1;
int       naxis2;

int      nsgp;
double    scale;

double    cd1_1;
double    cd1_2;
double    cd2_1;
double    cd2_2;
double    lonpole;



};

enum Mem_Status{Auto =0, Min, Full};

class Lambert_Dust{

public:

Lambert_Dust();
virtual ~Lambert_Dust();
void LoadNorth();
void LoadSouth();
void ClearNorth(){if(NorthMap != NULL) {delete NorthMap; NorthMap=NULL;}};
void ClearSouth(){if(SouthMap != NULL) {delete SouthMap; SouthMap=NULL;}};
void ClearAll(){ ClearSouth(); ClearNorth(); return;};
void Set_Memory_Status(Mem_Status status);

double Value_From_RaDec(double Ra,double Dec);
double Value_From_Gal(double Lon, double Lat);


private:
string	NorthMapName,
	SouthMapName;
Mem_Status Memory_Status;
FitsImage* NorthMap;
FitsImage* SouthMap;
Lambert_transfo* NorthTransfo;
Lambert_transfo* SouthTransfo;
bool transfo_north_ok;
bool transfo_south_ok;
bool Interp;
int Interp_Level;

};


#endif
