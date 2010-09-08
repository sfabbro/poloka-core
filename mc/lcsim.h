#ifndef SIMLC__H
#define SIMLC__H

#include <fstream>
#include <iostream>
#include <vector>
#include "fakeobject.h"
#include "fakesnia.h"
#include "reducedimage.h"

#include "subimage.h"
//#include "sestar.h"
#include "datacards.h"
#include "simulation.h"
#include "imagematch.h"
#include "imagemanager.h"
#include "host.h"

class to_random
{
public :
int Mode;
double A;
double B;
to_random(){};
to_random(int mode, double a, double b)
{
Mode=mode; 
if(mode==1){A=min(a,b);B=max(a,b);}
else {A=a;B=b;}
};
};

typedef enum HostMode{From_SE =0, From_photoz};


class lcsim_card {
public:
string Card_Name;
MethodeSim Method;
HostMode Host_Mode;// **************
int NumberOfFakes;


double imag_to_red_poly[5];
double k_min;
double k_max;
to_random H_ext;
to_random H_Rv;
to_random Stretch;
to_random Redshift;
to_random Dispersion;
to_random Color;
to_random Alpha;
to_random Beta;
to_random MW_Rv;

string instrument;
double omegam ;
double omegax ;
double w1 ;
double w0 ;
double H0 ;


int delta_x ;
int delta_y;
double deltaday_min;
double deltaday_max;
double deltaMag_min;
double deltaMag_max;


lcsim_card(){};
lcsim_card(string &card_name);
bool Read(string &card_name);
};




/********************* lcsim **********************/



class Sim {
private:
bool Is_Swarp;
vector<unsigned int> ccds;
string GeoRefFilter;

lcsim_card Sim_Data;

HostList Host_Gal;
GtransfoRef Georef_Pix2RaDec;
GtransfoRef Georef_RaDec2Pix;
double Georef_ZP;
double FirstMJD, LastMJD;

FakeSNIaList snlist;
Image_Manager  allsource;

ReducedImageRef GeoRef;
string MasterName; 	
void Construct_Sn_Random();
void Construct_Sn_Damier();
void Construct_SnList_WHost();
void Construct_SnList_AdaptedToHost();
void Randomize_fakesnList(int code=1023);
bool Is_Swarp_CCD(const unsigned int ccd);
void Get_Method_from_Datacard();
void SetGeoRef(const ReducedImage* image); // this is not  a constructor: image should be create with new before this function
public:
Sim();
Sim(const string &FileName);
string Master()const {return( MasterName);}
string GeoRefName()const {return( (*GeoRef).Name());}

//void DoIt();
void MakeGeoRef(const string &subname , bool overwrite = false);
void MakeSnList(int code = 2047);
void MakeHostList();
void MakeHostListFromSECat();
void MakeHostListFromPhotoZCat();
void SaveList(const string &Name){snlist.write(Name);}// save list as FakeSNIaList
void WriteShortList(const string &Name,const double  MJDate, const string &filter,SaltModel &sn_model, const string magsys="AB"); // save list as a FakeObjectList
void PrintSim(ostream& s = cout);
void MakeSubfiles();
void SetDataCard( string &DataName){Sim_Data.Read(DataName);}
// members to acces an old simulation
void ReadSnList(const string &Name){snlist.read(Name);}
void SetGeoRef(const string &Name); // must be an existing reduced image with catalog... Use with care
void KeepInsideSN();
};

#endif 
