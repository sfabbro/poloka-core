#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <minuit.h>
#include <algorithm>
#include <ctime>
#include <cmath>

#include "ccd.h"
#include "sestar.h"
#include "reducedimage.h"
#include "fitsimage.h"
#include "toadscards.h"
#include "hostabsorption.h"
#include "cosmology.h"
#include "instrument.h"
#include "saltmodel.h"
#include "minuit.h"
#include "lcsim.h"
#include "myrdm.h"
#include "wcsutils.h"
#include "fileutils.h"
#include "pathutils.h"
#include "lambert_dust.h"
#include "imagemanager.h"
#define c 299792.458

#ifndef  WAVELENGTH_STEP
/* 10 Angst */
#define WAVELENGTH_STEP 10 
#endif
#define max_dist  0.001

static double mpoly(double x,double *a,int degree = 4)
{

double res = 0.0;
for(int i=degree ; i>=0 ; i--) res = res*x + a[i];


return ( res);
}
static string IntToString(int nombre,int min_char = 0 )
{
	string res=""; 

	if( min_char>0)
		{
		int i;
		if(nombre != 0) i = min_char - int(log10(double(abs(nombre)))+1.0);
		else i = min_char - 1 ;
		while(i>0){res = '0' + res; i--;}
 		}
	if(nombre < 0){nombre= 0-nombre;res = '-'+res;}	
	char a = nombre%10 +'0';
	if (nombre/10 !=0) res = res + IntToString(nombre/10,0);
	res=res+a;
	return (res);
}

static string RefByFilter( const string filter, const string field ,const string card)
{

DataCards Data(card);
string out="";
string key= field +"_MASTER_";

if ( "i" == filter) key =key+"I";
	else if ("r" == filter ) key =key+"R";
	else if ("g" == filter ) key =key+"G";
	else if ("z" == filter ) key =key+"Z";
	else {cerr << "filter not recognized" << endl; exit(1);}
out=Data.SParam(key);	 
return (out);
}
/* random value 
mode 0	: constante return A
mode 1 : flat dsitribution betwen A and B
mode 2 : gaussian distribution sigma = A, mean  = B
*/
static double RangedRndm(double A , double B, int mode=1 )
{
switch (mode)
{
	case 1:
		if ( A == B) return A; 
		else return B +(A-B)*my_rndm();
		break;
	case 0:
		return A;
		break;
	case 2:
		return (A*NormalGaussRand()+B);
		break;

	default : {cerr << "unknown random mode " << endl; exit(1);}
}
}

static double RangedRndm(to_random &in){return RangedRndm(in.A,in.B,in.Mode);}




//Permet de recuperer le filtre en minuscule
static string Band (const string &ImageName)
{
	string band;
	ReducedImage current(ImageName);
	band=current.Band();
	if (band.size()==1)
	{
		if(band[0]> 64 && band[0]< 91) band[0]=band[0]+32;
	}
	else {cerr << "Band error" << endl; exit(1);} 

return (band);
}



 

static string first_word(const char *line, const char sep = ' ')
{
  int start;
  int end_line = strlen(line);

  while ( line[start] == sep && start < end_line) start++;
  int end = start;
  while ( line[end] != sep &&  end < end_line) end++;
  return string(line).substr(start,end);
}


static bool read_random_data(to_random   &out,const string &tag, DataCards &data)
{
const string tag_a = tag+"_A";
const string tag_b = tag+"_B";
const string tag_mode = tag+"_MODE";
bool res = true;
if (!(data.HasKey(tag_a))){cerr << "key "<< tag_a << " not present" << endl; res = false ;}
if (!(data.HasKey(tag_b))){cerr << "key "<< tag_b << " not present" << endl; res = false ;}
if (!(data.HasKey(tag_mode))){cerr << "key "<< tag_mode << " not present" << endl; res = false ;}

out.A=data.DParam(tag_a);
out.B=data.DParam(tag_b);
out.Mode=data.IParam(tag_mode);
return res;
}




static bool read_poly(double   *out,const string &tag,  DataCards &data, int degree = 4/* 4+1 factor*/)
{
bool res = true;
for(int i = 0 ; i<=degree; i++)
	{
	string n_tag = tag+IntToString(i);
	if (!data.HasKey(n_tag)){cerr << "key "<< n_tag << " not present" << endl; res = false ;}
	out[i]=data.DParam(n_tag);
	}
return res;
}




lcsim_card::lcsim_card(string &card_name)
{
Read(card_name);
}


bool 
lcsim_card::Read(string &DataCardFile)
{
bool res= true;
if (!FileExists(DataCardFile))
    {
      cerr << " DatSim::Read Cannot open  FileName " << endl;
      return(false);
    }
  Card_Name= DataCardFile;
  DataCards Data(Card_Name);

  Method = Random;
  string line = Data.SParam("GENERATION_METHOD");
  if (strstr(line.c_str(),"FAKES_INHOST"))Method = InHost;
  else if (strstr(line.c_str(),"FAKES_ADAPTED"))Method = AdaptedToHost; 
  else if (strstr(line.c_str(),"FAKES_RANDOM"))Method = Random; 
  else if (strstr(line.c_str(),"FAKES_ONGRID"))Method = Damier ;
  else { cerr << "wrong GENERATION_METHOD" << endl ; res = false;}
  
res = res && read_random_data( H_ext,"H_EXTINCTION",Data);
res = res && read_random_data( H_Rv,"H_RV",Data);
res = res && read_random_data( Stretch,"STRETCH",Data);
res = res && read_random_data( Redshift,"REDSHIFT",Data);
res = res && read_random_data( Dispersion,"DISPERSION",Data);
res = res && read_random_data( Color,"COLOR",Data);
res = res && read_random_data( Alpha,"ALPHA",Data);
res = res && read_random_data( Beta,"BETA",Data);
res = res && read_random_data( MW_Rv,"MW_RV",Data);
res = res && read_poly(imag_to_red_poly,"IMAG_TO_RED_A",Data);

if(Data.HasKey("K_MIN")) k_min =Data.DParam("K_MIN");
else {cerr << "key K_MIN not present" << endl; res = false ;}

if(Data.HasKey("K_MAX")) k_max =Data.DParam("K_MAX");
else {cerr << "key K_MAX not present" << endl; res = false ;}


if(Data.HasKey("DELTADAY_MIN")) deltaday_min =Data.DParam("DELTADAY_MIN");
else {cerr << "key DELTADAY_MIN not present" << endl; res = false ;}

if(Data.HasKey("DELTADAY_MAX")) deltaday_max =Data.DParam("DELTADAY_MAX");
else {cerr << "key DELTADAY_MAX not present" << endl; res = false ;}


if(Data.HasKey("MIN_DELTA_MAG")) deltaMag_min =Data.DParam("MIN_DELTA_MAG");
else {cerr << "key  MIN_DELTA_MAG not present" << endl; res = false ;}

if(Data.HasKey("MAX_DELTA_MAG")) deltaMag_max =Data.DParam("MAX_DELTA_MAG");
else {cerr << "key MAX_DELTA_MAG not present" << endl; res = false ;}


if(Data.HasKey("FAKE_DELTA_X")) delta_x =Data.IParam("FAKE_DELTA_X");
else {cerr << "key FAKE_DELTA_X not present" << endl; res = false ;}

if(Data.HasKey("FAKE_DELTA_Y")) delta_y =Data.IParam("FAKE_DELTA_Y");
else {cerr << "key FAKE_DELTA_Y not present" << endl; res = false ;}

if(Data.HasKey("INSTRUMENT")  )      instrument = Data.SParam("INSTRUMENT");
else {cerr << "key INSTRUMENT not present" << endl; res = false ;}

if(Data.HasKey("OMEGAM")  )      omegam =Data.DParam("OMEGAM");
else {cerr << "key OMEGAM not present" << endl; res = false ;}

if(Data.HasKey("OMEGAX")   )     omegax =Data.DParam("OMEGAX");
else {cerr << "key OMEGAX not present" << endl; res = false ;}

if(Data.HasKey("W1"))       w1 =Data.DParam("W1");
else {cerr << "key W1 not present" << endl; res = false ;}

if(Data.HasKey("W0"))       w0 =Data.DParam("W0");
else {cerr << "key W0 not present" << endl; res = false ;}

if(Data.HasKey("H0"))       H0 =Data.DParam("H0");
else {cerr << "key H0 not present" << endl; res = false ;}

if(Data.HasKey("NUMBER_OF_FAKES"))      NumberOfFakes=Data.IParam("NUMBER_OF_FAKES");
else {cerr << "key NUMBER_OF_FAKES not present" << endl; res = false ;}



  line = Data.SParam("HOSTMODE");
  if (strstr(line.c_str(),"PHOTO_Z")) Host_Mode= From_photoz;
  else if (strstr(line.c_str(),"SE_CAT"))Host_Mode = From_SE; 
  else { cerr << "wrong HOSTMODE" << endl ; res = false;}


return res;

};


/**********************************  lcsim  *****************************************/



Sim::Sim()
{
	MasterName="";
	FirstMJD=-1.0;
	LastMJD=-1.0;
	string temp = DefaultDatacards();
	Sim_Data.Read(temp);	
	Is_Swarp=false;
}
Sim::Sim(const string &FileName)
{
	Is_Swarp=false;
	string GeoRefName="";
	MasterName="";
	FirstMJD=-1.0;
	LastMJD=-1.0;
	string temp = DefaultDatacards();
	Sim_Data.Read(temp);

	
	FILE *file = fopen(FileName.c_str(),"r");
	if (!file)
	{
		cerr << " cannot open \"" << FileName << "\"" << endl;
		exit(1);
	return;
	}
	char line[512];
	bool inNew = false;
	bool isMaster = false;
	bool isGeoRef = false;
	while (fgets(line,512,file))
	{
		if (line[0] == '#') continue;
		if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
		
		if (strstr(line,"MASTER"))
			{inNew = false; isMaster = true; isGeoRef = false; continue;}
			
		if (strstr(line,"NEW"))
			{inNew = true; isMaster = false; isGeoRef = false; continue;}
			
		if (strstr(line,"GEOREF"))
			{inNew = false; isMaster = false; isGeoRef = true; continue;}
		if (strstr(line,"SWARP"))
		{
		int ccd;
		int i =0;
		string temp = string(line + strcspn(line,"1234567890"));
		bool all = false;
		while(sscanf(temp.c_str(),"%d", &ccd) == 1 && i <100)
			{ 	
			i++;
			temp = string(temp.c_str() + strcspn(temp.c_str()," \t"));
			temp = string(temp.c_str() + strcspn(temp.c_str(),"1234567890"));
			if ( ccd == 99 and ccds.size()==0)
				{ ccds.push_back(ccd); all= true; continue;}
			else if ( ccd >= 0 && ccd < 36 && all==false) ccds.push_back(ccd);
			else	{
				cerr << "SWARP syntax error : unrecognized ccd " << endl;
				exit(1);
				}
			
			}
		
			
		if ( ccds.size() == 0 ) {cerr << "SWARP syntax error : no ccd " << endl;exit(1);}
		Is_Swarp = true;
		continue;
		}
 if (!inNew && !isMaster && !isGeoRef)
	{
	cerr << " ERROR : unrecognised syntax in file " << FileName << endl ;
	cerr << line << endl ;
	cerr <<" stop here " << endl;
	exit(1);
	 return;
	}
 char *start_line = line+strspn(line," "); // skip spaces
 if (start_line[strlen(start_line)-1] == '\n' )
	{
	 start_line[strlen(start_line)-1] = '\0' ;
	}
 if (strlen(start_line) == 0) continue; // skip blank lines
 
 string currentName = first_word(start_line);
 RemovePattern(currentName," ");
 ReducedImage current(currentName);
 if (!current.IsValid())
	{
	 if (isMaster) 
	 { 
	 cerr << " cannot find Master Image : Abort" << endl;
	 exit(1);
	 }
	 else if (inNew)
	 { 
	cerr << " cannot find DbImage : \"" << currentName << "\"" << endl;
 	cerr << " Image skipped" << endl;
	continue;
	}
	else
	{
	cerr << " cannot find GeoRef Image : Abort" << endl;
	 exit(1);
	}

	}
 if (isMaster)
	{
	 if (MasterName != "")
	 {
	 	cerr << "more than one Master : abort" << endl;
		exit(1);
	 }
	 MasterName=currentName;
	 continue;
	}
 if (inNew)
	{
	double date;
	one_image current(currentName);
	if(current.Filtre()=="u") 
	{
	cout << " u filter!!!  skip image " << currentName << endl;
	continue;
	}
		date = current.Date();
		if (date < FirstMJD || FirstMJD<0.0) FirstMJD=date;
		if (date > LastMJD || LastMJD<0.0) LastMJD=date;
	

	allsource.push_back(current);
	cout << "Adding " << current.Name() << " to newlist" ;
	cout << " with date: " << date << endl;
	//cout << FirstMJD << endl;
	continue;
	}
  if (isGeoRef)
	{
	if (GeoRefName=="")
		{
		std::string remainder = string(start_line + strcspn(start_line," \t"));
	  	RemovePattern(remainder," ");
		if (remainder.size() == 0)
		{GeoRefName=currentName;continue;}
		
		char refName[128];
	      	int imin,imax,jmin,jmax;
		if (sscanf(remainder.c_str(),"%[^=]=[%d:%d,%d:%d]",refName, &imin,&imax,&jmin,&jmax) != 5)
		{
		  std::cerr << " can't decode end of line : " 
			    << remainder << std::endl
			    << line << std::endl
			    << " giving up " << std::endl;
		  exit(1);
		}
		SubImage ref(refName, currentName, Frame(imin,jmin,imax,jmax));
	      	ref.Execute(DoFits | DoWeight | DoSatur | DoCatalog );
	      	currentName = refName; // for AllInputImages
	      	GeoRefName = currentName;
		continue;
		}
	else { cerr << " more than on Georef in Simfile: abort" << endl; exit(1);}
		

 	}
	}
 if (allsource.size() == 0 )
 {
 cerr << " ERROR : Did not find any NEW image name in " << FileName << " : stop here " << endl;
 exit(1);
 }
 
 if (MasterName=="" && GeoRefName=="")
 {
 cerr << " ERROR : Did not find any Master or Georef image name in " << FileName << " : stop here " << endl;
 exit(1);
 }
 if (MasterName!="" && GeoRefName!="")
 {
 cerr << " ERROR : cant handle  Master AND GeoRef : stop here " << endl;
 exit(1);
 }
 
 fclose(file);
 
 if (MasterName!="")
 {
 	cout << MasterName << " as Master" << endl;
 	cout << " Read " << FileName << " successfully " << endl;
 	GeoRefFilter=Band(MasterName);
 }
 else SetGeoRef(GeoRefName);
 allsource.sort();
};

void
Sim::SetGeoRef(const string &Name)
{

ReducedImage *current = new ReducedImage(Name);
if (GeoRef) {cout << "Warning : GeoRef Already set ... GeoRef changed!"<< endl;}
 if (!current->IsValid())
	{
	 cerr << " cannot find GeoRef Image : Abort" << endl;
	 exit(1);
	 }

SetGeoRef(current);


 
}

void
Sim::SetGeoRef(const ReducedImage* image)
{
 GeoRef = image ;
 cout << "GeoRef set to " << GeoRef->Name()<< endl;
 GeoRefFilter=Band(GeoRef->Name());

 FitsHeader head(GeoRef->FitsName());
 Frame frame(head, WholeSizeFrame);
 Gtransfo *inPix2RaDec;
 if (!WCSFromHeader(head, inPix2RaDec))
 {
 cerr << " ERROR : cannot handle a large reference without a WCS " 
	 << endl;
 exit(1);
 }
 Gtransfo *inRaDec2Pix = inPix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, frame);
 Georef_RaDec2Pix=inRaDec2Pix;
 Georef_Pix2RaDec=inPix2RaDec;
 Georef_ZP = GeoRef->AnyZeroPoint();
if (Sim_Data.Method==InHost ||  Sim_Data.Method==AdaptedToHost) MakeHostList();
}

void
Sim::MakeGeoRef(const string &subname, bool overwrite  )
{
 {
 ReducedImage current(subname);
 if (current.IsValid() &&  !overwrite)
 { cerr << "subimage allready exists : abort"<< endl; exit(1);}
 }
 if (MasterName==""){cerr << "  No MasterName set : abort" << endl; exit(1);}
 
 
 /* first compute the union of Frame(largeImage)*Frame(anyOtherImage) */
 ReducedImage large(MasterName);
 FitsHeader largeHead(large.FitsName());
 Frame largeFrame(largeHead, WholeSizeFrame);
 Gtransfo *largePix2RaDec;
 if (!WCSFromHeader(largeHead, largePix2RaDec))
 {
 cerr << " ERROR : cannot handle a large reference without a WCS " 
	 << endl;
 exit(1);
 }
 Gtransfo *largeRaDec2Pix = 
 largePix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, largeFrame);
 delete largePix2RaDec;

 Frame toExtract, FrameInter; // =(0,0,0,0)
 for (ManagerIterator i = allsource.begin() ; i != allsource.end(); )
 {
 const string &currentName = i->Name();
 if (currentName == MasterName) {++i;continue;}
 ReducedImage current(currentName);
 FitsHeader currentHead(current.FitsName());
 CountedRef<Gtransfo> current2Large;
 if (HasLinWCS(currentHead))
	{
	 Gtransfo *currentPix2RaDec;
	 WCSFromHeader(currentHead, currentPix2RaDec);
	 current2Large = GtransfoCompose(largeRaDec2Pix, currentPix2RaDec);
	 delete currentPix2RaDec;
	}
 else // go the hard way
	{
	 CountedRef<Gtransfo> large2Current;
	 if (ImageListMatch(current, large, current2Large, large2Current))
	 {
	 cerr << " could not match " << current.Name() << " and " << large.Name() << endl;
	 cerr << " we forget " << current.Name() << endl;
	 i = allsource.erase(i);
	 //allsource.remove(current.Name());
	 continue;
	 }
	}
 cout << " transfo between " << large.Name() << ' ' << " and " 
	 << current.Name() << endl 
	 << *current2Large;
 Frame currentFrame(currentHead, WholeSizeFrame); // 
 currentFrame = ApplyTransfo(currentFrame,*current2Large,SmallFrame);
 /* if we are at the first image in the loop we have to initialize
	 "toExtract", else we just "increment" it so that we have the union
	 of all input images at the end */
 if (toExtract.Area() == 0.)
	{toExtract = (largeFrame*currentFrame);
	FrameInter = (largeFrame*currentFrame);}
 else
	{toExtract *= (largeFrame*currentFrame);
	FrameInter *= (largeFrame*currentFrame);}
 cout << " after adding " << current.Name() 
	 << ", extract frame = " << toExtract;
	 
	 cout << " after intersecting " << current.Name() 
	 << ", inter frame = " << FrameInter;
	 
	 
 ++i;
 }
 delete largeRaDec2Pix;
 SubImage *subImage = new SubImage(subname, MasterName, toExtract);
 // to be able to match, "make" image and catalog
 subImage->MakeFits(); // need image to make catalog!
 subImage->MakeWeight(); // need also weight to make catalog!
 subImage->MakeCatalog();

 SetGeoRef(subImage);
 
 }
 
 
 void
 Sim::MakeSnList(int code)
 {
 if (code & 512) switch (Sim_Data.Method)
    { 
    case  Random :
      cout << "Generating Random Fakes list." << endl;
      Construct_Sn_Random();
      break;
    case Damier :
      cout << "Generating Fakes list on grid " << endl;
      Construct_Sn_Damier();
      break;
    case InHost :  
      cout << "Generating Fakes in Host list." << endl;
      Construct_SnList_WHost();
      break;
    case AdaptedToHost :  
      cout << "Generating Fakes in Host list." << endl;
      Construct_SnList_AdaptedToHost();
      break;      
    default:
      cerr << "Unknown Method " << Sim_Data.Method  << " for Fake SN Generation " << endl ;
    }
     Randomize_fakesnList(code);
     cout << "exiting MakeSnList("<< code<< ")" << endl;
}

void Sim::Construct_Sn_Random()
{


	int Xref=GeoRef->XSize();
	int Yref=GeoRef->YSize();

  int compteur=snlist.size() ;
 
  if ( compteur == 0)  {
  cout  << "list empty: create a new list of " << Sim_Data.NumberOfFakes
	<< " fake sn" << endl;
  while (compteur<Sim_Data.NumberOfFakes)
  
	{
   	FakeSNIa *p = new FakeSNIa();
	snlist.push_back(p);
	++compteur;
	}
  }
  
  FakeSNIaIterator it= snlist.begin();
 while (it != snlist.end()) 
	{
	
      // tirage en coordonnees ref.
      double x = Xref * my_rndm();//coord image ref
      double y = Yref * my_rndm();//coord image ref
      double Ra, Dec;
      Georef_Pix2RaDec->apply(x,y,Ra,Dec);
      


			
     if ( 10 < x && x < Xref-10 && 10 < y && y < Yref-10 && snlist.NumberOfNeighbors(Ra,Dec,max_dist)==0)
	{
	  (*it)->Ra()=Ra;
	  (*it)->Dec()=Dec;
	  it++;
	}
    } 
cerr << "snList size : " << snlist.size() << endl;


 }



void Sim::Construct_Sn_Damier()
{

	
	int Xref=GeoRef->XSize();
	int Yref=GeoRef->YSize();
	int x = Sim_Data.delta_x ;
	int y = Sim_Data.delta_y ;
	double Ra,Dec;
	int compteur=snlist.size() ;
	 if ( compteur == 0)
	 {
		cout  << "list empty: create a new list of fake sn" << endl;
		while (y < Yref )
		{

			Georef_Pix2RaDec->apply(x,y,Ra,Dec);
			FakeSNIa *p = new FakeSNIa(Ra,Dec);
	 		snlist.push_back(p);
	
		x += Sim_Data.delta_x ;
		if (x > Xref){x = Sim_Data.delta_x ; y += Sim_Data.delta_y ;}
		}
		cerr << "snList size : " << snlist.size() << endl;
	}
	else
	{
		FakeSNIaIterator it = snlist.begin(); 
		while (y < Yref && it != snlist.end())
		{

			Georef_Pix2RaDec->apply(x,y,Ra,Dec);
			(*it)->Ra()=Ra;
			(*it)->Dec()=Dec;
	 		it++;
	
		x += Sim_Data.delta_x ;
		if (x > Xref){x = Sim_Data.delta_x ; y += Sim_Data.delta_y ;}
		}
	
	}		

}



void Sim::Construct_SnList_WHost()
{

		
	int Xref=GeoRef->XSize();
	int Yref=GeoRef->YSize();


  HostList BG;
  Host_Gal.CopyTo(BG); //coord as on image ref.
  if (BG.size() == 0)
    {
      cerr << "0 Model galaxies selected, no SNe generation" << endl ;
      exit(1);
    }

  // construction de la liste de supernovae ds gal
  // on tire une galaxie au hasard,on l'enleve de la liste
  // (pas 2 sn ds la meme galaxie !) puis on tire
  // 1 sn au hasard
  
  int nbmax =Sim_Data.NumberOfFakes;
  
  if ( snlist.size() != 0)  
  {
  	cerr  	<< "list already existing: abort" << endl;
	exit(1);

  }



  while (nbmax>0 && BG.size()>0)
    {
     	int i = 0;
	double x = -1;
	double y = -1;
	// tirage au sort de la galaxie en coord ref
      	double xref = Xref * my_rndm();
      	double yref = Yref * my_rndm();
      	
	double H_Ra, H_Dec;
	Georef_Pix2RaDec->apply(xref,yref,H_Ra,H_Dec);
      	Host* galref = BG.FindAndRemoveClosest(H_Ra,H_Dec);

	double  DeltaX, DeltaY ;
	double Ra, Dec;

	while (!( 10 < x && x < Xref-10 && 10 < y && y < Yref-10 && snlist.NumberOfNeighbors(Ra,Dec,max_dist)==0) && i<1000)
      	{
      		i++;
      		// Determination de la position de la Sn coord image ref
      		
		double teta = RangedRndm(0.0, 360.0); //angle avec le grand axe
		double rayon = sqrt(galref->A()*galref->A()*cos(teta)*cos(teta)
				  +galref->B()*galref->B()*sin(teta)*sin(teta));

		DeltaX=rayon*cos(teta+galref->Theta());
		DeltaY=rayon*sin(teta+galref->Theta());
		
		double k = RangedRndm(Sim_Data.k_min,Sim_Data.k_max);
		
		// Allow 20% without host

		if(Sim_Data.Host_Mode==From_SE)
		{
			double testv=my_rndm();
			if (testv<0.2) k=10.;
		}

      		// sn coord on image ref
      		double gal_x,gal_y;
		Georef_RaDec2Pix->apply(galref->Ra(),galref->Dec(),gal_x,gal_y);
      		// sn coord on image ref
      		x = gal_x + k*DeltaX;
      		y =  gal_y+ k*DeltaY;

		Georef_Pix2RaDec->apply(x,y,Ra,Dec);
		
		
	}
	
	if ( i > 999) {delete galref; continue;}
	
	FakeSNIa *p = new FakeSNIa();
	snlist.push_back(p);
	p->Ra()=Ra;
	p->Dec()=Dec;
	p->H_Ra()=galref->Ra();
	p->H_Dec()=galref->Dec();
	p->H_Imag()=galref->I_Mag();

	nbmax--;
	delete galref;
    }
  cerr << snlist.size() << " supernovae generated in host." << endl ;

}





void Sim::Construct_SnList_AdaptedToHost()
{
		
	int Xref=GeoRef->XSize();
	int Yref=GeoRef->YSize();


  HostList BG;
  Host_Gal.CopyTo(BG); //coord as on image ref.
  if (BG.size() == 0)
    {
      cerr << "0 Model galaxies selected, no SNe generation" << endl ;
      exit(1);
    }

  // construction de la liste de supernovae ds gal
  // on tire une galaxie au hasard,on l'enleve de la liste
  // (pas 2 sn ds la meme galaxie !) puis on tire
  // 1 sn au hasard
  
  int nbmax =Sim_Data.NumberOfFakes;
  
  if ( snlist.size() != 0)  
  {
  	cerr  	<< "list already existing: abort" << endl;
	exit(1);

  }



  while (nbmax>0 && BG.size()>0)
    {
     	int i = 0;
	double x = -1;
	double y = -1;
	// tirage au sort de la galaxie en coord ref
      	double xref = Xref * my_rndm();
      	double yref = Yref * my_rndm();
      	
	double H_Ra, H_Dec;
	Georef_Pix2RaDec->apply(xref,yref,H_Ra,H_Dec);
      	Host* galref = BG.FindAndRemoveClosest(H_Ra,H_Dec);

	double  DeltaX, DeltaY ;
	double Ra, Dec;

	while (!( 10 < x && x < Xref-10 && 10 < y && y < Yref-10 && snlist.NumberOfNeighbors(Ra,Dec,max_dist)==0) && i<1000)
      	{
      		i++;
      		// Determination de la position de la Sn coord image ref
      		
		double teta = RangedRndm(0.0, 360.0); //angle avec le grand axe
		double rayon = sqrt(galref->A()*galref->A()*cos(teta)*cos(teta)
				  +galref->B()*galref->B()*sin(teta)*sin(teta));

		DeltaX=rayon*cos(teta+galref->Theta());
		DeltaY=rayon*sin(teta+galref->Theta());
		
		double k = RangedRndm(Sim_Data.k_min,Sim_Data.k_max);
		
		// Allow 20% without host

		if(Sim_Data.Host_Mode==From_SE)
		{
			double testv=my_rndm();
			if (testv<0.2) k=10.;
		}

      		// sn coord on image ref
      		double gal_x,gal_y;
		Georef_RaDec2Pix->apply(galref->Ra(),galref->Dec(),gal_x,gal_y);
      		// sn coord on image ref
      		x = gal_x + k*DeltaX;
      		y =  gal_y+ k*DeltaY;

		Georef_Pix2RaDec->apply(x,y,Ra,Dec);
		
		
	}
	
	if ( i > 999) {delete galref; continue;}
	
	FakeSNIa *p = new FakeSNIa();
	snlist.push_back(p);
	p->Ra()=Ra;
	p->Dec()=Dec;
	p->H_Ra()=galref->Ra();
	p->H_Dec()=galref->Dec();
	p->H_Imag()=galref->I_Mag();
	p->Redshift()=RangedRndm(galref->Z()-galref->Err_Z(),galref->Z()+galref->Err_Z(),1);
	nbmax--;
	delete galref;
    }
  cerr << snlist.size() << " supernovae generated adapted to host." << endl ;

}







 void Sim::Randomize_fakesnList(int code)
{

if(Sim_Data.Method == AdaptedToHost) code = (code & 1019);//mask 000000100 111111011
	
Lambert_Dust my_dust;
for (FakeSNIaIterator it = snlist.begin(); it != snlist.end();++it) 
	{
	
	FakeSNIa *fake=*it;
	if(code & 1)fake->H_Extinction()=RangedRndm(Sim_Data.H_ext);
	if(code & 2)fake->Stretch()=RangedRndm(Sim_Data.Stretch);
	if(code & 4)fake->Redshift()=RangedRndm(Sim_Data.Redshift);
	if(code & 8)fake->D0() = RangedRndm(LastMJD+Sim_Data.deltaday_max,Sim_Data.deltaday_min+FirstMJD,1);
	if(code & 16)fake->Dispersion()=RangedRndm(Sim_Data.Dispersion);
	if(code & 32)fake->Color()=RangedRndm(Sim_Data.Color);
	if(code & 64)fake->Alpha()=RangedRndm(Sim_Data.Alpha);
	if(code & 128)fake->Beta()=RangedRndm(Sim_Data.Beta);
	if(code & 256)fake->H_Rv()=RangedRndm(Sim_Data.H_Rv);
	if(code & 512)fake->MW_Rv()=RangedRndm(Sim_Data.MW_Rv);
	fake->MW_Extinction()=my_dust.Value_From_RaDec(fake->Ra(),fake->Dec());
	
	}

}


bool 
Sim::Is_Swarp_CCD(const unsigned int ccd)
{
if (ccds.size() == 0){ cout << "No ccd to SwarpStack" << endl;  return false;}
for ( int i = 0; i<ccds.size(); i++) if (ccd == ccds.at(i) || ccds.at(i)==99 ) return  true;
return false;
}  


void
Sim::PrintSim(ostream& s )
{
s	<< "@MasterName " << MasterName << endl
        << "@GeoRefName " << (*GeoRef).Name() << endl
	<< "@Datacard " << Sim_Data.Card_Name << endl
	<< "@SnSpectra "<< LocateFileFromCard("SNSPECTRA") <<  endl;
if (Is_Swarp)
	{
	s << "@SWARP " ;
	for (int i =0 ; i < ccds.size(); i++) s << ccds.at(i) ;
	s << endl;


	Image_Manager manager(allsource);
	manager.sort();
	ManagerIterator miter;
	miter = manager.begin();

	while ( miter != manager.end()  )
		{
		int ccd = miter->CCD();
		if (!Is_Swarp_CCD(ccd)) { miter++;continue ;}
		
	
		int intdate;
		string filtre;
		s << "SUB="<< miter->Sub_Path() << endl;

		filtre= miter->Filtre();
		intdate=miter->int_Date();
		miter++;

		
		while( miter != manager.end() && intdate ==  miter->int_Date() && filtre == miter->Filtre() && ccd == miter->CCD()) miter++;
		}	
		
	}
	else
	{
	
	
	Image_Manager manager(allsource);
	ManagerIterator miter=manager.begin();
	while ( miter != manager.end()  )	

	   	{

		
		if(miter->Done =true) {miter++;continue ;}
		
		
		int intdate,ccd ;
		string filtre;
		
		ccd = miter->CCD();		
		filtre= miter->Filtre();
		intdate=miter->int_Date();

		
		s << "SUB="<< miter->Sub_Path() << endl;	


			
			for (ManagerIterator it = miter; it != manager.end();it++)
			{

				
			if (it->Done) continue;	
			if (intdate!= it->int_Date()) continue;
			if (filtre != it->Filtre()  ) continue;
			if (ccd != it->CCD()) continue;
			it->Done=true;
			
			}
		
  	   	}
        }
  
s 	<< endl 
	<< "@FirstMJD "<< FirstMJD << endl
	<< "@LastMJD "<< LastMJD << endl;
s	<< "@Methode ";
 	switch (Sim_Data.Method)
    { 
    case  Random :
      s << "Generating Random Fakes list." << endl;

      break;
    case Damier :
      s << "Generating Fakes list On Grid " << endl;

      break;
    case InHost : 
      s << "Generating Fakes List In Host" << endl;

      break;
    case AdaptedToHost :
      s << "Generating Fakes Adapted to Host"<< endl;

      break;
	default :
	s << endl;
    }

}







void
Sim::MakeSubfiles()
{	
	SaltModel sn_model;
	ofstream sub_log ("sub.log");
	if (Is_Swarp)
	{

	Image_Manager manager(allsource);
	manager.sort();
	ManagerIterator miter;
	miter = manager.begin();

	while ( miter != manager.end()  )
	{

		
		if(miter->Done ==true) {miter++;continue ;}
		int ccd = miter->CCD();		
		if (!Is_Swarp_CCD(ccd)) {miter->Done=true; miter++;continue ;}
		
		int intdate;
		string datestr,filtre,field,dirname,ccdstr;
		
		ccdstr= miter->str_CCD();
		filtre= miter->Filtre();
		field = miter->str_Field();
		intdate=miter->int_Date();
		datestr=miter->str_Date();
		
// construction du directory et du path de substraction		

		dirname= field;
		MKDir(dirname.c_str(),false);
		dirname = dirname +'/'+ datestr;
		MKDir(dirname.c_str(),false);
		dirname = dirname +'/'+ filtre;
		MKDir(dirname.c_str(),false);
		string dbim_dir = dirname+ "/dbim";
		MKDir(dbim_dir.c_str(),false);
		dirname = dirname +"/sub";
		MKDir(dirname.c_str(),false);
		dirname = dirname +"/"+ ccdstr;
		MKDir(dirname.c_str(),false);
		cout << dirname << endl;


		sub_log << "SUB " << dirname << " " << intdate << " " << filtre <<" " << ccd << endl;
	
// construction du subfile		
		
		string subfname = dirname +'/'+"subfile";
		ofstream pr (subfname.c_str());
		pr << "#MJD: " << intdate << endl;
		pr << "# " << dirname << endl;
		

			pr << "NEW SWARP_STACK" << endl;
			string master = RefByFilter(filtre,field,Sim_Data.Card_Name);
			ReducedImage ref(master);
			Toads_ccd::Ccds ccds(master);
			bool H,B,G,D,HG,HD,BG,BD;
			D = !(ccd%9 == 8);
			G = !(ccd%9 == 0);
			H = !(ccd/9 == 0);
			B = !(ccd/9 == 3);
			HG= (H&&G);
			HD= (H&&D);
			BG= (B&&G);
			BD= (B&&D);			
			for (ManagerIterator i= manager.First_with(intdate,filtre); i != manager.end(); ++i) 
			{
			

			if (intdate!= i->int_Date()) break;
			if (filtre != i->Filtre()  ) break;

			if (ccd == i->CCD()) {ccds.AddImage(i->Name()); i->Done=true; continue;}
			if (H  && i->CCD()== ccd- 9){ccds.AddImage(i->Name());continue;}
			if (B  && i->CCD()== ccd+ 9){ccds.AddImage(i->Name());continue;}
			if (D  && i->CCD()== ccd+ 1){ccds.AddImage(i->Name());continue;}
			if (G  && i->CCD()== ccd- 1){ccds.AddImage(i->Name());continue;}
			if (HD && i->CCD()== ccd- 8){ccds.AddImage(i->Name());continue;}
			if (HG && i->CCD()== ccd-10){ccds.AddImage(i->Name());continue;}
			if (BD && i->CCD()== ccd+10){ccds.AddImage(i->Name());continue;}
			if (BG && i->CCD()== ccd+ 8){ccds.AddImage(i->Name());continue;}		

			
			}

			int overlap = 20;
			int ny = 4;
			int nx = 9;
  			const Frame &refFrame = ccds.WholeRefFrame();
  			int xChunk = (int(refFrame.Nx())+(nx-1)*overlap)/nx+1;
  			int yChunk = (int(refFrame.Ny())+(ny-1)*overlap)/ny+1;
  			int xStep = xChunk-overlap;
  			int yStep = yChunk-overlap;

      			int i = ccd % 9;	//    ******************* calcul des coordonéée du ccd
      			int j = ccd / 9; 
      
			double ymax = refFrame.Ny() - j*yStep;
			double ymin = ymax - yChunk;;
			double xmin = i*xStep;
			double xmax = xmin+xChunk;
			
			// *************************************
			


			
			// ********************************
			
			
			Frame subFrame(Point(xmin,ymin),Point(xmax,ymax));
			subFrame *= refFrame; // on sides, it goes beyond large ref limits
			
			cout << "creating ShortList" << endl;
			
			if(FileExists((dbim_dir+'/'+"fakesn.list"))) 
			 {
			 cout << " Short list " << (dbim_dir+'/'+"fakesn.list") << " allready exist" << endl;
			 }
			else 
			{
			cout << "creating ShortList" << endl;
			WriteShortList((dbim_dir+'/'+"fakesn.list"),intdate,filtre,sn_model);
			}

			
			
				
			StringList insideList;
			ccds.SelectInsiders(subFrame, insideList);
			for (StringIterator it = insideList.begin(); it != insideList.end();)
			{
			pr      << *it << endl; 
			sub_log << *it << endl;
			it++;
			}
			
			
			pr << "REF" << endl;
			pr   << ref.Name()<< " ref=" 
			   << '[' << subFrame.xMin << ':' << subFrame.xMax << ','<< subFrame.yMin << ':' << subFrame.yMax << "]" << endl;
			cout << ref.Name()<< " ref=" 
			   << '[' << subFrame.xMin << ':' << subFrame.xMax << ','<< subFrame.yMin << ':' << subFrame.yMax << "]" << endl;
		
		pr.close();

	}	
	
	miter = manager.begin();

	while ( miter != manager.end()  )
		{
		
		int intdate;
		string datestr,filtre,field;
		

		filtre= miter->Filtre();
		field = miter->str_Field();
		intdate=miter->int_Date();
		datestr=miter->str_Date();
		
// construction du directory et du path de substraction		

		string image_list_name = field+'/'+ datestr+'/'+ filtre+ "/dbim/image.list";
		ofstream image_list (image_list_name.c_str());	
		while ( miter != manager.end()  )
		 
			{
			

			if (intdate!= miter->int_Date()) break;
			if (filtre != miter->Filtre()  ) break;
			image_list << miter->Name() << endl;
			miter++;


			
			}
		image_list.close();	
	
	
		}
	
	
	}
	else
	{
	
	
	Image_Manager manager(allsource);
	ManagerIterator miter=manager.begin();
	while ( miter != manager.end()  )	

	{

		
		if(miter->Done =true) {miter++;continue ;}
		
		
		int intdate,ccd ;
		string datestr,filtre,field,dirname,ccdstr;
		ccd = miter->CCD();		
		ccdstr= miter->str_CCD();
		filtre= miter->Filtre();
		field = miter->str_Field();
		intdate=miter->int_Date();
		datestr=miter->str_Date();
		
		
// construction du directory et du path de substraction		

		dirname= field;
		MKDir(dirname.c_str(),false);
		dirname = dirname +'/'+ datestr;
		MKDir(dirname.c_str(),false);
		dirname = dirname +'/'+ filtre;
		MKDir(dirname.c_str(),false);
		string dbim_dir = dirname+ "/dbim";
		MKDir(dbim_dir.c_str(),false);
		dirname = dirname +"/sub";
		MKDir(dirname.c_str(),false);
		dirname = dirname +"/ccd_"+ ccdstr;
		MKDir(dirname.c_str(),false);
		cout << dirname << endl;
// création de la shortlist
		cout << "creating ShortList" << endl;
		sub_log << "SUB " << dirname << " " << intdate << " " << filtre <<" " << ccd <<  endl;
		WriteShortList((dirname+'/'+"fakesn.list"),intdate,filtre,sn_model);		
// construction du subfile		
		
		string subfname = dirname +'/'+"subfile";
		ofstream pr (subfname.c_str());
		pr << "#Date:MJD " << int(FirstMJD) << " + " <<(intdate - int(FirstMJD)) << endl;
		pr << "# " << dirname << endl;
		
				
			pr << "NEW" << endl;
			pr << miter->Name() << endl;
			
			miter->Done=true;
			miter++;
			
			for (ManagerIterator it = miter; it != manager.end();it++)
			{

				
			if (it->Done) continue;	
			if (intdate!= it->int_Date()) continue;
			if (filtre != it->Filtre()  ) continue;
			if (ccd == it->CCD()) { pr << it->Name() << endl;sub_log << it->Name() << endl ; it->Done=true;continue;}	
				
				
				
				
				
			}
			pr << "REF" << endl;
			pr << RefByFilter(filtre,field,Sim_Data.Card_Name)<< "  SubImage" << endl;
		pr.close();

	}	
	}
sub_log.close();	
		
}



void
Sim::WriteShortList(const string &Name,const double  MJDate, const string &filter,SaltModel &sn_model, const string magsys) // save list as a FakeObjectList
{
	cout << "WriteShortList version 0.0" << endl;

	FakeObjectList ShortList;

	     
	GeneralCosmo cosmo_model(Sim_Data.omegam, Sim_Data.omegax, Sim_Data.w0 ,Sim_Data.w1);
	

  double d0 = 10; // pc, distance for absolute magnitude
  //double Mb = -19.46; // absolute B magnitude of a typical Ia

 

	// pour limité le temps de calcul on utilise pas la fonction de la class fakesnia
	//SaltModel sn_model;

	
	LcParam* m_pday = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("DayMax"));
	LcParam* m_stretch = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("Stretch"));
	LcParam* m_dflux = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("Dflux"));
	LcParam* m_ext = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("Tau"));
	LcParam* m_color = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("Color"));
	LcParam* m_rv = const_cast<LcParam *>(sn_model.params.LocateOrAddParam("Rv"));
	
	m_ext->val = 0.0; 
	m_pday->val = .0;
	m_stretch->val = 1;  
	m_color->val =0.0;
	m_rv->val=0.0;
	m_dflux->val = 1;	


  Instrument instrument(Sim_Data.instrument);
  const Filter& instrumental_filter = instrument.EffectiveFilterByBand(filter);
  
  // calcul du point zero
  const Sampled1DFunction &refstar4magsystem = RefSpectrumForMags(magsys); // memory leak here? Non décalaration en static et renvoie adresse si deja existant
  double ref_flux = IntegPhotons(refstar4magsystem,instrumental_filter,0.);
  double zp = 2.5*log10(ref_flux);

	for (FakeSNIaIterator it = snlist.begin(); it != snlist.end();++it) 
	{
	double phase  = (MJDate-(*it)->D0()) / (1+(*it)->Redshift());
	
	if ( phase  > 69.0  || phase < -19.0 ) continue;
	 double flux_at_dl,mag,dl;	
	
	sn_model.SetRedshift((*it)->Redshift());	
  	m_ext->val = (*it)->H_Extinction();
	m_pday->val = (*it)->D0();
	m_stretch->val = (*it)->Stretch();  
	m_color->val = (*it)->Color();		 
	m_rv->val= (*it)->H_Rv();

    	dl = cosmo_model.Dl((*it)->Redshift()); // (en unites de c/H0)
    	dl  = (c/Sim_Data.H0)*dl; // Mpc 
    	dl *= 1.0E6; // pc
    
    	m_dflux->val = pow(d0/dl,2);
	
	//cout << sn_model.Flux(&instrumental_filter,MJDate) << endl;

    /********************************************compute filter with MW absorption begin**************************************/
    
    
  Filter EffectiveFilterWithMilkyWay = instrumental_filter;
  if(fabs((*it)->MW_Extinction())>1.e-10) 
  {
    double data[2];
    data[0] = (*it)->MW_Extinction();
    data[1] = (*it)->MW_Rv(); 
    // mean Rv for the Galaxy (see E.L. Fitzpatrick astro-ph/9809387 , J.A. Cardelli et al., ApJ 425, 245 (1989) )
    Computed1DFunction HostGalacticExtinction(&FluxExtinctionCardelli,data);
    (General1DFunction &) EffectiveFilterWithMilkyWay = Product(WAVELENGTH_STEP,2,&instrumental_filter,&HostGalacticExtinction);   
   } 
    
    /********************************************compute filter with MW absorption end**************************************/	
	
	flux_at_dl = sn_model.Flux(&EffectiveFilterWithMilkyWay,MJDate)/(1+(*it)->Redshift());//*normalisation;
	mag = zp-2.5*log10(flux_at_dl);
 	//cout << flux_at_dl << endl;
	//
	mag = mag + (*it)->Dispersion() -(*it)->Alpha()*((*it)->Stretch()-1) + (*it)->Beta()*(*it)->Color();
	// ajouter une coupure??? 
	if ( mag > 30.0 ) continue;
	//cout << mag << endl;
	FakeObject *p = new FakeObject((*it)->Ra(), (*it)->Dec(),mag);
	ShortList.push_back(p);
		
	}
	ShortList.write(Name);
}

void
Sim::MakeHostList()
{

Host_Gal.clear();
switch(Sim_Data.Host_Mode)
	{
	case From_SE:
		MakeHostListFromSECat();
		break;
	case From_photoz:
		MakeHostListFromPhotoZCat();
		break;
	default:
	cerr << "error on host type" << endl;
	}

cout << Host_Gal.size() << " galaxies selected!" << endl;
Host_Gal.write("host.list");

}





void 
Sim::MakeHostListFromPhotoZCat()
{
cout << "Making Host list with Photo Z catalogue" << endl;

 FitsHeader head(GeoRef->FitsName());
 string field_cat_key = head.KeyVal("OBJECT");
 field_cat_key="HOST_CAT_" + field_cat_key;
 DataCards Data(Sim_Data.Card_Name);
 string cat_name = Data.SParam(field_cat_key);
 cout << "catalogue name : " << cat_name << endl;  
 HostList all_host(cat_name); 
 Frame frame(head, WholeSizeFrame);
 for(HostCIterator it =all_host.begin(); it !=all_host.end() ;it++)
 {
 double x,y;
 Georef_RaDec2Pix->apply((*it)->Ra(),(*it)->Dec(),x,y);
 if(!frame.InFrame(x, y)) continue;
 if(Sim_Data.Method==AdaptedToHost && ((*it)->Z()>Sim_Data.Redshift.B || (*it)->Z()<Sim_Data.Redshift.A)) continue;
 Host_Gal.push_back(*it);

 
 }

}




void
Sim::MakeHostListFromSECat()
{
cout << "Making Host list with SE catalogue" << endl;

  SEStarList liste(GeoRef->CatalogName());


  cout << "Selecting host galaxies " << endl;
  RemoveNonOKObjects(liste); 
  SEStarList lgal;
  // all objects FAINTER than datsim.LimGalMag are galaxies.
  
    // take all galaxies with mag < maxGalMag
  double i_maxGalMag = 26 ;  
  double i_LimGalMag = 23 ;


  
  GalaxyFinder(liste, lgal, i_LimGalMag, Georef_ZP);


   if(Sim_Data.Redshift.Mode !=1 && Sim_Data.Method ==  AdaptedToHost ) 
   { cerr <<" adapted to host method need a redshift mode set to 1 : abort!!" << endl; exit(1);}

	  
   for (SEStarIterator it= lgal.begin(); it != lgal.end();++it ) 
    { 
      SEStar *star = *it;
      double imag =  Georef_ZP-2.5*log10(star->flux);
      if ( imag >  i_maxGalMag) continue;
	


	 double red_min = Sim_Data.Redshift.A;
	 double red_max = Sim_Data.Redshift.B;
	 if ( Sim_Data.Method ==  AdaptedToHost )
         {
	 double tmp_min = mpoly((imag+Sim_Data.deltaMag_min),Sim_Data.imag_to_red_poly);
	 double tmp_max = mpoly((imag+Sim_Data.deltaMag_max),Sim_Data.imag_to_red_poly);

	 red_min = max(Sim_Data.Redshift.A,tmp_min);
	 red_max = min(Sim_Data.Redshift.B,tmp_max);
	 }
	 double redshift,redshift_error;
	 
	 switch ( Sim_Data.Redshift.Mode )
	 	{
		case 1:
		if(red_min>red_max)  continue;	
	 	redshift = (red_min+red_max)/2.0;
	 	redshift_error =(red_max-red_min)/2.0;
	 	break;
		case 0:
		redshift = Sim_Data.Redshift.A;
	 	redshift_error =0.0;
	 	break;
		
		case 2:
		redshift = Sim_Data.Redshift.B;
	 	redshift_error =Sim_Data.Redshift.A;
	 	break;
		default :
		cerr << "Make Host from SECatalogue : unknown Redshift Mode : abort!! " << endl;
		exit(35);
		
		}
		
	 
	 double Ra ,Dec;
	 Georef_Pix2RaDec->apply((*it)->X(),(*it)->Y(),Ra,Dec);
	 Host *p =  new Host(Ra,Dec);
	 p->Z()=redshift;
	 p->Err_Z()=redshift_error;
	 p->Gal_Type()=-1.0;
	 p->A()=(*it)->A();
  	 p->B()=(*it)->B();
	 p->Theta()=(*it)->Gyr_Angle();
	 p->I_Mag()=imag;
	 Host_Gal.push_back(p);
	    
	    
	

    }
}
 
 
void
Sim::KeepInsideSN()
 {
 	int Xref=GeoRef->XSize();
	int Yref=GeoRef->YSize(); 
        cout << "Selecting SN in georef" << endl;
 	for (FakeSNIaIterator it = snlist.begin(); it != snlist.end();) 
	{
		double x,y;
		Georef_RaDec2Pix->apply((*it)->Ra(),(*it)->Dec(),x,y);
		if(10 > x || x > Xref-10 || 10 > y || y < Yref-10 )  it=snlist.erase(it);
		else ++it;
	
	}
 cout <<  snlist.size() <<" SN Selected " << endl;
 
 }

