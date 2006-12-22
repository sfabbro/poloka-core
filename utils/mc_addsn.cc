#include <iostream>
//#include <string>
#include "fakeobject.h"
#include "mcimage.h"
#include "gtransfo.h"
#include "simsnstar.h"
#include "wcsutils.h"
#include "fileutils.h"
#include "pathutils.h" 
#include "fitsimage.h"

 
 
 static void usage(const char * prog)
{
  cerr << " usage : " << endl;
  cerr << prog << " -l <fakelist>  -Z <ZP>  -i <source_image>" << endl ; 

  exit(1);
}

int main(int argc, char **argv)
{  
string list_name ="";
double zero_point = 99.0;
string source_name ="";
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-')  usage(argv[0]);
      switch (arg[1])
	{
	case 'h' : usage(argv[0]);
	case 'i' : i++; if (i >= argc) usage(argv[0]); source_name = argv[i]; break;
	case 'Z' : i++; if (i >= argc) usage(argv[0]);  zero_point= atof(argv[i]); break;
	case 'l' : i++; if (i >= argc) usage(argv[0]); list_name = argv[i]; break;
	default : cerr << " do not understand " << arg << endl; usage(argv[0]);
	}
    }

if(list_name =="" || zero_point == 99.0 || source_name=="") usage(argv[0]);



 string fak = "MC";

 ReducedImageRef source = new ReducedImage(source_name) ; 
 Gtransfo *RaDecToPix;		      
		      

 FitsHeader head(source->FitsName(),RW);
 head.AddOrModKey("MCZP",zero_point);
 Frame largeFrame(head, WholeSizeFrame);
 Gtransfo *Pix2RaDec;
 if (!WCSFromHeader(head, Pix2RaDec))
 {
 	cerr 	<< " ERROR : cannot handle a large reference without a WCS " 
 		<< endl;
 	exit(1);
 }
 RaDecToPix = Pix2RaDec->InverseTransfo(0.01 /* accuracy in pixels*/, largeFrame);
 delete Pix2RaDec;
 	      
		      
		      
		      
// transformation de la liste  en simsnstar avec coordonné sur l'image
	
		      
		      
 if(!FileExists(list_name))
 {
	cerr << list_name << " not find : abort" << endl;
	exit(1);
 }
 FakeObjectList snlist(list_name);
 cout << "snlist.size " << snlist.size()<< endl;
 SimSNStarList SNList;	
 for ( FakeObjectIterator it =snlist.begin();it!=snlist.end();it++)
 {
 	
 	double x,y;
  	double sn_mag= (*it)->Mag() ;
  	double sn_flux = pow(10,0.4*(zero_point-sn_mag));
  	RaDecToPix->apply((*it)->Ra(),(*it)->Dec(),x, y);
	if(largeFrame.InFrame(x,y) || largeFrame.MinDistToEdges(Point(x,y))< 7.0 )
	{
  		SimSNStar *p = new SimSNStar(x,y,sn_flux);
  		p->Mag_SN() = sn_mag ;
  		SNList.push_back(p);
	}
 } 
 cout << "number of added SN : " << SNList.size()<< endl;		      
	/// création de l'image	      
		      
		      
 string nomfake = fak + source_name ;
 cout << "creating " << nomfake << endl;
		      
 MCImage  mcim(nomfake,source, &SNList, WDaoPsf);	     
 bool ok = mcim.Execute(DoFits | DoCatalog | DoSatur | DoWeight);
 if (!ok)
 {
	cerr << "Failing to create " << nomfake << endl ;
	abort();
 }
delete 	RaDecToPix;		
			
}
