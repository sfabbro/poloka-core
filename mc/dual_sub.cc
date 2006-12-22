#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>  

#include "fileutils.h"
#include "dual_sub.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "simulation.h"
#include "imagesum.h"
#include "transformedimage.h"
#include "gtransfo.h"
#include "mcimage.h"

#include "listmatch.h"
#include "imagematch.h"
#include "imageutils.h"
#include "fakeobject.h"


static string ToBuildName()
{
  string s= "MC" ;
  return s;
}



string DualSub::MCResultName() const
{
  string si = ToBuildName();
  string name = MCDir() + "/mcresult" + si + ".list" ;
  return(name);
}


DualSub::DualSub(Sub const & ASub, string image_path , string list_name): Sub(ASub), original_globnewname(ASub.GlobalNewName()), original_globsubname(ASub.GlobalSubName()), Original_AllNewNames(ASub.AllNewNames()) 
{
  //ca marche que si geometric reference = reference
  

  cerr << "Entering DualSub " <<endl ;
  MKDir(MCDir().c_str());
  cerr << "Directory " << MCDir() << " created." << endl ;
Dual_Image_Path = image_path;

  // keeping old names
  Original_GlobSubName() =  GlobalSubName();
  Original_GlobNewName() =  GlobalNewName();

  GlobalNew->dump();

  // keeping old sub
  Original_Sub = GlobalSub ; 
  // keeping old new
  Original_New = GlobalNew ; 

  // MC names
  string fak = ToBuildName();
  GlobalSubName() = fak + Original_GlobSubName();
  GlobalNewName() = fak + Original_GlobNewName();


  MakeFakeList(list_name);



  // checking if it contains the precious psfmatch
  if (! Original_Sub->FitExists())
    {
      cerr << "kernel fit without fake absent (MC mode only working if sub re-computed), exiting." << endl ;
      abort();
    }
  // on decroche la sub et la new
  ReducedImageRef vide ;
  ImageSubtractionRef subvide ;
  GlobalNew = vide ;
  GlobalSub = subvide ;
  cout << " cerating new stack " << endl;
  // on fait la meme chose pour les new_i stacks
  for (unsigned i=0; i< AllNew.size(); ++i)
    {
      NewStack &stacki = AllNew[i]; // c'est le stack concernant la newi
      string namei = fak + stacki.name ;
      // 1 possibilites: 1 stack swarp


      // should be ok: the components "numbers" such as photratio, backvar, are computed ccording the fake_alignedimage header or catalog which is == to the alignedimage. 
      //In the ImageSum constructor,
      // MakeFits will be called on the fake_alignedimages, thus not re-aligning them. MakeCatalog called on components in ImageSum.


	  NewStack stackmci;
	  stackmci.stackType = SwarpKind ;
	  stackmci.name = namei ;
	  stackmci.original_sub = stacki.sub;
	  stackmci.original_newStack = stacki.newStack;
	  //on recupere les composantes.
	  // on construit les images mc
	  // on les fourre ds stackmci
	  ReducedImageList unalignedImages(stacki);
	  for (ReducedImageIterator i = unalignedImages.begin(); i!= unalignedImages.end(); ++i)
	    {
	      ReducedImageRef ri = *i;
	      string nomfake = Dual_Image_Path + ri->Name() ;
	      // add some check here	      
	      stackmci.push_back(nomfake);
	    }
	  // on ecrase stacki
	  stacki = stackmci ;
	

      if (! stacki.original_sub->FitExists() )
	{
	  cerr << "No kernel fit present for : " << stacki.name <<endl ;
	  cerr << "MC mode only working when the subs are re-computed, exiting." << endl ;
	  abort();
	}

    }
  GlobalMCResultName() = "";


  cout << " end constructor " << endl;
}








void DualSub::MakeFakeList(string list_name)
{
  cout <<" Entering MakeFakeList" << endl ;

SNList = new SimSNStarList();
 double zerop = GlobalSub->ZZZeroP();
 Gtransfo *RaDecToPix;		      
 FitsHeader head(GeometricReference->FitsName());
 Frame largeFrame(head, WholeSizeFrame);
 Gtransfo *Pix2RaDec;
 if (!WCSFromHeader(head, Pix2RaDec))
 {
 	cerr 	<< " ERROR : cannot handle a large reference without a WCS " 
 		<< endl;
 	exit(1);
 }
 RaDecToPix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, largeFrame);
 delete Pix2RaDec;
 	      
		      
		      
		      
// transformation de la liste  en simsnstar avec coordonné sur l'image
	
		      
		      
 if(!FileExists(list_name))
 {
	cerr << "fakesn.list not find : abort" << endl;
	exit(1);
 }
 FakeObjectList fakelist(list_name);
 FakeObjectList inlist;
 cout << "snlist.size " << fakelist.size()<< endl;	
 for ( FakeObjectIterator it =fakelist.begin();it!=fakelist.end();it++)
 {
 	
 	double x,y;
  	double sn_mag= (*it)->Mag() ;
  	double sn_flux = pow(10,0.4*(zerop-sn_mag));
  	RaDecToPix->apply((*it)->Ra(),(*it)->Dec(),x, y);
	if(largeFrame.InFrame(x,y))
	{

		SimSNStar *p = new SimSNStar(x,y,sn_flux);
  		p->Mag_SN() = sn_mag ;
  		SNList->push_back(p);
		inlist.push_back(*it);
	}
 } 
  cout <<"Saving List" << endl ;	  

      SNList->write("mc/fakes.list");
      inlist.write("fakesn.list");
	  
  delete RaDecToPix;  
  cout <<" Ending MakeFakeList" << endl ;
  } 



void
DualSub::MatchDetectionsWithFakes()
{
  cerr <<"Entering MatchDetectionsWithFakes" << endl ;
  BaseStarList * ldet;
  MatchedDetectionList matchedDetections;
  DetectionList detsOnGlobal;
  if (AllNew.size() > 1)
    {
      cerr <<"Reading " << GlobalSub->MatchedCatName() << endl ;
      matchedDetections.read(GlobalSub->MatchedCatName());
      ldet  = (BaseStarList *) &matchedDetections;
    }
  else
    {
      detsOnGlobal.read(GlobalSub->DetectionsName());
      ldet = (BaseStarList *) &detsOnGlobal;
    }
 double maxdist = 1000.0; // pixels
  // get the actual fake flux on Sub.
  //SNList->dump();
  double zerop = GlobalSub->ZZZeroP();
  cerr << "Point zero sub: " <<  zerop << endl ; 
  SNList->NewFlux(zerop);
  //SNList->dump();

  GtransfoIdentity identity;
  BaseStarList * lfake = (BaseStarList*) SNList;
  StarMatchList *matches = 
    ListMatchCollect(*lfake, *ldet,  &identity, maxdist);
  

  cout << " For " << MCResultName() << " we have " << lfake->size() 
       << " generated fakes and recovered " 
       << matches->RecoveredNumber(2.) 
       << " ("   << matches->size() << ")." << endl;
    
  matches->write(MCResultName());
  if ( GlobalMCResultName() != "" )
    {
      if ( ! FileExists(GlobalMCResultName()) )
	matches->write(GlobalMCResultName());
      else
	{
	  ofstream pr(GlobalMCResultName().c_str(), ios::app);
	  matches->write_wnoheader(pr);
	  pr.close();
	}
    }
  delete matches;
}

void DualSub::DoIt()
{
  Sub::DoIt();
  cerr <<"fin de DoIt" << endl ;
  string fakedetname = MCDir()+"fakedet.list";
  MakeDetectionsWithFakes( fakedetname);
  MatchDetectionsWithFakes();
}
 
void ProcessDualSub(Sub const & ASub, string image_path, string list_name)
{

      DualSub submc = DualSub(ASub, image_path , list_name);	  
      submc.DoIt(); 
      char command[512];
      sprintf(command,"l2tup %s",submc.MCResultName().c_str());
      system(command); 

}

void DualSub::MakeDetectionsWithFakes(string &ListName)
{
  BaseStarList* positions = (BaseStarList*) SNList ;
  DetectionList detsOnGlobal;
  GlobalSub->RunDetection(detsOnGlobal,positions,ListName,true);


}
