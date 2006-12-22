#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>  

#include "fileutils.h"
#include "mcsub.h"
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

//#define N_MCLOOP 3


#include <starlist.cc>

#include "toadscards.h"
static AddMethod ReadMethod()
{
  DataCards Data(DefaultDatacards());
  AddMethod meth = WModel ;
  string line = Data.SParam("ADDITION_METHOD");
  if (strstr(line.c_str(),"WITH_MODEL"))
    meth = WModel ;
  if (strstr(line.c_str(),"WITH_GAUSSIAN"))
    meth = WGaussian ;
  if (strstr(line.c_str(),"WITH_DAOPSF"))
    meth = WDaoPsf ;
  return(meth);
}
static int ReadNLOOP()
{
  DataCards Data(DefaultDatacards());
  int n = Data.IParam("MC_N_LOOP");
  return(n);
}
  
  


// when looping over N_MCLOOP, in loop number i,
// all created images will be named MCi.......
static string ToBuildName(int i)
{
  string s="" ;
  if (i >= 0) 
    {
      char cc[8] ; 
      sprintf(cc,"MC%d",i); 
      s = cc ;
    }
  else //i=-1 quand on ne boucle pas
    {
      s = "MC" ;
    }
  return s;
}

string MCSub::FakeListName() const
{
  string si = ToBuildName(I_MC);
  string name = MCDir()+"/fakes" + si + ".list";
  return(name);
}

string MCSub::MCResultName() const
{
  string si = ToBuildName(I_MC);
  string name = MCDir() + "/mcresult" + si + ".list" ;
  return(name);
}


MCSub::MCSub(Sub const & ASub, const int imc): Sub(ASub), original_globnewname(ASub.GlobalNewName()), original_globsubname(ASub.GlobalSubName()), Original_AllNewNames(ASub.AllNewNames()) 
{
  //ca marche que si geometric reference = reference
  

  cerr << "Entering MCSub " << imc << endl ;
  MKDir(MCDir().c_str());
  cerr << "Directory " << MCDir() << " created." << endl ;
  I_MC = imc ;
  Addition_Method = ReadMethod(); //WModel par defaut
  cerr << "Addition_Method " << Addition_Method << endl ;

  // keeping old names
  Original_GlobSubName() =  GlobalSubName();
  Original_GlobNewName() =  GlobalNewName();

  GlobalNew->dump();

  // keeping old sub
  Original_Sub = GlobalSub ; 
  // keeping old new
  Original_New = GlobalNew ; 

  // MC names
  string fak = ToBuildName(I_MC);
  GlobalSubName() = fak + Original_GlobSubName();
  GlobalNewName() = fak + Original_GlobNewName();


  MakeFakeList();



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

  // on fait la meme chose pour les new_i stacks
  for (unsigned i=0; i< AllNew.size(); ++i)
    {
      NewStack &stacki = AllNew[i]; // c'est le stack concernant la newi
      string namei = fak + stacki.name ;
      // 2 possibilites: 1 stack classique ou 1 stack swarp
      // stack classique: on recupere les images originales
      // on fabrique les FakeImages, on les refoure dans le stack
      // avec la transfo.
      // on fabrique l'image sum.
      if ( stacki.stackType == RegularKind)
	{
	  // on met de cote la sub et la new
	  stacki.original_newStack = stacki.newStack ;
	  stacki.original_sub = stacki.sub ;
	  stacki.name = namei ;

	  ReducedImage *rr = stacki.newStack;
	  ImageSum *sum = dynamic_cast<ImageSum*>(rr);
	  if (sum)// si  c'est bien un *ImageSum
	    {
	      ReducedImageList alignedimages = sum->Components() ;
	      ReducedImageList fake_images ;
	      for (ReducedImageIterator c = alignedimages.begin(); 
		   c != alignedimages.end(); ++c)
		{
		  ReducedImage *ri = *c ;
		  TransformedImage *ti = dynamic_cast<TransformedImage*>(ri);
		  // si ce sont bien des TransformedImage
		  if (ti)
		    {
		      string nomfake = fak + ti->SourceName() ;
		      ReducedImageRef source = new ReducedImage(ti->SourceName()) ;
		      MCImage mcim(nomfake,source, ti->FromRef(), SNList, Addition_Method);
		      bool ok = mcim.Execute(DoFits | DoCatalog | DoSatur | DoWeight);
		      if (!ok)
			{
			  cerr << "Failing to create " << nomfake << endl ;
			  abort();
			}
		      string nomfake_tf = TransformedName(nomfake,GeometricReference->Name());
		      cerr << " Creating " << nomfake_tf  << endl ;
		      TransformedImage *tffake = new TransformedImage(nomfake_tf, mcim, ti->Transfo());
		      fake_images.push_back(tffake);
		    }
		  else
		    {
		      cerr << " The new stack" << i << " component " << ri->Name() << " is not a transformed image: exiting " 
			   << endl ;
		      exit(1);
		    }
		}
	   
	      // on decroche l'ancien stack et on met a la place le nouveau
	      // stack avec les fake.
	      stacki.newStack = new ImageSum(stacki.name, fake_images , RefStack);
	      // on decroche la sub.
	      ImageSubtractionRef vide ;
	      stacki.sub = vide ;
	    }
	  else //stackmc.newStack n'est pas une *ImageSum
	    {
	      cerr << " The new stack " << i << " is not a sum: exiting " << endl ;
	      abort();
	    }
	}
      // should be ok: the components "numbers" such as photratio, backvar, are computed ccording the fake_alignedimage header or catalog which is == to the alignedimage. In the ImageSum constructor, MakeFits will be called on the fake_alignedimages, thus not re-aligning them. MakeCatalog called on components in ImageSum.

      if( stacki.stackType == SwarpKind)
	{
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
	      string nomfake = fak + ri->Name() ;
	      //	      MCImage mcim(nomfake,GeometricReference, ri,SNList ,Addition_Method); 
	      MCImage mcim(nomfake,GeometricReference, ri, SNList,Addition_Method);  
	      mcim.Execute(DoFits| DoCatalog | DoWeight|DoSatur);
	      stackmci.push_back(nomfake);
	    }
	  // on ecrase stacki
	  stacki = stackmci ;
	}

      if (! stacki.original_sub->FitExists() )
	{
	  cerr << "No kernel fit present for : " << stacki.name <<endl ;
	  cerr << "MC mode only working when the subs are re-computed, exiting." << endl ;
	  abort();
	}

    }
  GlobalMCResultName() = "";
  // les resultats obtenus lors de la boucle sur imc seront 
  // mis en append dans GlobalMCResultName.

  if ( I_MC >=  0) // I_MC = -1 si on ne boucle pas
    {
      GlobalMCResultName() = MCDir()+"/global_mcresult.list";
      // lors du 1er passage ds la boucle on l'enleve si il en existe deja un.
      if (I_MC == 0) 
	if ( FileExists(GlobalMCResultName()))
	  remove(GlobalMCResultName().c_str());
    }
}




bool MCSub::CheckTransfoHomogeneity()
{
  ReducedImageList allnewimages(Original_AllNewNames);
   if (allnewimages.size() <= 1) return true ;

  ReducedImageIterator c = allnewimages.begin();
  ReducedImage *r1 = *c;
  
  // skip first image
  c++;
  for ( ; c != allnewimages.end(); ++c)
    {
      ReducedImage *ri = *c ;
      cerr << "Checking transfo from : " 
	   << r1->Name() << " to " 
	   << ri->Name() << " : " << endl ;
      // r1 joue le role de la reference
      ImageGtransfo imTransfo(*r1,*ri);
      imTransfo.FromRef()->dump();

      const GtransfoLin *tf1il = imTransfo.FromRef(); //dynamic cast par la class template de countedref
      if (!tf1il) 
	{
	  cerr << " Guessed transfo is not linear " << endl ;
	}
      else
	{
	  double seuil = 0.001 ;
	  cerr << "Linear transformation : " << endl ;
	  tf1il->dump() ;
	  cerr << "A12, A21 : " << tf1il->A12() << " " << tf1il->A21() << endl ;
	  if ( (tf1il->A12() > seuil)  || (tf1il->A21() > seuil) )
	  {
	    cerr <<"To big a rotation between new images : " << r1->Name() 
		 << " and " << ri->Name() <<  endl ;
	    cerr <<"Addition with Model Star unadapted" << endl ;
	  }
	}
    }
  return true ;
}

// Get the transfo from ref to image 1.

ImageGtransfoRef  MCSub::GetFirstTransfo() 
{
  cerr <<"Getting transfo from ref to 1st image of new stack " << endl ;
  ReducedImageList allnewimages(Original_AllNewNames);
  ReducedImageIterator c = allnewimages.begin();
  ReducedImage *first = *c;
  ImageGtransfoRef imTransfo = new  ImageGtransfo(*GeometricReference,*first);
  return (imTransfo) ;
} 


void MCSub::MakeFakeList()
{
  cout <<" Entering MakeFakeList" << endl ;
  if (FileExists(FakeListName()))
    {
      cout << "Reading fake supernovae list: " << FakeListName() << endl;
      SNList->read(FakeListName());	
    }
  else
    {
      if (Addition_Method == WModel)
	{
	  
	  cout << "Checking transfo Homogeneity." << endl ;
	  bool ok1 = CheckTransfoHomogeneity();
	  if (! ok1)
	    {
	      cerr << "Problem with the transfo between new components.  aborting " << endl;
	      abort() ;
	    }
	  ImageGtransfoRef tf = GetFirstTransfo() ;

	  SEStarList stlref(GeometricReference->CatalogName());
	  SEStarList BellesEtoiles;
	  {FitsImage ref(GeometricReference->FitsName());
	  double satur = ComputeSaturation(ref);
      
	  cout << "Saturation for " << GeometricReference->Name() 
	       << " : " << GeometricReference->Saturation() << endl ;
	  cout << "Computed Saturation for " << GeometricReference->Name() 
	       << " : " << satur << endl ;
	  int nsat = FlagSaturatedStars(stlref, min(satur,GeometricReference->Saturation()));
	  cerr << nsat  << " stars flagged as saturated. " << endl ;

	  ImageSum *sum = (ImageSum*) Original_New;//dynamic cast ds
	  // la class template de countedref.
	  if (!sum)// si stack.newStack est bien un *ImageSum
	    {
	      cout << " The new stack is not a sum" << endl ;
	      abort() ;
	    }
	  ReducedImageList alignedimages = sum->Components() ;
	  bool ok2 = SelectModelStars(alignedimages, stlref,
				      BellesEtoiles, ref);
      
	  if (! ok2)
	    {
	      cout << "Failing selecting model stars for simulation, aborting " << endl;
	      abort() ;
	    }}
	  ForSimWModel forsim_wm(*GeometricReference,BellesEtoiles,tf);
	  SNList = (SimSNStarList*) (forsim_wm.MakeListSNWModel());
	}
      else
	{
	  ForSim forsim(*GeometricReference); 
	  SNList = forsim.MakeListSN();
	}

      cout <<"Writing fake list in " << FakeListName() << endl ;
      SNList->write(FakeListName());
	  
    } 
} 



void
MCSub::MatchDetectionsWithFakes()
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

void MCSub::DoIt()
{
  Sub::DoIt();
  cerr <<"fin de DoIt" << endl ;
  string fakedetname = MCDir()+"fakedet.list";
  MakeDetectionsWithFakes( fakedetname);
  MatchDetectionsWithFakes();
}

void MCProcess(Sub const & ASub)
{
  int N_MCLOOP =  ReadNLOOP();
  if ( N_MCLOOP <=  1 )
    {
      MCSub submc = MCSub(ASub);	  
      submc.DoIt(); 
      char command[512];
      sprintf(command,"l2tup %s",submc.MCResultName().c_str());
      system(command); 
    } 
  else
    {
      cerr << "Looping " <<  N_MCLOOP << " times for MC" << endl ;
      string globalname ;
      for (int i = 0; i < N_MCLOOP ; i++)
	{
	  MCSub submc = MCSub(ASub,i);	  
	  submc.DoIt();
	  globalname = submc.GlobalMCResultName();
	}
      char command[512];
      sprintf(command,"l2tup %s", globalname.c_str());
      system(command); 
    }
}

void MCSub::MakeDetectionsWithFakes(string &ListName)
{
  BaseStarList* positions = (BaseStarList*) SNList ;
  DetectionList detsOnGlobal;
  GlobalSub->RunDetection(detsOnGlobal,positions,ListName,true);


}
