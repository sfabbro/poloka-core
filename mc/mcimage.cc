#include <iostream>

#include "mcimage.h"

#include "gtransfo.h"
#include "wcsutils.h"
#include "simsnstar.h"
#include "fitsimage.h"
#include "reducedimage.h"
#include "daophotpsf.h"
#include "daophotutils.h"
#include "agaussian.h"


MCImage::MCImage(string name, const ReducedImageRef source, const GtransfoRef Transfo, SimSNStarList* List, AddMethod a_method  ) : ReducedImage(name)
{
  if (!FileExists(name))
    Create("here");
  Source = source ;
  Transfo_fromref = Transfo ;  
  SNList = List;
  Methode = a_method ;
  
}
  
MCImage::MCImage(string name,const ReducedImageRef source, SimSNStarList* List, AddMethod a_method  ) : ReducedImage(name)
{
  if (!FileExists(name))
    Create("here");
  Source = source ;
  SNList = List;
  Methode = a_method ;
  Transfo_fromref = new GtransfoIdentity();
}
  
MCImage::MCImage(string name, const ReducedImageRef Ref, const ReducedImageRef source, SimSNStarList* List, AddMethod a_method ) : ReducedImage(name)
{
  if (!FileExists(name))
    Create("here");
  Source = source ;

  Gtransfo* Pix2RaDecref=0; 
  FitsHeader headerref(Ref->FitsName()); 
  WCSFromHeader(headerref, Pix2RaDecref);

  Gtransfo* Pix2RaDec=0;  
  FitsHeader header(Source->FitsName());
  WCSFromHeader(header, Pix2RaDec); 
  Frame globalframe = Source->UsablePart();
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.001,globalframe ); 

  // permet d'aller de xref,yref -> ra,dec -> x,y
  Transfo_fromref = GtransfoCompose(RaDec2Pix,Pix2RaDecref );
   

  SNList = List;
  Methode = a_method ;
}
  

SimSNStarList*  MCImage::PrepareList()
{
  //SNList doit imperativement etre une  SimSNStarList.
  SimSNStarList *lsntf = new SimSNStarList;
  SNList->CopyTo(*lsntf); // si c'est une SimSNWModelStarList 
  // ca plante.
  lsntf->ApplyTransfo(*Transfo_fromref);
  lsntf->NewFlux(Source->ZP());
  return(lsntf);
}
  

bool MCImage::MakeFits()
{
  
  FitsImage *psat = NULL ;
  if (Source->HasSatur())
    {
      cerr << "Opening : " << Source->FitsSaturName() << endl ;
      psat = new FitsImage(Source->FitsSaturName());
    }
  double saturlev =  Source->Saturation();
  bool print_debug = true ;


  {cerr << "Opening : " << Source->FitsName() << endl ;
   FitsImage img(Source->FitsName());
 switch (Methode)
   { 
    case WModel :
      {
	// pas de verification du type du pointer
 	SimSNWModelStarList *lwmodel = (SimSNWModelStarList *) SNList;
	// si ce n'est pas le cas, on est mal.
	cerr <<"Adding FakeSupernova WModel" << endl ;
	lwmodel->AddWModelToImage(img, Transfo_fromref, psat,saturlev,print_debug);
	break ;
      }
   case WGaussian :
     {
       cerr <<"Adding FakeSupernova WGaussian" << endl ;
       SimSNStarList* liste1 = PrepareList();
       if (liste1)
	 {
	   liste1->AddWGaussianToImage(Source->Seeing(),Source->Seeing(),0.,
				     img, psat, saturlev);
	   delete liste1 ;
	 }
       else
	 return(false);
       break ;
     }
   case WDaoPsf :
     {
       
       cerr <<"Adding FakeSupernova WDaoPsf" << endl ;
       SimSNStarList* liste = PrepareList();
       string psfName = Source->ImagePsfName();
       if (!FileExists(psfName.c_str()))
	 {
	   MakeDaoPsf(*Source); 
	   if (!FileExists(psfName.c_str()))
	     {
	       cerr << "Failing to obtain : " << psfName << endl ;
	       return(false);
	     }
	 }
       if (liste)
	 {	 
	   DaoPsf daopsf(Source->ImagePsfName());
	   AddListWDaoPsfToImage(daopsf, (BaseStarList*) liste, img, psat, saturlev);
	 }
       else
	 return(false);
	 
       break ;
     }
    default:
      {
	cerr << "Unknown Method " << Methode  << " for Adding Fake " << endl ;
	break;
      }
   }

  {cerr << "Saving : " << FitsName() << endl ;
  FitsImage imgwsn(FitsName(),img,img); //imgwsn.SetWriteAsFloat();
  }
  }
 // la fits image est sauvee.
  if (psat) 
    {
      {cerr << "Saving : " << FitsSaturName() << endl ;
      FitsImage outsat(FitsSaturName(), *psat, *psat);
      }
      delete psat ;    
    } 
  return(true);
}



//on recopie le catalogue de la source.
bool  MCImage::MakeCatalog() 
{
   string fileName = ReducedImage::CatalogName();
   // si deja fait
  if (FileExists(fileName)) return true;
  SEStarList inList(Source->CatalogName());
  if (FileExists(FitsSaturName()))
    {
      FitsImage sat(FitsSaturName());
      FlagSatFromImage(inList,sat);
    }
  inList.write(fileName);
  return true;
}
  
// recopie, satur est separe de weight de ttes facons
bool MCImage::MakeWeight() 
{ 
  string fileName = FitsWeightName();
  if (FileExists(fileName)) return true;
  if (!Source->HasWeight() && ! Source->MakeWeight())
    {
      cerr <<" cannot make weight for " << Name() << " using " 
	   << Source->Name() << endl;
      return false;
    }
    
  MakeRelativeLink(Source->FitsWeightName().c_str(),
		   FitsWeightName().c_str());
  return(true);
}
// recopie car les weight ne sont pas changes.
bool MCImage::MakeBad() 
{ 
  string fileName = FitsBadName();
  if (FileExists(fileName)) return true;
  if (!Source->HasBad() && ! Source->MakeBad())
    {
      cerr <<" cannot make weight for " << Name() << " using " 
	   << Source->Name() << endl;
      return false;
    }
    
  MakeRelativeLink(Source->FitsBadName().c_str(),
		   FitsBadName().c_str());
  return(true);
}
// recopie car les weight ne sont pas changes.

bool MCImage::MakeSatur() 
{
  if (FileExists(FitsSaturName())) return true;
  else
    return(MakeFits());
}

