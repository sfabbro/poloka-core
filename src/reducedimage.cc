#include <iostream>

#include "fitsimage.h"
#include "fileutils.h"
#include "reducedimage.h"
#include "astroutils.h"
#include "frame.h"
#include "fitstoad.h"
#include "sestar.h"
#include "cluster.h"
#include "gtransfo.h"


/* it seems that a constructor cannot call another constructor:
   hence put a routine that both call */

void ReducedImage::init()
{
  actuallyReduced = (DbImage::IsValid() && FileExists(FitsName()) && FileExists(CatalogName()));
}

ReducedImage::ReducedImage(const DbImage &a_DbImage) : DbImage(a_DbImage)
{
  init();
}

ReducedImage::ReducedImage(const string &Name) : DbImage(Name)
{
  init();
} 


bool ReducedImage::MakeFits()
{
  if (FileExists(FitsName())) return true;
  cerr << " No way yet to make " << FitsName() << endl; 
  return false;
}


#include "sextractor_box.h" /* for SEStarListMake */
#include "toadscards.h"
#include "seeing_box.h" /*pour le calcul du seeing lors du catalogue */


/* Pour remplir le carnet d'ordres de SExtractor */
/* si le fond est deja soustrait, on le sauve pas, et on prend = cste
   (SExtractor le recalcule quand meme) sinon, SExtractor le
   calculera, l'utilisera, et on le sauvera.
*/

void 
ReducedImage::FillSExtractorData(ForSExtractor & data, 
				 bool  fond_deja_soustrait, bool sauver_fond,
				 bool use_sigma_header)
{
  
  data.saturation = (int) (SATUR_COEFF*Saturation());
  data.sigma_back  = -1 ;
  if (use_sigma_header)
    {
      if ( SigmaBack() > 1.e-10 )
	{
	  data.sigma_back  = SigmaBack() ;
	  cerr << "Using the background sigma written in header to compute detctions levels: " << data.sigma_back << endl ;
	}
    }
  data.FitsFileName = FitsImageName(Calibrated).c_str();

  if (fond_deja_soustrait)
    {
      data.backmean=0.0;
      data.back_type_manual = true ;
      cout << "Background image taken as constant = 0. " << endl ;
    } 
  if (sauver_fond)
    {
      data.FitsBackName = FitsBackName();
      data.FitsMiniBackName = FitsMiniBackName();
    }
  else
    {
      string vide ;
      data.FitsBackName = vide ;
      data.FitsMiniBackName = vide ;
    }
  if (HasWeight())
    {
      cout << "Weighting from " << FitsWeightName() << endl ;
      data.FitsWeightName = FitsWeightName();
      MakeBad();
      if ( HasBad() )
	data.FitsMaskName = FitsBadName(); // flag image for SExtractor process
      else
	data.FitsMaskName = "" ;
    }
  else 
    {     
      data.FitsWeightName = "";
      if (HasDead())
	{
	  data.FitsMaskName = FitsDeadName(); 
	}
      else
	{
	  cout << "could not open dead image for " <<Name() << endl;
	  cout << "Processing without dead zones image. " << endl ;
	  data.FitsMaskName =  "";
	}
    }


  
}


bool 
ReducedImage::RecoverBack(bool add_to_image) {
  if ( FileExists(FitsBackName()) )
   {
     cerr << "background image existing for image " << FitsName()
	  << endl ;
     if (add_to_image)
       {
	 FitsImage img(FitsName(), RW);
	 FitsImage back(FitsBackName());
	 img += back;
       }
     return true;
   } 
  if ( ! FileExists(FitsMiniBackName()) )
   {
     cerr << "Mini background image NOT existing for image " << FitsName()
	  << endl ;
     return false;
   } 
  FitsImage miniback(FitsMiniBackName());
  int back_meshx = 64 ;
  int back_meshy = 64 ;
  if (miniback.HasKey("SEXBKGSX") && miniback.HasKey("SEXBKGSY"))
    {
      back_meshx = miniback.KeyVal("SEXBKGSX");
      back_meshy = miniback.KeyVal("SEXBKGSY");
    }
  else
    {
      cerr << "No key SEXBKGSX and SEXBKGSY in MiniBack header for " << Name() << ". Can't build Back Image. " << endl ;
      return false;
    }
  Image * imageback = BackFromMiniBack(miniback, XSize(), YSize(), 
				       back_meshx, back_meshy);
  if (add_to_image)
    {
      FitsImage img(FitsName(), RW);
      img += *imageback;
    }
  FitsHeader head(FitsName());
  {FitsImage back(FitsBackName(), head, *imageback);}
  delete imageback ;
  return true ;
  
}

bool
ReducedImage::ReAddBackground_and_ResetKeys()
{
  bool add_to_image = true ;
  if (  RecoverBack(add_to_image) )
    {
      SetSaturation(OriginalSaturation(),"Saturation level corrected from sky subtraction"); 
      SetBackSub(false,"Background was re-added");
      // removing all the keys that were set by the back subtraction
      RemoveBackLevel();
      // removing all the keys concerning the background that were computed by sextractor
      RemoveSESky() ;
      RemoveSESigma();
      return true ;
    }
  else
    {
      cerr << "Failed to Re Add Background to image " << FitsName()
	   << endl ;
      return false;
    }
}



bool
ReducedImage::MakeCatalog(bool redo_from_beg, 
			  bool overwrite, bool savemasksat,
			  bool pas_sub_fond, 
			  bool use_sigma_header)
{
  if (redo_from_beg) overwrite = true ;
   if ( !overwrite )
    { 
      if ( FileExists(ImageCatalogName(SExtractor)))
	{
	  cerr << "catalog already done for image, exiting " << FitsName()
	       << endl ;
	  return true;
	}
    }

   // if the background has already been subtracted,then its value will be taken as 0. (but SExtractor will re-compute the sigma)and it won't in any case be subtracted again.
  bool fond_deja_soustrait = BackSub();
  if (redo_from_beg)
    {
      cerr << "Trying to reprocess image " << Name() << " from beginning" << endl ;
      if (BackSub())
	{
	  cerr <<"Trying to re-add background ...: " ;
	  if (ReAddBackground_and_ResetKeys())
	    {
	      cerr << "success." << endl ;
	      fond_deja_soustrait = BackSub(); //should be false
	      if (HasBack() )
		{
		  cerr <<"Removing Back map for " << Name() << endl ;
		  remove(FitsBackName().c_str());
		}
	      if (HasMiniBack() )
		{
		  cerr <<"Removing MiniBack map for " << Name() << endl ;
		  remove(FitsMiniBackName().c_str());
		}
	    }
	  else
	    cerr << "FAILED." << endl ;
	}
      RemoveSeeing() ;
      if (HasWeight() )
	{
	  cerr <<"Removing weight map for " << Name() << endl ;
	  remove(FitsWeightName().c_str());
	}
      if (HasCosmic() )
	{
	  cerr <<"Removing cosmic map for " << Name() << endl ;
	  remove(FitsCosmicName().c_str());
	}
      if (HasSatellite() )
	{
	  cerr <<"Removing satellite map for " << Name() << endl ;
	  remove(FitsSatelliteName().c_str());
	}
      if (HasBad())  
	{
	  cerr <<"Removing bad map for " << Name() << endl ;
	  remove(FitsBadName().c_str()); 
	}
    }
  // fabrication de la carte de poids
   if (! HasWeight() )
     {
       MakeWeight();
     }

  ForSExtractor data ;
  

   // the background is subtracted if we do want to subtract it AND  if it wasn't subtracted before. 
  bool soustraire_fond = !(fond_deja_soustrait) && !pas_sub_fond;

  // the background map is saved only if it is subtracted.
  bool sauver_fond = soustraire_fond ;
  FillSExtractorData(data,fond_deja_soustrait,sauver_fond, use_sigma_header);
  data.Print();
  double Fond = 0., SigmaFond = 0. ;
  FitsImage *pmasksat = NULL ;
  string nommasksat ;
  if (savemasksat)
    {
      FitsHeader imgheader(FitsName());
      pmasksat = new FitsImage(FitsSaturName(), imgheader);
      pmasksat->AddOrModKey("BITPIX",8);
      pmasksat->EnableWrite(false); // in case something goes wrong
    }
  SEStarList stlse ;
  int status = SEStarListMake(data, stlse, Fond, SigmaFond, pmasksat);
  if (status == 0) {if (pmasksat) delete pmasksat; return 0;}

  // we re-flag for saturation here
  int nsat = FlagSaturatedStars(stlse, data.saturation);
  cerr << nsat << " stars flagged as saturated. " << endl ;
  stlse.write("check_se.list");

  if (pmasksat) 
    { 
    pmasksat->EnableWrite(true); 
    cerr << "Writing Saturated stars pixels map in " 
	 << FitsSaturName() << endl ;
    delete pmasksat;
    }

  SetSESky(Fond, "SExtractor computed background");
  SetSESigma(SigmaFond, "SExtractor computed sigma on background"); 
  // mise a zero starliste s le fond est soustrait
  // et correction de la saturation 
  if (soustraire_fond)
    {
      // mise a zero du fond des etoiles
      SetStarsBackground(stlse,0.);
      // on update les cles
      SetOriginalSaturation(Saturation(),"Original saturation level before sky subtraction"); 
      SetSaturation(Saturation()-Fond,"Saturation level corrected from sky subtraction"); 
      SetBackLevel(0.," SExtractor Background subtracted"); // activates BackSub 
      // subtracting background and saving image
      if(true)
	{
	  cout << "TOADS: Subtracting Image Background" << endl ;
	  FitsImage back(FitsBackName(), RO);
	  FitsImage img(FitsName(), RW);
	  img -= back ;
	  cout << " removing the background image computed by sextractor : " 
	       << endl << FitsBackName() << endl;
	  remove(FitsBackName().c_str());
	}
    }


  if (HasBad())  
    {
      cout << " removing binary mask " << FitsBadName() << endl;
      remove(FitsBadName().c_str()); 
    }

  // Calcul du Seeing
  DatSeeing datsee(data.saturation) ;
  SortieSeeing sortiese;
  SEStarList seestar;
  CalculeSeeingSE(datsee, sortiese, stlse, seestar );
  SetSeeing(sortiese.seeing, "SExtractor computed seeing (pixels, sigma)"); 
  sortiese.Print();

  MakeWeight();

  stlse.write(CatalogName());

  return(status);
}


bool ReducedImage::MakeCatalog() 
{ 
  bool ok = MakeCatalog(/*redo_from_beg=*/ false, 
			/* overwrite = */ false, 
			/*savemasksat= */ true,
			/*pas_sub_fond= */false,
			/*use_sigma_header= */  false);
  MakeCosmic();
  return(ok);
}




bool ReducedImage::IsToadsWeight()
{
  if ( HasWeight() )
    {
      FitsHeader head(FitsWeightName());
      if ( !(head.HasKey("DEADPIXS")) && !(head.HasKey("COSMPIXS"))
	   && !(head.HasKey("FLATPIXS")) && !(head.HasKey("VARPIXS")) 
	   && !(head.HasKey("SATEPIXS")) )
	{
	  return false ;
	} 
    }
  return true ;
}

bool ReducedImage::MakeWeight()
{
  if ( HasWeight() && (! IsToadsWeight()) )
    {
      cerr << "Weight map : " << FitsWeightName() 
	   << " from un-known source, exiting from MakeWeight " << endl ;
      return true ;
    }
  
  FitsImage *pweights = NULL;
  bool destroy_weights = false ; // wether to remove existing weight map or not
  bool adddead = true;
  bool addcosmic = true;
  bool addsatellite = true;
  bool addflat = true;
  bool updatevar  = true;
  double oldvariance = 1 ;
  double newvariance = 1./(SigmaBack()*SigmaBack());
  if (HasWeight())
    {
      if (true)
	// to keep FitsHeader of weight as local, 
	// avoiding conflict opening of the same file
	{
	  FitsHeader head(FitsWeightName());
	  adddead = !(head.HasKey("DEADPIXS") && head.KeyVal("DEADPIXS"));
	  addcosmic = !(head.HasKey("COSMPIXS") && head.KeyVal("COSMPIXS"));
	  addflat = !(head.HasKey("FLATPIXS") && head.KeyVal("FLATPIXS"));
	  addsatellite = !(head.HasKey("SATEPIXS") && head.KeyVal("SATEPIXS"));

	  updatevar = false ;
	  // many different cases:
	  //1)no key VARPIXS: old map, remove and start from beginning
	  if ( !(head.HasKey("VARPIXS")) ) // no key VARPIXS
	    {
	      // on enleve la carte, on la recree.
	      destroy_weights = true ;
	      adddead = true ;
	      addcosmic = true ;
	      addflat = true ;
	      updatevar = true ;
	      oldvariance = 1. ;
	    }
	  //2) key VARPIXS = true, but no value for the inverse of variance INVERVAR specified: idem as above
	  
	  if ( head.HasKey("VARPIXS") && head.KeyVal("VARPIXS") && 
	       ! ( head.HasKey("INVERVAR") ) )
	    {
	      destroy_weights = true ;
	      adddead = true ;
	      addcosmic = true ;
	      addflat = true ;
	      updatevar = true ;
	      oldvariance = 1 ;
	    }
	  //3) key VARPIX at false: update for variance
	  if ( head.HasKey("VARPIXS") && !(head.KeyVal("VARPIXS")) )
	    {
	      updatevar = true ; 
	      oldvariance = 1. ;
	    }
	  
	  //4) key VARPIX at true, value for the inverse of 
	  // variance INVERVAR specified. update only 
	  // if INVERVAR different from the new value
	  if ( head.HasKey("VARPIXS") && head.KeyVal("VARPIXS") && 
	       head.HasKey("INVERVAR") )
	    {
	      oldvariance = head.KeyVal("INVERVAR");
	      if (fabs(oldvariance - newvariance)/newvariance > 1.e-5)
		{
		  updatevar = true ;
		  cerr<< "Old var: " << oldvariance << ", new var:  " 
		      << newvariance << "  for " 
		      << FitsWeightName()   << endl ; 
		}
	    }
	}
      if ( !updatevar && !adddead && !addcosmic && !addflat && !addsatellite) 
	return true;
      if (destroy_weights)
	{
	  cout << " removing " << FitsWeightName() << endl;
	  remove(FitsWeightName().c_str());
	}
      else
	{
	  pweights = new FitsImage(FitsWeightName(), RW);
	  cout << " Updating weight map :" << pweights->FileName() << endl;
	}
    }
  
  if (! HasWeight() ) // wasn't there or was destroy in last loop.
    {
      pweights = new FitsImage(FitsWeightName(), FitsHeader(FitsName())); 
      cout << " Creating  weight map :" << pweights->FileName() << endl;
      // impose that zeros are preserved after FITS write/read
      pweights->PreserveZeros(); 
      *pweights += 1. ; 
    }

  FitsImage &weights = *pweights; // does nothing!


  // check if we have dead, cosmics and flat frames and use them.
  // satur is out of the game because we wish to keep it separate.
  if (updatevar)
    {
      //fill with the inverse of the sky variance
      double s = newvariance / oldvariance ;
      weights *= s ; 
      weights.AddOrModKey("VARPIXS",true);
      weights.AddOrModKey("INVERVAR",newvariance );
      cout << "accounting for sky variance in  " << FitsWeightName() 
	   << " old, new variance^-1: " <<  oldvariance << " " 
	   << newvariance << endl;
      cout << "accounting for sky variance in  " << FitsWeightName() << endl;
    }
  if (adddead && HasDead())
    {
      FitsImage dead(FitsDeadName());
      weights *= (1.- dead);
      weights.AddOrModKey("DEADPIXS",true);
      cout << " zeroing dead pixels in " << FitsWeightName() << endl;
    }
  if (addcosmic && HasCosmic())
    {
      FitsImage cosmic(FitsCosmicName());
      weights *= (1.- cosmic);
      weights.AddOrModKey("COSMPIXS",true);
      cout << " zeroing cosmic pixels in " << FitsWeightName() << endl;
    }

  if (addsatellite && HasSatellite())
    {
      FitsImage satellite(FitsSatelliteName());
      weights *= (1.- satellite);
      weights.AddOrModKey("SATEPIXS",true);
      cout << " zeroing satellite pixels in " << FitsWeightName() << endl;
    }

  // flat contribution
  if (addflat && HasFlat())
    {
      FitsImage flat(FitsFlatName());
      /* We want to account for flat variations, for a flat normalized to 1.,
        not flat values which decribe gain differences due e.g. to 
	different amplifiers */
      FitsHeader head(FitsName());
      Frame illuFrame = TotalIlluRegion(head);
      // sometimes we get flats which are not trimmed
      if ((flat.Nx() != weights.Nx()) || (flat.Ny() != weights.Ny()))
	{
	  cerr << " warning : " << flat.FileName() 
	       << " was seemingly not trimmed !. " << endl
	       << "trimming it  in memory " << endl;
	  flat.Trim(illuFrame);
	}
      int namp = head.KeyVal("TOADNAMP");
      for (int iamp = 1; iamp<= namp; ++iamp)
	{
	  Frame thisAmpFrame = IlluRegion(head,iamp);
	  Pixel flatAverage, sigma;
	  flat.SkyLevel(thisAmpFrame, &flatAverage, &sigma);
	  double factor = 1/flatAverage;
	  int jmin = int(thisAmpFrame.yMin-illuFrame.yMin+0.5);
	  int imin = int(thisAmpFrame.xMin-illuFrame.xMin+0.5);
	  int jmax = int(thisAmpFrame.yMax-illuFrame.yMin+0.5); // beyond last
	  int imax = int(thisAmpFrame.xMax-illuFrame.xMin+0.5);
	  for (int j=jmin; j<jmax; ++j)
	    for (int i=imin; i<imax; ++i)
		flat(i,j) *= factor;
	} 
      weights.MultiplyBySquare(flat) ;
      weights.AddOrModKey("FLATPIXS",true);
      cout << " accounting for flat variations in " << FitsWeightName() << endl;
    }
  // print some statistics
  {
    FitsImage im(FitsName());
    if (im.IsValid())
      cout << Name() << " Image/Weight stats : " 
	   << ImageAndWeightError(im,weights) << endl;
  }
  if (pweights) delete pweights;
  return(true);
}


bool ReducedImage::MakeBad()
{ 
  if (HasWeight())
    {
      cout << " making " << FitsBadName() << endl;
      FitsImage weight(FitsWeightName());

      Image &bad = weight;
      bad.Simplify(0,0,1);// set to 0 what is > 0 and to 1 otherwise
      FitsHeader head(FitsWeightName());      
      FitsImage imBad(FitsBadName(), head, bad);
      imBad.ModKey("BITPIX",8);
      return true ;
    }
  else
    return false ;
}

bool ReducedImage::MakeSatur()
{
  cerr << "  ReducedImage::MakeSatur() should in principle never be called .... " << endl;
// because it is done when making the catalog.. we could however have a rescue routine..
  return false;
}

bool ReducedImage::MakeDead()
{
  cerr << "  ReducedImage::MakeDead() should in principle never be called .... " << " Name : " << Name() << endl;
  // because it is in fact shared between images sharing the same flat. no way do build it.
  return false;
}



#include "fastfinder.h"

void ReducedImage::FlagCosmicsInCatalog(const Image &CosmicImage,
					const double dist)
{
  cout << " Updating catalog" << endl;
  SEStarList stars(CatalogName());

  FastFinder finder(*SE2Base(&stars));
  int nx = CosmicImage.Nx();
  int ny = CosmicImage.Ny();
  int nobj = 0;
  for (int j=0; j<ny; ++j) for (int i=0; i<nx; ++i) 
    if (CosmicImage(i,j) > 0) 
      {
	Point where(i,j);
	const BaseStar *b = finder.FindClosest(where, dist);
	SEStar *cosmic = (SEStar *) b;
	if (cosmic && (cosmic->Distance(where) < dist)) 
	  {
	    nobj++; 
	    cosmic->FlagAsCosmic();
	  }
      }
  cout << " Number of flagged objects: " << nobj << endl;
  stars.write(CatalogName());
}

bool ReducedImage::MakeCosmic()
{
  if (FileExists(FitsCosmicName())) return true;
    clock_t tstart = clock();
  {
    FitsHeader head(FitsName());
    if (head.HasKey("COMBINET"))
      {
	string combinet = head.KeyVal("COMBINET");
	if (strstr(combinet.c_str(),"MEDIAN"))
	  {
	    cout << " image " << head.FileName() 
		 << " was obtained by median combination" << endl;
	    cout << " no cosmic frame needed " << endl;
	    return false;
	  }
      }
    if (head.HasKey("STACKMET")) // this test could be refined
      {
	cout << " image " << head.FileName() 
	     << " is in fact an ImageSum, cosmic free in principle " << endl;
	return false;
      }
  }// close the fits file.
  FitsImage Img(FitsName());
  Pixel Mean, Sigma;
  Img.SkyLevel(&Mean, &Sigma);
  double seeing = Seeing();
    {
      FitsHeader head(FitsName());
      FitsImage cosmic(FitsCosmicName(), head);
      Image &CosmicImage = cosmic;
      Img.Cosmics(Sigma, Mean, seeing, CosmicImage);  
      cosmic.AddOrModKey("BITPIX",8);

      if (HasCatalog()) 
	FlagCosmicsInCatalog(cosmic);
      else
	cerr << "ERROR : " << CatalogName() << "does not exist" << endl;
    }

  // update weights:
    MakeWeight();
  // CPU : 
  clock_t tend = clock();
  cout << " CPU for cosmics + update weights " 
       << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
  return(true);
}



bool ReducedImage::MakeSatellite()
{
  if (FileExists(FitsSatelliteName())) return true;
  clock_t tstart = clock();
  {
    FitsHeader head(FitsName());
    if (head.HasKey("COMBINET"))
      {
	string combinet = head.KeyVal("COMBINET");
	if (strstr(combinet.c_str(),"MEDIAN"))
	  {
	    cout << " image " << head.FileName() 
		 << " was obtained by median combination" << endl;
	    cout << " no satellite frame needed " << endl;
	    return false;
	  }
      }
    if (head.HasKey("STACKMET")) // this test could be refined
      {
	cout << " image " << head.FileName() 
	     << " is in fact an ImageSum, satellite free in principle " << endl;
	return false;
      }
  }// close the fits file.
  FitsImage Img(FitsName());
  {
    ClusterList clustList(*this);
    clustList.Cut();
    Image mask = clustList.Mask();
    FitsHeader header(FitsName());
    FitsImage satellite(FitsSatelliteName(), header, mask);
    satellite.AddOrModKey("NBRSATEL",double(clustList.size()));
    satellite.AddOrModKey("BITPIX",8);
  }
  
  // update weights:
  MakeWeight();
  // CPU : 
  clock_t tend = clock();
  cout << " CPU for satellite + update weights " 
       << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
  return(true);
}





/* The two following routines have to do with reloading saved ReducedImage's
and inheriters. They would disappear if we switched to root for I/O's .
but in fact what we have at the moment just does not work. i.e. most
of the reduced images cannot be put back in the same state as when
they were created just using their name. */

#define TYPE_FILE_NAME "type.name"

string ReducedImage::TypeFileName() const
{
return Dir()+"type.name";
}

bool ReducedImage::SetTypeName(const string &TypeName)
{
string fileName = TypeFileName();
FILE *file = fopen(fileName.c_str(),"w");
if (!file) 
  {
    cerr << " could not open " <<  fileName << endl;
    return false;
  }
fprintf(file,"%s",TypeName.c_str());
fclose(file);
return true;
}
			       

string ReducedImage::TypeName() const
{
char name[256];
string fileName = TypeFileName();
FILE *file = fopen(fileName.c_str(),"r");
if (!file)
  {
    cerr << "cannot open in write mode " << fileName << endl;
    return "NoType";
  }
fscanf(file,"%s",name);
fclose(file);
return string(name);
}


bool ReducedImage::ActuallyReduced() const
{
return actuallyReduced;
}

bool  ReducedImage::Execute(const int ToDo)
{
  bool status = true;
  if (ToDo & DoFits) status &= MakeFits();
  if (ToDo & DoCatalog) status &= MakeCatalog();
  if (ToDo & DoDead) status &= MakeDead();
  if (ToDo & DoSatur) status &= MakeSatur();
  if (ToDo & DoCosmic) status &= MakeCosmic();
  if (ToDo & DoSatellite) status &= MakeSatellite();
  if (ToDo & DoWeight) status &= MakeWeight();
  return status;
}

  

#define UNDEFINED -1

// hidden routines



void 
ReducedImage::remove_key(const char *KeyName, const string &RoutineName)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      if (header.IsValid() && header.HasKey(KeyName) )
	header.RmKey(KeyName);
    }
  else
    {
      cerr << " ReducedImage::"<< RoutineName << " : cannot Read key " << KeyName << " for image " << Name() << " file " << fileName << endl;
    }
}


bool
ReducedImage::has_key(const char *KeyName, const string &RoutineName) const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName);
      if (header.IsValid())
	if ( header.HasKey(KeyName) )
	  return true;
	else
	  return false ;
    }
  cerr << " ReducedImage::"<< RoutineName 
       << " : cannot Read key " << KeyName 
       << " for image " << Name() << " file " 
       << fileName << endl;
    
  return false ;
}




int ReducedImage::read_int_key(const char *KeyName, const string &RoutineName) const
{
  string fileName = FitsName();
    if (FileExists(fileName))
    {
      return int (FitsHeader(fileName).KeyVal(KeyName));
    }
  else
    {
      cerr << " ReducedImage::"<< RoutineName << " : cannot Read key " << KeyName << " for image " << Name() << " file " << fileName << endl;
    }
return UNDEFINED;
}

double ReducedImage::read_double_key(const char *KeyName, const string &RoutineName) const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      return double (FitsHeader(fileName).KeyVal(KeyName));
    }
  else
    {      cerr << " ReducedImage::"<< RoutineName << " : cannot Read key " << KeyName << " for image " << Name() << " file " << fileName << endl;
    }
return UNDEFINED;
}


bool ReducedImage::set_double_key(const double &Value, const char *KeyName, 
                                  const string &RoutineName,
                                  const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (header.IsValid() && header.AddOrModKey(KeyName, Value, Comment));
    }
  else
    {
      cerr << " ReducedImage::" << RoutineName << " : could not write Key " << KeyName << " for file " << fileName << endl;
      return false;
    }
}

string ReducedImage::read_string_key(const char *KeyName, const string &RoutineName) const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      return string(FitsHeader(fileName).KeyVal(KeyName));
    }
  else
    {
      cerr << " ReducedImage::"<< RoutineName << " : cannot read key " << KeyName << " for image " << Name() << " filename " << fileName << endl;
    }
  return "ERROR";
}


bool ReducedImage::set_string_key(const string &Value, const char *KeyName, 
                                  const string &RoutineName,
                                  const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (header.IsValid() && header.AddOrModKey(KeyName, Value.c_str(), Comment));
    }
  else
    {
      cerr << " ReducedImage::" << RoutineName << " : could not write Key " 
	   << KeyName << " for file " << fileName << endl;
      return false;
    }
}

bool ReducedImage::set_int_key(const int &Value, const char *KeyName, 
                                  const string &RoutineName,
                                  const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (header.IsValid() && header.AddOrModKey(KeyName, Value, Comment));
    }
  else
    {
      cerr << " ReducedImage::" << RoutineName << " : could not write Key " 
	   << KeyName << " for file " << fileName << endl;
      return false;
    }
}
      
bool ReducedImage::set_bool_key(const bool &Value, const char *KeyName, 
                                  const string &RoutineName,
                                  const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (header.IsValid() && header.AddOrModKey(KeyName, Value, Comment));
    }
  else
    {
      cerr << " ReducedImage::" << RoutineName << " : could not write Key " 
	   << KeyName << " for file " << fileName << endl;
      return false;
    }
}
      
//**************routines to read/write FitsKeys************
int ReducedImage::XSize() const
{
  return read_int_key("NAXIS1","XSize()");
}

int ReducedImage::YSize() const
{
  return read_int_key("NAXIS2","YSize()");
}

#define READ_ROUTINE(ROUTINE_NAME, VALUETYPE, FITS_KEY_NAME)\
VALUETYPE ReducedImage::ROUTINE_NAME() const\
{\
  return read_##VALUETYPE##_key(FITS_KEY_NAME,#ROUTINE_NAME);\
}

#define  SET_ROUTINE(ROUTINE_NAME, VALUETYPE, FITS_KEY_NAME) \
bool ReducedImage::Set##ROUTINE_NAME(const VALUETYPE &Value, const string Comment)\
{\
  return set_##VALUETYPE##_key(Value,FITS_KEY_NAME,"Set"#ROUTINE_NAME, Comment);\
}

#define  REMOVE_ROUTINE(ROUTINE_NAME, FITS_KEY_NAME) \
void ReducedImage::Remove##ROUTINE_NAME()\
{\
  remove_key(FITS_KEY_NAME,"Remove"#ROUTINE_NAME);\
}

#define  HAS_ROUTINE(ROUTINE_NAME, FITS_KEY_NAME) \
bool ReducedImage::Has##ROUTINE_NAME() const \
{\
  return( has_key(FITS_KEY_NAME,"Has"#ROUTINE_NAME));\
}

#define READ_AND_SET_ROUTINES(ROUTINE_NAME, VALUETYPE, FITS_KEY_NAME) \
        READ_ROUTINE(ROUTINE_NAME, VALUETYPE, FITS_KEY_NAME)\
        SET_ROUTINE(ROUTINE_NAME, VALUETYPE, FITS_KEY_NAME)





READ_AND_SET_ROUTINES(Exposure,double,"TOADEXPO") 
READ_AND_SET_ROUTINES(Epoch,double,"TOADEQUI") 
READ_AND_SET_ROUTINES(PixelSize,double,"TOADPIXS") // could actually use WCS routine PixelSize 
READ_AND_SET_ROUTINES(Gain,double,"TOADGAIN") 
READ_AND_SET_ROUTINES(ReadoutNoise,double,"TOADRDON") 
READ_AND_SET_ROUTINES(Band,string,"TOADBAND") 
READ_AND_SET_ROUTINES(Filter,string,"TOADFILT") 
READ_AND_SET_ROUTINES(Chip,int,"TOADCHIP") 
READ_AND_SET_ROUTINES(Date,string,"TOADDATE") 
READ_AND_SET_ROUTINES(TimeObs,string,"TOADUTIM") 
READ_AND_SET_ROUTINES(PhotomReference,string,"PHOTOREF") 
READ_AND_SET_ROUTINES(Target,string,"TOADOBJE") 
READ_AND_SET_ROUTINES(Airmass,double,"TOADAIRM") 
READ_AND_SET_ROUTINES(SignalToNoise23,double,"USNOSB23") 

/* calibration */

// different stories for different Zero Point ...
// as measured with USNO Catalog
READ_AND_SET_ROUTINES(ZeroPoint,double,"ZEROUSNO")
// as computed from instrument specifications
READ_AND_SET_ROUTINES(Zerop,double,"TOADPZPT") 

// a Zero Point that was once thought to be good ......
READ_AND_SET_ROUTINES(ZP0, double,"ZP0") 
HAS_ROUTINE(ZP0,"ZP0") 

// This one is supposed to be better
READ_AND_SET_ROUTINES(ZP, double,"ZP") 
HAS_ROUTINE(ZP,"ZP") 

// as decided once for 1 image (photometric ref) and then propagated according to the relationship to this image
// this key is purely internal to TOADS
READ_AND_SET_ROUTINES(ZZZeroP, double,"ZPTOADS") 
HAS_ROUTINE(ZZZeroP,"ZPTOADS") 
REMOVE_ROUTINE(ZZZeroP,"ZPTOADS") 


double ReducedImage::AnyZeroPoint() const
{
  
  if ( HasZZZeroP() )
    {
      cerr << "AnyZeroPoint for image : " << FitsName() 
	   << " uses internal TOADS zero point ZZZeroP. " << endl ;
      return ( ZZZeroP() );
    }

  if ( HasZP() )
    {
      cerr << "AnyZeroPoint for image : " << FitsName() 
	   << " uses ZP key. " << endl ;
      return ( ZP() );
    }

  if ( HasZP0() )
    {
      cerr << "AnyZeroPoint for image : " << FitsName() 
	   << " uses ZP0 key. " << endl ;
      return ( ZP0() );
    }

  cerr << "No regular Zero Point (ZP or ZP0) in image " 
       << FitsName() << ", will use what I can (TOADPZPT) " << endl ;
  return (Zerop());
    
}





SET_ROUTINE(OldGain, double, "OLDGAIN");

/* KEYS added after SExtractor Catalog is done */ 
REMOVE_ROUTINE(Seeing,"SESEEING") ;
READ_AND_SET_ROUTINES(Seeing,double,"SESEEING") ;

REMOVE_ROUTINE(SESky, "SEXSKY");
SET_ROUTINE(SESky, double, "SEXSKY");

REMOVE_ROUTINE(SESigma, "SEXSIGMA");
SET_ROUTINE(SESigma, double, "SEXSIGMA");
/*  ******************************** */


/* KEYS added after back subtraction */ 
bool  ReducedImage::BackSub() const
{

  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (head.HasKey("BACK_SUB")) return bool(head.KeyVal("BACK_SUB"));//
    }
  return(false);
}
SET_ROUTINE(BackSub, bool, "BACK_SUB");


REMOVE_ROUTINE(BackLevel,"BACKLEV");
bool ReducedImage::SetBackLevel(const double &Value, const string Comment)
{
  bool ok_back = set_double_key(Value,"BACKLEV","BackLevel",Comment) ;
  if  ( fabs(Value) < 1.e-10 )
    {
      FitsHeader head(FitsName(),RW);
      bool ok_head = head.AddOrModKey("BACK_SUB",true,"Background subtracted");
      ok_back = ok_back && ok_head;
    }
  return(ok_back);
}

double ReducedImage::BackLevel() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      // used to tell if the background map computed by sextractor was subtracted off
      if (head.HasKey("BACK_SUB") && bool(head.KeyVal("BACK_SUB"))) return 0;
      if (head.HasKey("BACKLEV")) return double(head.KeyVal("BACKLEV")); // set by transformedimage and such
      if (head.HasKey("SEXSKY")) return double(head.KeyVal("SEXSKY"));   // set when making the catalog.
      if (head.HasKey("SKYLEV")) return double(head.KeyVal("SKYLEV"));   // set when flatfielding
      cerr << " no way to figure out BackLevel in " << Name() << endl;
      return 0;
    }
  else //
    {
      cerr << " ReducedImage::BackLevel  could not read file " << FitsName() << endl;
      return 0;

    }
}

bool ReducedImage::IsSkySub() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {// DH: only BackSub indicates that the back has been subtracted....... this seems to be a patch accounting for bad key in the swarp reference ........
      if ((head.HasKey("BACK_SUB") && bool(head.KeyVal("BACK_SUB"))) ||
	  (head.HasKey("BACKLEV") && int(head.KeyVal("BACKLEV")) == 0))
	return true;
    }
  return false;	  
}

/*  ******************************** */


double ReducedImage::SigmaBack() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (head.HasKey("SKYSIGMA")) return double(head.KeyVal("SKYSIGMA")); // set when transforming the image
      if (head.HasKey("SEXSIGMA")) return double(head.KeyVal("SEXSIGMA")); // set when making the catalog.
      if (head.HasKey("SKYSIGEX")) return double(head.KeyVal("SKYSIGEX")); // set when flatfielding
      cerr << " no way to figure out SigmaBack in " << Name() << endl;
      return 0;
    }
  else
    {
      cerr << " ReducedImage::SigmaBack  could not read file " << FitsName() << endl;
      return 0;
    }
}

SET_ROUTINE(SigmaBack, double, "SKYSIGMA");

double ReducedImage::NoisePow() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (!head.HasKey("NOISEPOW")) return 0;
      return double(head.KeyVal("NOISEPOW"));
    }
  return 0;
}

SET_ROUTINE(NoisePow, double, "NOISEPOW");


double ReducedImage::OriginalSkyLevel() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (!head.HasKey("SKYLEV")) return BackLevel();
      return read_double_key("SKYLEV","OriginalSkyLevel");
    }
  return 0;
}

SET_ROUTINE(OriginalSkyLevel, double, "SKYLEV");

double ReducedImage::OriginalSaturation() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (!head.HasKey("ORIGSATU")) return Saturation();
      return read_double_key("ORIGSATU","OriginalSaturation");
    }
  return 0;
}

SET_ROUTINE(OriginalSaturation, double, "ORIGSATU");



double ReducedImage::Saturation() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (head.HasKey("SATURLEV")) return double(head.KeyVal("SATURLEV"));//set by flatfielding and catalog
      if (head.HasKey("WELLDEPT")) return double(head.KeyVal("WELLDEPT"));// LBL already flatfielded style
      if (head.HasKey("SATURATE") )return double(head.KeyVal("SATURATE"));//
    }
  double gain = head.KeyVal("TOADGAIN");
  return 60000*gain;
}

SET_ROUTINE(Saturation, double, "SATURLEV");

double ReducedImage::FlatFieldNoise() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (head.HasKey("FLATNOIS")) 
	return double(head.KeyVal("FLATNOIS"));
      /*double sigma = SigmaBack();
      double sky = OriginalSkyLevel();
      double ro = ReadoutNoise();
      return sqrt(fabs((sigma*sigma-ro*ro-sky)))/sky;*/
    }
  return 0;
}

SET_ROUTINE(FlatFieldNoise, double, "FLATNOIS");

double ReducedImage::ProfileError() const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      if (head.HasKey("SIGPSF")) return double(head.KeyVal("SIGPSF"));
      if (head.HasKey("SESIGX") && head.HasKey("SESIGY"))
	{
	  double sigmax = head.KeyVal("SESIGX");
	  double sigmay = head.KeyVal("SESIGY");
	  return 0.02/(sigmax * sigmay); // 2% error on PSF
	}
      if (head.HasKey("SESEEING"))
	{
	  double sigma = head.KeyVal("SESEEING");
	  return 0.02/(sigma * sigma); // 2% error on PSF
	}
      return 0;
    }
  return 0;
}

SET_ROUTINE(ProfileError, double, "SIGPSF");

bool ReducedImage::GetPsfShapeParams(double &SigmaX, double &SigmaY, double &ThetaXY) const
{
  FitsHeader head(FitsName());
  if (head.IsValid())
    {
      bool status = false;
      if (head.HasKey("SESIGX"))  {SigmaX = head.KeyVal("SESIGX");status=true;}  // set by sextractor
      if (head.HasKey("SESEEING"))  {SigmaX = head.KeyVal("SESEEING");status=true;}  // set by sextractor
      if (!status) cerr << " no way to figure out SigmaX in " << Name() << endl;
      status = false;
      if (head.HasKey("SESIGY"))  {SigmaY = head.KeyVal("SESIGY");status=true;}  // set by sextractor
      if (head.HasKey("SESEEING"))  {SigmaY = head.KeyVal("SESEEING");status=true;}  // set by sextractor
      if (!status) cerr << " no way to figure out SigmaY in " << Name() << endl;
      status = false;
      if (head.HasKey("SERHO"))  {ThetaXY = head.KeyVal("SERHO");status=true;}  // set by sextractor
      if (head.HasKey("SESEEING"))  {ThetaXY = 0;status=true;}

      return status;
    }
  else //
    {
      cerr << " ReducedImage::GetPsfShapeParams could not read file " << FitsName() << endl;
      return false;
    }

  cerr << " ReducedImage::GetPsfShapeParams could not read file " << FitsName() << endl;
  return false;

}

bool ReducedImage::SetPsfShapeParams(const double &SigmaX, const double &SigmaY, 
				     const double &ThetaXY, const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (
	      header.IsValid() && 
	      header.AddOrModKey("SESIGX", SigmaX, Comment) &&
	      header.AddOrModKey("SESIGY", SigmaY, Comment) &&
	      header.AddOrModKey("SERHO", ThetaXY, Comment));
    }
  cerr << " ReducedImage::SetShapeParams : could not write at least one of"
       << " shape parameters keys for file " << fileName << endl;
  return false;
}

bool ReducedImage::GetRaDecEpoch(double &Ra, double &Dec, double &Epoch) const
{
  FitsHeader header(FitsName());
  if (header.IsValid())
    {
      bool status = false;
      if (header.HasKey("TOADRASC")) {Ra = RaStringToDeg(header.KeyVal("TOADRASC"));status=true;}
      if (!status) cerr << " no way to figure out Ra in " << Name() << endl;
      status = false;
      if (header.HasKey("TOADDECL")) {Dec = RaStringToDeg(header.KeyVal("TOADDECL"));status=true;}
      if (!status) cerr << " no way to figure out Dec in " << Name() << endl;
      status = false;
      if (header.HasKey("TOADEQUI")) {Epoch = header.KeyVal("TOADEQUI");status=true;}
      if (!status) cerr << " no way to figure out Epoch in " << Name() << endl;
      return status;
    }
  cerr << " ReducedImage::GetRaDecEpoch could not read file " << 
    FitsName() << endl;
  return false;
}



bool ReducedImage::SetRaDecEpoch(const double &Ra, const double &Dec,
		   const double &Epoch, const string Comment)
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader header(fileName, RW);
      return (
	      header.IsValid() && 
	      header.AddOrModKey("TOADRASC", Ra, Comment) &&
	      header.AddOrModKey("TOADDECL", Dec, Comment) &&
	      header.AddOrModKey("TOADEQUI", Epoch, Comment));
    }
  cerr << " ReducedImage::SetRaDecEpoch : could not write at least one of"
       << " shape parameters keys for file " << fileName << endl;
  return false;
}


//  *********specific routines****************
double ReducedImage::RaDeg2000() const
{
  FitsHeader header(FitsName());
  double radeg, decdeg;
  RaDec2000(header, radeg, decdeg);
  return radeg;
}

bool ReducedImage::SetRaDeg2000(const double &Value, const string Comment)
{
  SetEpoch(2000);
  return set_string_key(RaDegToString(Value),"TOADDECL","SetDecDeg2000", Comment);
}

double ReducedImage::DecDeg2000() const
{
  FitsHeader header(FitsName());
  double radeg, decdeg;
  RaDec2000(header, radeg, decdeg);
  return decdeg;
}

bool ReducedImage::SetDecDeg2000(const double &Value, const string Comment)
{
  SetEpoch(2000);
  return set_string_key(DecDegToString(Value),"TOADDECL","SetDecDeg2000", Comment);
}

double ReducedImage::JulianDate() const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader head(FitsName());
      if (head.HasActualKey("JD")) return double(head.KeyVal("JD"));
      return JulianDay(head);
    }
  else
    {
      cerr << " ReducedImage::JulianDate() cannot compute JulianDay for file " 
	   << fileName << endl;
    }
  return UNDEFINED;
}

bool ReducedImage::SetJulianDate(const double &Value, const string Comment)
{
  return set_double_key(Value,"JD","SetJulianDate", Comment);
}

bool ReducedImage::SetUsablePart(const Frame &NewFrame)
{
  FitsHeader head(FitsName()); 
  NewFrame.WriteInHeader(head);
  return true;
}

Frame ReducedImage::UsablePart() const
{
  FitsHeader head(FitsName());
  if (head.IsValid()) return Frame(head);
  else return Frame();
}

Frame ReducedImage::PhysicalSize() const
{
  FitsHeader head(FitsName());
  if (head.IsValid()) return Frame(head, WholeSizeFrame);
  else return Frame();
}

string ReducedImage::Telescope() const
{
  FitsHeader head(FitsName());
  if (head.IsValid()) return TelescopeName(head);
  else return "telescope_prout";
}

string ReducedImage::Instrument() const
{
  FitsHeader head(FitsName());
  if (head.IsValid()) return InstrumentName(head);
  else return "ccd_prout";
}

bool ReducedImage::SameChipFilterInst(const ReducedImage &Other,const bool Warn) const
{
  return FitsHeader(FitsName()).SameChipFilterInst(FitsHeader(Other.FitsName()), Warn);
}

bool ReducedImage::SameChipFilterInstNight(const ReducedImage &Other,const bool Warn) const
{
  return ((SameChipFilterInst(Other, Warn)) && 
	  (fabs(JulianDate()-Other.JulianDate()) ) < 0.5);
}

bool ReducedImage::SamePhysicalSize(const ReducedImage &OtherImage) const
{
  return (PhysicalSize() == OtherImage.PhysicalSize());
}



#include "wcsutils.h"
double ReducedImage::OverlapArcmin2(const ReducedImage& Other) const
{
  return Arcmin2Overlap(FitsHeader(FitsName()), FitsHeader(Other.FitsName()));
}

Gtransfo *ReducedImage::RaDecToPixels() const
{
  Gtransfo *pix2radec(0);
  FitsHeader head(FitsName());
  if (!WCSFromHeader(head, pix2radec))
    {
      cerr << " ReducedImage::RaDecToPixels() did not find the expected WCS in FITS header " 
	   << FitsName() << endl;
      return 0;
    }  
  Gtransfo *inv =  pix2radec->InverseTransfo(0.01, Frame(head));
  delete pix2radec;
  return inv;
}

Gtransfo *ReducedImage::PixelsToRaDec() const
{
  Gtransfo *pix2radec(0);
  FitsHeader head(FitsName());
  if (!WCSFromHeader(head, pix2radec))
    {
      cerr << " ReducedImage::PixelsToRaDec() did not find the expected WCS in FITS header " 
	   << FitsName() << endl;
      return 0;
    }
  return pix2radec;
}


void ReducedImage::dump(ostream & s) const 
{
  s  << "*******ReducedImage ******" << endl;
  s << " Image:  " << Name() << endl;
  s << " Band " << Band() << endl;
  s << " Target " << Target() << endl;  
  s << " Dead pixels: " << FitsDeadName() << endl ;  
  s << " Saturated stars pixels: " << FitsSaturName() << endl ;  
  s << " Stars list: " << CatalogName() << endl ;
  s << " Seeing = " << Seeing() << endl;
  s << " Background level = "<< BackLevel() << endl;
  s << " Background sigma = " << SigmaBack() << endl;
  s << " Saturation = " << Saturation() << endl;
}

bool ReducedImage::IsGoodImage() const
{
  if (!actuallyReduced) 
    {
      cerr << " " << Name() << " is not reduced " << endl;
      return false;
    }

  if (BackLevel() >= Saturation())
    {
      cerr << " BackLevel > Saturation for " << Name() << endl;
      return false;
    }

  if (OriginalSkyLevel() >= OriginalSaturation())
    {
      cerr << " OriginalSkyLevel > OriginalSaturation for " << Name() << endl;
      return false;
    }

  if (OriginalSkyLevel() <= 0)
    {
      cerr << " OriginalSkyLevel is undefined for " << Name() << endl;
      return false;
    }

  if (Seeing() <= 0)
    {
      cerr << " Seeing is undefined for " << Name() << endl;
      return false;
    }
  if (Saturation() <= 0)
    {
      cerr << " Saturation is undefined for " << Name() << endl;
      return false;
    }

  if (PixelSize() <= 0)
    {
      cerr << " PixelSize is undefined for " << Name() << endl;
      return false;
    }

  if (Exposure() <= 0)
    {
      cerr << " Exposure is undefined for " << Name() << endl;
      return false;
    }

  return true;
}


ReducedImage::~ReducedImage()
{
  /* the call also exists in the base class destructor (~DbImage), but
  when we reach it, the object in hand is no longer a derived object, and hence
  the  called Streamer (in the root framework) is DbImage::Streamer.
  This is not a bad C++ feature, since when we reach DbImage::~DbImage, 
  the private data of ReducedImage (if any) is no longer meaningful :
  if this destructor contained delete satements, they would be called 
  before ~DbImage. hope this is clear, but just remember that if you want your
  own inheritant of ReducedImage save its private data, you *DO* have to call
  writeEverythingElse in its destructor. This routine is protected against
  multiple calls for the same object. */
  writeEverythingElse();
}


//! sniffs if gain multiplied (by simply checking if gain is 1) and 
//! if not, multiplies it 
double ReducedImage::MultiplyGain()
{
  double gain = Gain();
  if (gain==1)
    {
      cout << "ReducedImage::GainMultiply Gain == 1 !!!! Nothing done"<< endl;
      return gain;
    }
  else
    {
      cout << "ReducedImage::GainMultiply => multiply the image"<< Name() <<" by the gain " << gain << endl;
      FitsImage img(FitsName(), RW);
      img *= gain;
    }
  SetOldGain(gain);
  SetGain(1.0);
  return gain;
}


  

// this routine  with a (very) awful signature is done this way
// because it is called for both dead and satur maps. 
// It may be used by other classes


bool BoolImageOr(ReducedImageList &List,
                 string (ReducedImage::*InFitsFileName)() const,
		 bool (ReducedImage::*MakeInFits)() ,
		 const string OutFitsName)
{
  if (List.size() == 0) return false;
  vector<string> fitsNames(List.size());
  Image *sum = NULL;
  string firstName;
  int count = 0;
  for (ReducedImageIterator i=List.begin(); i != List.end(); ++i)
    {
      ReducedImage *reducedImage = *i;
      string currentName = (reducedImage->*InFitsFileName)();
      if (!FileExists(currentName)) (reducedImage->*MakeInFits)();
      FitsImage currentFits(currentName);
      if (!currentFits.IsValid())
	{
	  cerr << " cannot find " << currentName << " for ReducedImage " 
	       << reducedImage->Name() << endl;
	  continue;
	}
      if (sum == NULL)
	{
	  sum = new Image(currentFits);
	  firstName = currentName;
	  count = 1;
	}
      else 
	{
	  *sum += currentFits;
	  ++count;
	}
    }
  if (!sum) 
    {
      cerr << " could not make " << OutFitsName << endl;
      return false;
    }
  if (count != int(List.size()))
    {
      cerr << " could not find all ingredients for " << OutFitsName  << endl;
    }
  sum->Simplify(0.1);
  FitsImage out(OutFitsName, FitsHeader(firstName), *sum);
  delete sum;
  out.ModKey("BITPIX",8);
  out.AddOrModKey("IMCOUNT",count);
  return true;
}

bool BoolImageAnd(ReducedImageList &List,
                 string (ReducedImage::*InFitsFileName)() const,
		 bool (ReducedImage::*MakeInFits)() ,
		 const string OutFitsName)
{
  if (List.size() == 0) return false;
  vector<string> fitsNames(List.size());
  Image *sum = NULL;
  string firstName;
  int count = 0;
  for (ReducedImageIterator i=List.begin(); i != List.end(); ++i)
    {
      ReducedImage *reducedImage = *i;
      string currentName = (reducedImage->*InFitsFileName)();
      if (!FileExists(currentName)) (reducedImage->*MakeInFits)();
      FitsImage currentFits(currentName);
      if (!currentFits.IsValid())
	{
	  cerr << " cannot find " << currentName << " for ReducedImage " 
	       << reducedImage->Name() << endl;
	  continue;
	}
      if (sum == NULL)
	{
	  sum = new Image(currentFits);
	  firstName = currentName;
	  count = 1;
	}
      else 
	{
	  *sum += currentFits;
	  ++count;
	}
    }
  if (!sum) 
    {
      cerr << " could not make " << OutFitsName << endl;
      return false;
    }
  if (count != int(List.size()))
    {
      cerr << " could not find all ingredients for " << OutFitsName  << endl;
    }
  double cut = List.size()-0.5;
  sum->Simplify(cut);
  FitsImage out(OutFitsName, FitsHeader(firstName), *sum);
  delete sum;
  out.ModKey("BITPIX",8);
  out.AddOrModKey("IMCOUNT",count);
  return true;
}


bool IncreasingSeeing(const ReducedImage* one, const ReducedImage* two)
{ return (one->Seeing() < two->Seeing());}

bool IncreasingJulianDate(const ReducedImage* one, const ReducedImage* two)
{ return (one->JulianDate() < two->JulianDate());}

bool IncreasingResolution(const ReducedImage* one, const ReducedImage* two)
{
  return (one->PixelSize() < two->PixelSize());
}

bool DecreasingArea(const ReducedImage *one, const ReducedImage *two)
{ 
  FitsHeader head1(one->FitsName());
  FitsHeader head2(two->FitsName());
  return (Arcmin2Area(head1) > Arcmin2Area(head2));
}


/************************************ ReducedImageList *******************/

//#include "imagelist.cc"
//template class ImageList<ReducedImage> ; // to force instanciation 


Frame CommonFrame(ReducedImageList &RedList)
{
  sort(RedList.begin(), RedList.end(), DecreasingArea);
  Frame frame_common(RedList.front()->UsablePart());
  for (ReducedImageCIterator im = RedList.begin(); im != RedList.end(); ++im)
    frame_common *= Frame((*im)->UsablePart());
  cout << " Common frame for ReduceImageList: " << frame_common;
  return frame_common;
}


#ifdef USE_ROOT
/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class ReducedImage;
LINKDEF_CONTENT : #pragma link C++ class vector<ReducedImage*>;
LINKDEF_CONTENT : #pragma link C++ class ImageList<ReducedImage>-;
LINKDEF_CONTENT : #pragma link C++ function ReducedImageRead(const char *);
LINKDEF_CONTENT : #pragma link C++ class string;
LINKDEF_CONTENT : ostream& operator << (ostream& , const string &);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream& , const string&);


there is a problem with vector<RedducedImage*>::iterator, as well as 
the trick using a typedef.


*/

ClassImp(ReducedImage);



#include "root_dict/reducedimagedict.cc"

#include <TFile.h>
#include <TKey.h>
#endif /* USE_ROOT */



ReducedImage* ReducedImageRead(const char *Name)
{
  return ReducedImageRead(string(Name));
}
			  

ReducedImage* ReducedImageRead(const string &Name)
{
#ifdef USE_ROOT
  DbImage ri(Name);
  if (!ri.IsValid()) 
    {cerr << Name << " : No such ReducedImage " << endl; return NULL;}
  string fileName = ri.EverythingElseFileName();
  if (!FileExists(fileName)) 
    { cerr << Name << " : does not have a " << fileName << " file " << endl; 
    return NULL;
    }
  TFile tfile(fileName.c_str());
  TIter nextkey(tfile.GetListOfKeys());
  TKey *key = (TKey*)nextkey();
  ReducedImage *p;
  if (key) 
    {
      p = dynamic_cast<ReducedImage*>(key->ReadObj());
      if (!p) 
	{
	  cerr << " cannot find a ReducedImage in " << fileName << endl;
	  return p;
	}
      if (p->Name() != Name) 
	{
	  cerr << " something weird happened when reloading DbImage " 
	       << Name << endl
	       << " name in the file : " << p->Name() << endl
	       << " name requested " << Name << endl;
	}
      p->init_from_name();
    }
  else 
    {
      cerr << " could not read : " << fileName << endl;
    }
  tfile.Close(); 
  return p;
#else
  cerr << "for image :" << Name << ", ReducedImageRead only works with USE_ROOT  defined " << endl;
  return NULL;
#endif
}

  



