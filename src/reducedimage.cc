#include <iostream>
#include <memory> // for auto_ptr

#include "fitsimage.h"
#include "reducedimage.h"
#include "astroutils.h"
#include "fitstoad.h"
#include "sestar.h"
#include "cluster.h"
#include "gtransfo.h"
#include "wcsutils.h"
#include "imageutils.h" // ConvolveSegMask
#include "fitsexception.h"


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


ReducedImage *ReducedImage::Clone() const
{
  string storedTypeName = TypeName();
  if (TypeName() != storedTypeName)
    {
      std::cerr << " you are Cloning() " << Name() << " of type " 
		<< storedTypeName << " using ReducedImage::Clone() " 
		<< std::endl
		<< " May be class " << storedTypeName 
		<< " misses a Clone() method?" << std::endl;
    }
  return new ReducedImage(*this);
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


#define SATUR_COEFF 0.98


/* Pour remplir le carnet d'ordres de SExtractor */
/* si le fond est deja soustrait, on le sauve pas, et on prend = cste
   (SExtractor le recalcule quand meme) sinon, SExtractor le
   calculera, l'utilisera, et on le sauvera.
*/

void ReducedImage::FillSExtractorData(ForSExtractor & data, 
				      bool  fond_deja_soustrait, 
				      bool sauver_fond,
				      bool use_sigma_header)
{
  // what concerns on-the-flight decompression
  const char *tmpdir = getenv("IMAGE_TMP_DIR");
  if (tmpdir) data.TempDir = string(tmpdir); else data.TempDir = Dir();
  data.UniqueName = Name();
    
  // real stuff now.
  data.saturation = (int) (SATUR_COEFF*Saturation());
  data.sigma_back  = -1 ;
  if (use_sigma_header)
    {
      if ( SigmaBack() > 1.e-10 )
	{
	  data.sigma_back  = SigmaBack() ;
	  cerr << "Using the background sigma written in " 
	       << Name() << " header to compute detctions levels: " 
	       << data.sigma_back << endl ;
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

      /* The "BadImage" flags 0 weight pixels. SExtractor ignore 0 weight
	 pixels for computing object scores, but does NOT flag objects 
	 affected by 0 weight pixels. But it flags objects which make use
	 of pixels flagged in the flag images. This is why we build this 
	 "BadImage", just for SExtractor. It is deleted when we are done.
      */
      if (MakeBad())
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


// Pour le double mode simplifie : pas de carte de bad.
// a appeler sur l'image de photometrie (ie la 2eme image)
// fond en mode RELATIVE (dans prefs.h) == AUTO (dans datacard.sex) == carte de fond re-calculee// a appeler sur image de photometrie
// j'ai enleve les bad non necessaires pour la comparaison a Terapix
// fond en mode RELATIVE (dans prefs.h) == AUTO (dans datacard.sex) == carte de fond re-calculee
void ReducedImage::FillSExtractorData_2(ReducedImage & rim_det, 
					ForSExtractor & data)
{
  // what concerns on-the-flight decompression
  const char *tmpdir = getenv("IMAGE_TMP_DIR");
  if (tmpdir) 
    {// risque mais si c'est ce qu'on veut
      data.TempDir_0 = string(tmpdir); 
      data.TempDir_1 = string(tmpdir); 
    }
  else {
    data.TempDir_0 = rim_det.Dir();
    data.TempDir_1 = Dir();
  }
  data.UniqueName_0 = rim_det.Name();
  data.UniqueName_1 = Name();
    
  // saturation from measurement image
  data.saturation = (int) (SATUR_COEFF*Saturation());
  // pour comparaison avec datacard terapix
  //data.saturation = 40000. ;
  data.sigma_back  = -1 ;
  data.FitsFileName_0 = rim_det.FitsImageName(Calibrated).c_str();
  data.FitsFileName_1 = FitsImageName(Calibrated).c_str();
  
  data.backmean=0.0;
  // a l'air d'etre comme ca dans datacard terapix: data.back_type_manual = false ; 
  data.back_type_manual = false ; 

  //on ne sauve pas le fond car on ne le soustrait pas
  //data.FitsBackName = FitsBackName();
  data.FitsMiniBackName = FitsMiniBackName();

  if (HasWeight())
    {
      cout << "Weighting for photometry from " << FitsWeightName() << endl ;
      data.FitsWeightName_1 = FitsWeightName();
      // pas de flag image (mask image ) ici
    }
  else 
    {     
      data.FitsWeightName_1 = "";
    }

  if (rim_det.HasWeight())
    {
      cout << "Weighting for detection from " << rim_det.FitsWeightName() << endl ;
      data.FitsWeightName_0 = rim_det.FitsWeightName(); 
      // pas de flag image ici  
    }
  else 
    {     
      data.FitsWeightName_0 = "";
    }


  cerr << "Images Detection : " << data.FitsFileName_0 << " Mesure : " << data.FitsFileName_1 << endl ;
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

#define READ_IF_EXISTS(VAR,TAG,TYPE) \
if (cards.HasKey(TAG)) VAR=cards.TYPE(TAG)

#include "imageback.h"

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

   if (!HasImage())
     {
       cerr << " cannot make catalog without image " << endl;
       return false;
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
       if (!MakeWeight()) return false;
     }




  DataCards cards(DefaultDatacards());
  int use_poloka_back = 1;
  READ_IF_EXISTS(use_poloka_back, "USE_POLOKA_BACK", IParam);

  
   // the background is subtracted if we do want to subtract it AND  if it wasn't subtracted before. 
  bool soustraire_fond = !(fond_deja_soustrait) && !pas_sub_fond;

  /*
  The background map is saved only if it is subtracted, and we use SExtractor back (rather than Poloka back)
  */

  {

    // DEBUG
    // fix wrong saturation levels (occurs for short exposure times)
    FitsHeader head(FitsName());

    if (IsOfKind<Megacam>(head) && Saturation()< 50000)
      {
	FitsHeader flat(FitsFlatName());
	double fla = flat.KeyVal("FLATSCLA"); // A amplifier
	double flb = flat.KeyVal("FLATSCLB"); // B amplifier
	// actually the inverse of average_flat.
	double flat_mean = min(fla,flb);
	cout << " WARNING : changing saturation  from " << Saturation() << " to ";
	SetSaturation(65536.*flat_mean);
	cout << Saturation() << endl;
      }
  }

  bool sauver_fond = (use_poloka_back == 0);

  ForSExtractor data ;
  FillSExtractorData(data,fond_deja_soustrait,sauver_fond, use_sigma_header);
  // assign the name of the "SEGMENTATION" image. Temporary name.
  data.FitsSegmentationName = Dir()+"seg.fits";
  data.Print();
  double Fond = 0., SigmaFond = 0. ;
  FitsImage *pmasksat = NULL ;
  string nommasksat;
  if (savemasksat)
    {
      FitsHeader imgheader(FitsName());
      pmasksat = new FitsImage(FitsSaturName(), imgheader);
      pmasksat->AddOrModKey("BITPIX",8);
      pmasksat->EnableWrite(false); // in case something goes wrong
    }
  SEStarList stlse ;
  int status = SEStarListMake(data, stlse, Fond, SigmaFond, pmasksat);
  if (status == 0) {return 0;}
  // DEBUG
  if (!FileExists(data.FitsSegmentationName))
    {
      cout << " ReducedImage::MakeCatalog : cannot find object mask " <<data.FitsSegmentationName << endl;
    }
  //END DEBUG


  // we re-flag for saturation here
  int nsat = FlagSaturatedStars(stlse, data.saturation);
  cout << nsat << " (extra) stars flagged as saturated. " << endl ;

  if (pmasksat) 
    { 
    pmasksat->EnableWrite(true);
    cout << "Writing Saturated stars pixels map in " 
	 << FitsSaturName() << endl ;
    cout << " we have " << pmasksat->SumPixels() / pmasksat->NPix()*100
	 << "% pixels flagged as saturated" << endl;
    delete pmasksat;
    pmasksat = NULL;
    }

  SetSESky(Fond, "SExtractor computed background");
  SetSESigma(SigmaFond, "SExtractor computed sigma on background"); 
  // set back to zero in starlist and correct saturation in image

  //compression of the object mask
  FitsImage *segmentationMask = NULL;
  if (!FileExists(data.FitsSegmentationName))
    {
      cout << " ReducedImage::MakeCatalog : cannot find object mask " <<data.FitsSegmentationName << endl;
    }
  else
    {
      // "process" the segmentation mask :  it gets compressed on disk.
      FitsImage sexOutput(data.FitsSegmentationName);
      segmentationMask = new FitsImage(FitsSegmentationName(), sexOutput, sexOutput );
      /* the "segmentation image" has a span depending on the number of objects in the catalog.
	 So we adapt the disk representation to the actual number of entries in the catalog */
      int bitpix = (stlse.size() > 32000)? 32 : 16;
      segmentationMask->AddOrModKey("BITPIX",bitpix); 

      segmentationMask->AddOrModKey("INT_DATA",true, " this is integer data ");
      segmentationMask->Write(); // write it before we mess it up
      // and remove the mask provided  by Sextractor.
      remove(data.FitsSegmentationName.c_str());
    }


  if (soustraire_fond)
    {
      if (use_poloka_back)
	{
	  FitsImage *weight = NULL;
	  int poloka_back_object_mask_border = 5;
	  if (segmentationMask)
	    {
	      READ_IF_EXISTS(poloka_back_object_mask_border,"POLOKA_BACK_OBJECT_MASK_BORDER",IParam);
	      
	      cout << " convolving object mask with bordersize = " << poloka_back_object_mask_border << endl;
	      ConvolveSegMask(*segmentationMask, *segmentationMask, poloka_back_object_mask_border);
	      {
		string cvName = CutExtension(segmentationMask->FileName())
		  +".cv."+FileExtension(segmentationMask->FileName());
		FitsImage toto(cvName,*segmentationMask,*segmentationMask);
	      }
	      // build a weight image with objects masked

	      // "invert" the mask
	      Pixel *end = segmentationMask->end();
	      for (Pixel *pm = segmentationMask->begin(); pm < end; ++pm) if (*pm == 0) *pm=1; else *pm = 0;
	      weight = new FitsImage(FitsWeightName());	  
	      (*weight) *= (*segmentationMask);
	      delete segmentationMask; segmentationMask = NULL;		
	    }
	  else
	    {
	      cout << " could not mask object pixels when computing sky of " << Name() << endl;
	      return false;
	    }
	  if (true)
	    {
	      // read mesh size from datacards
	      int poloka_back_mesh_sizex = 256;
	      int poloka_back_mesh_sizey = poloka_back_mesh_sizex;
	      if (cards.HasKey("POLOKA_BACK_MESH_SIZE"))
		{
		  poloka_back_mesh_sizex = poloka_back_mesh_sizey = cards.IParam("POLOKA_BACK_MESH_SIZE");
		}
	      else
		{
		  READ_IF_EXISTS(poloka_back_mesh_sizex,"POLOKA_BACK_MESH_SIZEX", IParam);
		  poloka_back_mesh_sizey = poloka_back_mesh_sizex;
		  READ_IF_EXISTS(poloka_back_mesh_sizey,"POLOKA_BACK_MESH_SIZEY", IParam);
		}
	      // simple  isn't it?
	      cout << "TOADS: Computing Image Background with mesh = " 
		   <<  poloka_back_mesh_sizex << ',' << poloka_back_mesh_sizey << endl;
	      FitsImage im(FitsName(), RW);
	      
	      ImageBack back(im, poloka_back_mesh_sizex, poloka_back_mesh_sizey, weight);
	      delete weight;
	      { // save the mini back
		FitsImage miniBack(FitsMiniBackName(), back.BackValue());
		miniBack.AddOrModKey("SEXBKGSX", poloka_back_mesh_sizex);
		miniBack.AddOrModKey("SEXBKGSY", poloka_back_mesh_sizey);
		miniBack.AddOrModKey("BITPIX",-32);
		miniBack.AddOrModKey("EXTRAPIX",poloka_back_object_mask_border,
				     " by how many pixels we enlarged the segmentation sex patches");
		miniBack.AddCommentLine("Poloka Computed Miniback (class ImageBack)");
	      }
	      cout << "TOADS: Subtracting Image Background" << endl ;
	      Image *largeBack = back.BackgroundImage();
	      im -= *largeBack;
	      delete largeBack;
	      // DEBUG/TEST
              // im.SetWriteAsFloat();
	    }
	  
	  // set back ("Fond() ") of stars to zero.
	  SetStarsBackground(stlse,0.);
	  // update saturation level in image header
	  SetOriginalSaturation(Saturation(),"Original saturation level before sky subtraction"); 
	  SetSaturation(Saturation()-Fond,"Saturation level corrected from sky subtraction"); 
	  SetBackLevel(0.," Poloka background (with masked objects) subtracted"); // activates BackSub 
	}
      else // subtract sextractor background
	{
	  // set back ("Fond() ") of stars to zero.
	  SetStarsBackground(stlse,0.);
	  // update saturation level in image header
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
	}// end of sex back subtraction
    } // end of back subtraction
  if (segmentationMask) delete segmentationMask; segmentationMask = NULL;

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

  MakeWeight(); // updates a few things

  stlse.write(CatalogName());

  
  // dump some statistics
  cout << "@NUMBER_OF_OBJECTS " << stlse.size() << endl;
  int n_saturated_objects = 0;
  int n_ok = 0;
  int n_flag0 = 0;
  int n_flagbad0 = 0;
  
  for (SEStarIterator it= stlse.begin(); it!=stlse.end(); it++) {
    if( (*it)->IsSaturated() ) n_saturated_objects++;
    if( (*it)->IsOK() ) n_ok++;
    if( (*it)->Flag()==0 ) n_flag0++;
    if( (*it)->FlagBad()==0 ) n_flagbad0++;
    
    
  }
  cout << "@NUMBER_OF_SATURATED_OBJECTS " << n_saturated_objects << endl;
  cout << "@NUMBER_OF_OK_OBJECTS " << n_ok << endl;
  cout << "@NUMBER_OF_FLAG0_OBJECTS " << n_flag0 << endl;
  cout << "@NUMBER_OF_FLAGBAD0_OBJECTS " << n_flagbad0 << endl;
  
  if (stlse.size()==0) // if it's empty, it's a failure
    return false;

  return(status);
}


// version simplifiee pour 2 images : detection + mesure
// pas de recuparation du masque de saturation

bool
ReducedImage::MakeCatalog_2images(ReducedImage & rim_det, bool overwrite, 
				  bool weight_from_measurement_image, 
				  string catalog_name, bool do_segmentation)
{
  int status = 0 ;
   if ( !overwrite )
    { 
      if ( FileExists(ImageCatalogName(SExtractor)))
	{
	  cerr << "catalog already done for image, exiting " << FitsName()
	       << endl ;
	  return true;
	}
    }
  // fabrication de la carte de poids
   if (! HasWeight() )
     {
       MakeWeight();
     }

   if (! rim_det.HasWeight() )
     {
       rim_det.MakeWeight();
     }

  ForSExtractor data ;
  

  FillSExtractorData_2(rim_det,data);
  if (do_segmentation) 
    data.FitsSegmentationName = rim_det.Dir()+"/seg.fits";
  data.Print();
  double Fond_0 = 0., SigmaFond_0 = 0., Fond_1 = 0., SigmaFond_1 = 0. ;
  SEStarList stlse ;
  status = SEStarListMake_2(data, stlse, Fond_0, SigmaFond_0, Fond_1, 
			    SigmaFond_1,  weight_from_measurement_image);

   cerr << "Nombre d'objet detectes : " << stlse.size() << endl ;
   //SetSESky(Fond_1, "SExtractor computed background");
   //SetSESigma(SigmaFond_1, "SExtractor computed sigma on background"); 
  stlse.write(catalog_name);


  return(status);
}



bool ReducedImage::MakeCatalog() 
{
  if (HasCatalog()) return true;
  // if there is a satur map, it may be the result of e.g
  // stacking. Then there is no point in remaking one.
  bool savemasksat = !HasSatur();
  bool ok = MakeCatalog(/*redo_from_beg=*/ false, 
			/* overwrite = */ false, 
			savemasksat,
			/*pas_sub_fond= */false,
			/*use_sigma_header= */  false);

  // why should we call MakeCosmic() here?
  // MakeCosmic();
  return(ok);
}
#include "apersestar.h"
#include "datacards.h"
#include "toadscards.h"

/* routine steps

   1) copy SEstarList to AperSEStarList
   
   2) compute positions and shapes using a gaussian filtering

   3) compute seeing and read apertures

   4) compute aperture fluxes

   we do not do PositionsAndShapes in the same loop as 
   aperture fluxes, because it may be of interest to set
   zero weights for saturated pixels in the flux computation
   (so that the number of bad pixels accounts for saturation),
   but this cannot be done for shape computation


   This routine could now be split in 2 routines to allow the aperture
   fluxes being computed at positions obtained by other means than the
   ones implemented here (e.g.  the catalog from an other image...)

*/


static double sq(const double &x) {return x*x;}



#include "histopeakfinder.h"
#include "fitsslice.h" // for FitsParallelSlices

bool ReducedImage::MakeAperCat()
{
  string aperCatName = AperCatalogName();
  if (FileExists(AperCatalogName())) return true;

  if (!(MakeFits() && MakeWeight() && MakeCatalog()))
    {
      cerr << " cannot make pre-requisite for aper catalogue for " << Name() << endl;
      return false;
    }

  cout << " Entering MakeAperCat() for " << Name () << endl;

  SEStarList seList(CatalogName());
  if (seList.size() == 0)
    {
      cout << " empty SExtractor catalog for " << Name() 
	   << ", no aper catalog " << endl;
      return false;
    }

  //try to figure out a gain (needed for photometry and to compute pos error)
  double gain = 1;
  double toadpixs = 1;
  { // open a block to close the file after use

  FitsHeader  I(FitsName());
  gain = I.KeyVal("TOADGAIN");
  toadpixs = I.KeyVal("TOADPIXS");
  /* cannot go through the the average/variance trick to figure out
     the gain if the image comes from Swarp, which subtracts the
     sky of input images (probably with good reasons)
  */     
  string inst = I.KeyVal("TOADINST");
  bool useGain = I.HasKey("USEGAIN") && bool(I.KeyVal("USEGAIN"));
  if (useGain)
    cout << " INFO : found the USEGAIN key and using TOADGAIN " << endl;
  else if (inst != "Swarp" && I.HasKey("SEXSKY") && I.HasKey("SEXSIGMA"))
    {
      double sexsky = I.KeyVal("SEXSKY");
      double sexsigma = I.KeyVal("SEXSIGMA");
      if (sexsky > 0) gain = sexsky/(sexsigma*sexsigma);
    }
  if (gain == 0)
    {
      cout << " cannot figure out a gain" << endl;
      gain = 1;
    }
  cout << " Assuming a gain of " << gain << endl;
  }

  //this code now "slices" the input images in order to accomodate monsters

  int ySliceSize = 500;
  int sliceOverlap = 110;

  AperSEStarList apercat(seList); // copy of SExtractor catalog.

  // first pass on the pixels : compute shapes (and positions)
  {

    FitsParallelSlices slices(ySliceSize,sliceOverlap);

    slices.AddFile(FitsName());

    if (!bool(slices[0]->KeyVal("BACK_SUB")))
      {
	cout << " ReducedImage::MakeAperCat : this code cannot accomodate images with background left" << endl;
	return false;
      }
    
    slices.AddFile(FitsWeightName());
    
    
    double yStarMin = 0;
    do  // loop on image slices
      {
	double yStarMax = (slices.LastSlice())? slices.ImageJ(ySliceSize): slices.ImageJ(ySliceSize)-55;
	
	double offset = slices.ImageJ(0);
	for (AperSEStarIterator i = apercat.begin(); i != apercat.end(); ++i)
	  {
	    AperSEStar &s = **i;
	    if (s.y < yStarMax && s.y >= yStarMin)
	      {
		s.y -= offset;
		s.ComputeShapeAndPos(*slices[0], *slices[1], gain);
		s.y += offset;
	      }
	  }
	yStarMin = yStarMax;
      }
    while(slices.LoadNextSlice());
  } // end of first pass


  // compute the "seeing" from the star cluster locus in object shape space
  double seeing = Seeing();
  double xSize=0, ySize=0, corr = 0;
  StarScoresList scores;
  if (FindStarShapes(apercat, 20, xSize, ySize, corr, scores))
    {
      double mxx = sq(xSize);
      double myy = sq(ySize);
      double mxy = corr*xSize*ySize;
      seeing = pow(mxx*myy-sq(mxy),0.25);
    }
  else
    {
      cout << " MakeAperCatcould not figure out a seeing " << endl;
      apercat.write(Dir()+"/aperse_head.list");
      throw(PolokaException("MakeAperCat could not figure out a seeing "));
      return false;
    }
  cout << Name() << " star shapes " << xSize << ' ' << ySize << ' ' 
       << corr  << endl;
  cout << Name () << ": old seeing " <<  Seeing() << " new " << seeing << endl;
  
  SetGFSeeing(seeing," Seeing from a Gaussian fit to objects");
    
  // get aperture radius, either from datacards or provide defaults
  bool fixed_aper_rads = false;
  vector<double> rads;
  DataCards cards(DefaultDatacards());
  
  // FIXED APERTURE RADIUS, in ARCSECONDS ... 
  // may be used for the calibration.
  if (cards.HasKey("FIXED_APER_RADS")) 
    {
      cout << " using FIXED APER RAD values: ";
      fixed_aper_rads = true;
      int n = cards.NbParam("FIXED_APER_RADS");
      for(int i=0;i<n;i++) 
	{
	  double r = cards.DParam("FIXED_APER_RADS",i);
	  r /= toadpixs;
	  cout << " " << r;
	  rads.push_back(r);
	}
      cout << "\n";
    }
  else if (cards.HasKey("APER_RADS")) 
    {
      int n = cards.NbParam("APER_RADS");
      for (int i=0; i < n ; ++i) rads.push_back(cards.DParam("APER_RADS",i));
    }
  else
    {
      cout << " no APER_RADS card in " << DefaultDatacards()  << endl
	   << " resorting to internal defaults " << endl;
      double def[] = {2.,2.5,3.,3.5,4.,5.,7.5,10.,15.,20.};
      // canadian rads
      //      double def[] = {2.705, 4.057, 5.409, 6.762, 8.114, 10.819, 13.524, 16.228, 18.933, 21.638, 27.047};
      for (unsigned i =0; i < sizeof(def)/sizeof(def[0]); ++i)
	rads.push_back(def[i]);
    }

  unsigned nrads = rads.size();
  if(!fixed_aper_rads) 
      for (unsigned k=0; k < nrads; ++k) rads[k] *= seeing;
  
  double maxNeighborDist = rads[nrads-1];

  apercat.SetNeighborScores(*AperSE2Base(&apercat),maxNeighborDist);

  // we are all set for second traversal of pixels : aperure photometry
  {
    int safetyMargin = int(rads[nrads-1])+2;

    FitsParallelSlices slices(ySliceSize,2*safetyMargin);
    slices.reserve(6); // so that slices[x] does not change after a push_back.... used below


    slices.AddFile(FitsName());
    Image &I = (*slices[0]);
    slices.AddFile(FitsWeightName());
    Image &W = (*slices[1]);

    // set weight = 0 for saturated pixels ?
    // yes. If not, saturated spikes within apertures may contribute
    // to the flux without any notice.
    bool hasSatur = HasSatur();
    Image *satur = NULL;
    if (hasSatur)
      {
	slices.AddFile(FitsSaturName());
	satur = slices.back();
      }

    // nasty temporary hack allowing not to use the slow and buggy cosmic filter:
    // either it is already done and we use the output, or it is not done and we just go
    // without cosmics.
    // To be "exception safe", we declare an auto_ptr to clear the allocated memory if any:
    auto_ptr<Image> del_cosmic(NULL);
    Image *cosmic = NULL;
    if (FileExists(FitsCosmicName()))
      {
	cerr << "Using cosmic image " << endl ;
	slices.AddFile(FitsCosmicName());
	cosmic = slices.back();
      }
    else
      {
	cosmic = new Image(slices[0]->Nx(), slices[0]->Ny()); 
	del_cosmic.reset(cosmic);  // will delete cosmic when del_cosmic goes out of scope
      }
    

    // we may also miss the "segmentation image"
    // same trick with exception safe call to delete using auto_ptr
    auto_ptr<Image> del_segmentation(NULL);
    Image *segmentation = NULL; 
    if (FileExists(FitsSegmentationName()))
      {
	slices.AddFile(FitsSegmentationName());
	segmentation = slices.back();
      }
    else
      {
	segmentation = new Image(slices[0]->Nx(), slices[0]->Ny()); 
	del_segmentation.reset(segmentation);
// will delete segmentation when del_segmentation goes out of scope
      }

    // for the y limits take something 
    // larger than the image so that all objects will be processed
    double yStarMin = -100;
    do // loop on image slices
      {
	double yStarMax = (slices.LastSlice())? slices.ImageJ(ySliceSize): slices.ImageJ(ySliceSize)-safetyMargin;

	if (slices.LastSlice()) yStarMax += 100;

	// saturated pixels are considered as bad pixels.
	if (hasSatur) W *= (1.-(*satur));

	double offset = slices.ImageJ(0);
	for (AperSEStarIterator i = apercat.begin(); i != apercat.end(); ++i)
	  {
	    AperSEStar &s = **i;
	    if (s.y < yStarMax && s.y >= yStarMin)
	      {
		s.y -= offset;
	       for (unsigned k = 0; k < nrads; ++k) s.ComputeFlux(I,W,*cosmic,*segmentation,gain,rads[k]); 
		s.y += offset;
	      }
	  }
	yStarMin = yStarMax;
      }
    while(slices.LoadNextSlice());

  } // end of second  pass through the pixels. End of block releases image buffers.

    // check that each star was processed only once
    for (AperSEStarIterator i = apercat.begin(); i != apercat.end();)
      {
	const AperSEStar &s = **i;
	if (s.apers.size() == nrads)
	  {
	    ++i; continue;
	  }
	cout << " ERROR : ReducedImage::MakeAperCat : the star at " << Point(s) 
	     << " has " << s.apers.size() << " apertures instead of " << nrads << endl;
	if (s.apers.size() == 0)
	  {
	    cout << " dropping this object " << endl;
	    i = apercat.erase(i);
	  }
	else
	  {
	    cout << " aborting ! " << endl;
	    ++i;
	    abort();
	  }
      } // end loop on checking number of aper meaurements.


  // write some global stuff (radius and seeing, and neighbor cut)
  
  if (!apercat.empty())
    {
      vector<double> rads;
      const vector<Aperture> &a = apercat.front()->apers;
      unsigned naper = a.size();
      for (unsigned k =0; k < naper; ++k) rads.push_back(a[k].radius);
            GlobalVal &glob = apercat.GlobVal();
      glob.AddKey("RADIUS",rads);
      glob.AddKey("SEEING",seeing);
      vector<double> sh; 
      sh.push_back(xSize); sh.push_back(ySize); sh.push_back(corr);
      glob.AddKey("STARSHAPE", sh);
      glob.AddKey("MAXNEIGHBOR_DIST",maxNeighborDist);
    }

  apercat.write(aperCatName);
  return true;
}

		      
bool ReducedImage::MakeStarCat()
{
  if (!MakeAperCat()) 
    {
      cout << " MakeStarList : miss aper catalog for image " << Name() << endl;
      return false;
    }
  string starCatName = StarCatalogName();
  if (FileExists(starCatName)) return true;

  double xSize, ySize, corr;
  StarScoresList scores;
  AperSEStarList apercat(AperCatalogName());
  if (!FindStarShapes(apercat, 20, xSize, ySize, corr, scores))
    {
      cout << " MakeStarList : could not find the star locus for " 
	   << Name() <<endl;
      return false;
    }

  DataCards cards(DefaultDatacards());
  double sigCut = 5;
  if (cards.HasKey("STARLIST_SIG_CUT")) sigCut = 
    cards.DParam("STARLIST_SIG_CUT");

  // recycling apercat magically copies the global keys to output.
  apercat.clear();
  for (StarScoresCIterator i = scores.begin(); i != scores.end(); ++i)
    {
      const StarScores &scores = **i;
      if (scores.nSig<sigCut) // it lies within the star cluster
	{
	  const AperSEStar &a = dynamic_cast< const AperSEStar &>(*scores.star);
	  apercat.push_back(&a);
	}
    }
  GlobalVal &glob = apercat.GlobVal();
  glob.AddKey("STARSIGCUT", sigCut);


  cout << "MakeStarCat() : writing " << starCatName << endl;
  apercat.write(starCatName);

  return true;
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


/* what makes MakeWeight complicated is that it is the product
   of many factors :
   - sky variance
   - flat variations
   - satellites (or planes)
   - dead pixels
   - cosmics
   - (NOT satur, which is handled separately)

  this would be simple if  it could be done right away. But :
  - The catalog routine makes use of weight.
  - cosmic finding requires the seeing (got from catalog)
  - the sky variance is best estimated after sky subtraction done
    while making the catalog

  and perhaps other constraints. 

  So the choosen approach is to write in the weight image what was done
  and to make a routine that updates according to what is missing
  and what is available
*/ 



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
  bool addlow = true;
  double oldvariance = 1 ;
  double newvariance = 1./sq(SigmaBack());
  if (HasWeight())
    {
      if (true)
	// to keep FitsHeader of weight local, 
	// avoiding conflict opening of the same file in RW
	{
	  FitsHeader head(FitsWeightName());
	  adddead = !(head.HasKey("DEADPIXS") && head.KeyVal("DEADPIXS")) 
	    && HasDead();
	  addcosmic = !(head.HasKey("COSMPIXS") && head.KeyVal("COSMPIXS"))
	    && HasCosmic();
	  addflat = !(head.HasKey("FLATPIXS") && head.KeyVal("FLATPIXS")) 
	    && HasFlat() ;
	  addsatellite = !(head.HasKey("SATEPIXS") && head.KeyVal("SATEPIXS"))
	    && HasSatellite();
	  // to cut on too low pixels ,we wait for background subtraction:
	  addlow = !(head.HasKey("LOWPIXS") && head.KeyVal("LOWPIXS"))
	    && BackSub();

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
      if ( !updatevar && !adddead && !addcosmic && !addflat && 
	   !addsatellite && !addlow) 
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
      FitsHeader head(FitsName());
      if (!head.IsValid())
	{
	  cout << " ReducedImage::MakeWeight : cannot Make Weight without original image ... " << endl;
	  return false;
	}
      pweights = new FitsImage(FitsWeightName(), head); 
      cout << " Creating  weight map :" << pweights->FileName() << endl;
      // impose that zeros are preserved after FITS write/read
      pweights->PreserveZeros(); 
      *pweights += 1.; 
    }


  // use auto_ptr to trigger the call to FitsImage destructor
  auto_ptr<FitsImage> no_use(pweights);

  FitsImage &weights = *pweights; // does nothing!
  weights.EnableWrite(false); // prevent writing of partial files.

  /* kill pixels which are too low (usually unidentified bad columns)
     or central area of saturated stars (on some instruments) */
  if (addlow)
    {
      FitsHeader head(FitsName()); // important: don't open RW here
      // wait for image to be back subtracted:
      if (head.HasKey("BACK_SUB") && bool(head.KeyVal("BACK_SUB")))
	{
	  FitsImage image(FitsName()); // important: don't open RW here
	  Pixel mean,sigma;
	  image.MeanSigmaValue(&mean, & sigma);
	  double threshold = mean - 10. * sigma;
	  // only alters the memory image as long as we don't write it.
	  image.Simplify(threshold,0,1);
	  int nlow = int(image.SumPixels());
	  weights *= (1-image);
	  {
	    std::string lowName = Dir()+"/low.fits.gz";
	    remove(lowName.c_str());
	    FitsImage low(lowName, (const Image &) image);
	    low.AddOrModKey("BITPIX",8);
	    CopyWCS(image,low); // for MIBI to display correctly
	  }
	  weights.AddOrModKey("LOWPIXS",true, 
			      "zeroed weight of pixels too low in image");
	  weights.AddOrModKey("LOWTHRE", threshold, 
			      " under this (actual) image value, kill pixel");
	  std::cout << " zeroing weight of low (" << nlow 
		    << ") image pixels in " 
		    << FitsWeightName()  << std::endl;
	}
    }
      
  // check if we have dead, cosmics and flat frames and use them.
  // satur is out of the game because we wish to keep it separate.

  if (updatevar)
    {
      //fill with the inverse of the sky variance
      double s = newvariance / oldvariance ;
      weights *= s ; 
      weights.AddOrModKey("VARPIXS",true," this weight accounts for sky variance");
      weights.AddOrModKey("INVERVAR",newvariance );
      cout << "accounting for sky variance in  " << FitsWeightName() 
	   << " old, new variance^-1: " <<  oldvariance << " " 
	   << newvariance << endl;
    }
  if (adddead && HasDead())
    {
      FitsImage dead(FitsDeadName());
      weights *= (1.- dead);
      weights.AddOrModKey("DEADPIXS",true, " this weight accounts for dead pixels");
      cout << " zeroing dead pixels in " << FitsWeightName() << endl;
    }
  if (addcosmic && HasCosmic())
    {
      FitsImage cosmic(FitsCosmicName());
      weights *= (1.- cosmic);
      weights.AddOrModKey("COSMPIXS",true," this weight accounts for (identified) cosmics");
      cout << " zeroing cosmic pixels in " << FitsWeightName() << endl;
    }

  if (addsatellite && HasSatellite())
    {
      FitsImage satellite(FitsSatelliteName());
      weights *= (1.- satellite);
      weights.AddOrModKey("SATEPIXS",true, " this weight accounts for satelitte tracks");
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
      /* Next line is a BUG !!! : it should read weights *= flat; */
      weights.MultiplyBySquare(flat) ;
      weights.AddOrModKey("FLATPIXS",true," this weight accounts for flat variations");
      cout << " accounting for flat variations in " << FitsWeightName() << endl;

    }

  // end of weight ingredients arithmetics


  // print some statistics
  {
    FitsImage im(FitsName());
    if (im.IsValid())
      cout << Name() << " Image/Weight stats : " 
	   << ImageAndWeightError(im,weights) << endl;
  }
  weights.EnableWrite(true);
  return(true);
}


bool ReducedImage::MakeBad()
{ 
  if (MakeWeight())
    {
      cout << " making BadImage " << FitsBadName() << endl;
      // use "slicing" code in order to accomodate large images:
      FitsInOutParallelSlices inOut(FitsWeightName(),
				    FitsBadName());
      // BITPIX, BSCALE, BZERO are copied by default from in to out
      // So overwrite them here:
      inOut.out.AddOrModKey("BITPIX",8);
      inOut.out.AddOrModKey("BSCALE",1); 
      inOut.out.AddOrModKey("BZERO",0);

      do
	{
	  Image &out = inOut.out; // handler
	  out = inOut.in; // copy
	  out.Simplify(0,0,1); // set to 0 what is > 0 and to 1 otherwise
	} while (inOut.LoadNextSlice());      
      return true ;
    }
  else
    return false ;
}

bool ReducedImage::MakeSatur()
{
  if (HasSatur()) return true;
  cerr << "  ReducedImage::MakeSatur() should in principle never be called .... " << endl;
  std::cerr << " it was called for image " << Name () << std::endl;
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
	if (cosmic && (cosmic->Distance(where) < dist) && !cosmic->IsCosmic())
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
    if (!head.IsValid())
      {
	cout << " cannot Make cosmics without Image " << endl;
	return false;
      }
    if ("Swarp" == string(head.KeyVal("TOADINST")))
      {
	{
	    cout << " image " << head.FileName() 
		 << " was assembled by Swarp : it already has cosmics in the weights" << endl;
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
  if (!MakeWeight()) return false;
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
  if (ToDo & DoAperCatalog) status &= MakeAperCat() && MakeStarCat();
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
      return (header.IsValid() && header.HasKey(KeyName));
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

// seeing from gaussian fits to the objects+ star clump finding 
// (MakeAperCat routine)
REMOVE_ROUTINE(GFSeeing,"GFSEEING") ;
SET_ROUTINE(GFSeeing,double,"GFSEEING") ;

double ReducedImage::GFSeeing() const
{
  FitsHeader head(FitsName());
  if (head.HasKey("GFSEEING")) return double(head.KeyVal("GFSEEING"));
  GlobalVal glob(AperCatalogName());
  if (glob.HasKey("SEEING")) return glob.getDoubleValue("SEEING");
  throw PolokaException(" GFSeeing requested but absent both from fit image and apercat :"+Name());
}


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
      if (head.HasKey("BACK_SUB")) return bool(head.KeyVal("BACK_SUB"));
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

double ReducedImage::BackLevelNoSub() const
{
  FitsHeader head(FitsName());
  if (head.HasKey("BACKLEV")) return double(head.KeyVal("BACKLEV")); // set by transformedimage and such
  if (head.HasKey("SEXSKY")) return double(head.KeyVal("SEXSKY"));   // set when making the catalog.
  if (head.HasKey("SKYLEV")) return double(head.KeyVal("SKYLEV"));   // set when flatfielding
  cout << " ERROR : no way to figure out BackLevel in " << Name() << endl;  return 0;
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
      cerr << " no information about SigmaBack in header of " << Name() << ", compute it" << endl;
      {
	FitsImage image(FitsName());
	Pixel sky,sig;
	image.SkyLevel(&sky, &sig);
	return sig;
      }
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
      if (head.HasKey("SKYLEV"))
	return read_double_key("SKYLEV","OriginalSkyLevel");
      else if (head.HasKey("SEXSKY"))
	return read_double_key("SEXSKY","OriginalSkyLevel");
      else
	return BackLevel();
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
  throw(PolokaException("ReducedImage::Saturation() : cannot find any saturation indication in fits header " ));
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

double ReducedImage::ModifiedModifiedJulianDate() const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader head(FitsName());
      return ModifiedModifiedJulianDay(head);
    }
  else
    {
      cerr << " ReducedImage::ModifiedModifiedJulianDate() cannot compute JulianDay for file " 
	   << fileName << endl;
    }
  return UNDEFINED;
}

double ReducedImage::ModifiedJulianDate() const
{
  string fileName = FitsName();
  if (FileExists(fileName))
    {
      FitsHeader head(FitsName());
      return head.KeyVal("TOADMJD");
    }
  else
    {
      cerr << " ReducedImage::ModifiedJulianDate() cannot compute JulianDay for file " 
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
  FitsHeader head(FitsName());
  FitsHeader otherhead(FitsName());
  return head.SameChipFilterInst(otherhead, Warn);
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



double ReducedImage::OverlapArcmin2(const ReducedImage& Other) const
{
  FitsHeader head(FitsName());
  FitsHeader otherhead(Other.FitsName());
 return Arcmin2Overlap(head, otherhead);
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
  FitsHeader fhead(firstName);
  FitsImage out(OutFitsName, fhead, *sum);
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

  FitsHeader fhead(firstName);
  FitsImage out(OutFitsName, fhead, *sum);
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
  cout << " Common frame for ReduceImageList: " << frame_common << endl;
  return frame_common;
}


ReducedImageCIterator ReducedImageList::Find(const string &Name) const
{
  ReducedImageCIterator i=begin();
  for ( ; i != end(); ++i)
    if ((*i)->Name() == Name) break;
  return i;
}


//ReducedImageRef ReducedImageRead(const char *Name)
//{
//  return ReducedImageRead(string(Name));
//}
	



  



