#include <algorithm> // min_element

#include <fastfinder.h>
#include <reducedutils.h>
#include <fstream>
#include "simfit.h"
#include "simfitphot.h"



/*
static bool IncSeeing(const CountedRef<ReducedImage> one, const CountedRef<ReducedImage> two)
{ return (one->Seeing() < two->Seeing());}

 
static void init_phot(LightCurveList& Fiducials)
{  
#ifdef FNAME
cout << " > init_phot(LightCurveList& Fiducials)" << endl;
#endif
if (Fiducials.size() == 0)  {
      cerr << " init_phot() : no objects to initialize! \n";
      return;
    }

  // if no given reference, make it as best seeing image
  if (!Fiducials.RefImage) 
    {
      
      Fiducials.RefImage = *min_element(Fiducials.Images.begin(), Fiducials.Images.end(), IncSeeing);

      for_each(Fiducials.Objects.begin(), Fiducials.Objects.end(),
	       bind2nd(mem_fun(&Fiducial<PhotStar>::AssignImage), Fiducials.RefImage));
    }
 
  cout << " init_phot() : initalizing photometry for " << Fiducials.size() << " objects " << endl;

  const double maxDist = 2.;
  const SEStarList sexRefCat(Fiducials.RefImage->CatalogName());
  const BaseStarList *refCat(SE2Base(&sexRefCat));

  for (ReducedImageCIterator it=Fiducials.Images.begin(); it != Fiducials.Images.end(); ++it)
    {
      const ReducedImage *itRed  = *it;
      const SEStarList sexCat(itRed->CatalogName());
      const BaseStarList *cat(SE2Base(&sexCat));
      const FastFinder sexFinder(*cat);
      const Frame thisFrame = itRed->UsablePart();
      const double ratio = MedianPhotomRatio(*cat, *refCat);

      cout << " " << itRed->Name() << " " << Fiducials.RefImage->Name() 
	   << " median ratio = " << ratio << endl;

      for (LightCurveList::iterator si=Fiducials.begin(); si != Fiducials.end(); ++si)
	{
	  if (!thisFrame.InFrame(*si->Ref))
	    {
	      cerr << " init_phot() : Warning : " 
		   << si->Ref->x << " " << si->Ref->y << " is outside of image " << itRed->Name() << endl;
	      continue;
	    }

	  PhotStar *fidPhot = new PhotStar(BaseStar(si->Ref->x, si->Ref->y, 0.));
	  const SEStar *fidSex = (SEStar *) sexFinder.FindClosest(*si->Ref, maxDist);
	  if (fidSex) 
	    { 
	      delete fidPhot;
	      fidPhot = new PhotStar(*fidSex);
	      fidPhot->flux    /= ratio;
	      fidPhot->varflux /= ratio*ratio;
	    }
	  si->push_back(*it, fidPhot);
	}
    }
}
*/
SimFitPhot::SimFitPhot(LightCurveList& Fiducials)
{
#ifdef FNAME
  cout << " > SimFitPhot::SimFitPhot(LightCurveList& Fiducials)" << endl;
#endif
  //init_phot(Fiducials);
#ifdef DEBUG
  cout << "  zeFit.reserve ... " << endl;
#endif
  zeFit.reserve(Fiducials.Images.size());
  
  zeFit.VignetRef = new SimFitRefVignet(Fiducials.RefImage); //  no data will be read cause no star is defined
#ifdef DEBUG
  SimFitRefVignet* toto = zeFit.VignetRef;
  cout << " in SimFitPhot::SimFitPhot,  zeFit.VignetRef     = " << toto << endl;
#endif
  for (ReducedImageCIterator it=Fiducials.Images.begin(); it != Fiducials.Images.end(); ++it)
    {
      SimFitVignet *vig = new SimFitVignet(*it,zeFit.VignetRef);
      zeFit.push_back(vig);
    }
  // cannot do this since we do not have any data loaded yet 
  //zeFit.FindMinimumScale(worstSeeing);
  
}


void SimFitPhot::operator() (LightCurve& Lc)
{
#ifdef DEBUG
    cout << " ============= SimFitPhot::operator() zeFit.Load =============" << endl;
#endif  
  zeFit.Load(Lc);
  
  switch (Lc.Ref->type)
    {
    case 0: // sn+galaxy
      break;
    case 1: //star 
      zeFit.SetWhatToFit(FitFlux);
      zeFit.UseGalaxyModel(false);
      zeFit.DoTheFit();
      zeFit.SetWhatToFit(FitFlux  | FitPos);
      zeFit.UseGalaxyModel(false);
      break;
    case 2: //galaxy 
      zeFit.SetWhatToFit(                   FitGal); 
      break;
    default: cerr << " SimFitPhot::operator() : Error : unknown star type :" 
		  << Lc.Ref->type << endl; return;
    }
  
  string dir = Lc.Ref->name;
  if(!IsDirectory(dir))
    MKDir(dir.c_str());
  
  if(Lc.Ref->type == 0) {
#ifdef DEBUG
    cout << " ============= SimFitPhot::operator() zeFit.FitInitialGalaxy =============" << endl;
#endif
    zeFit.FitInitialGalaxy(); // first fit inital galaxy
    zeFit.write("sn_init",dir,WriteGalaxy);
#ifdef DEBUG
    cout << " ============= SimFitPhot::operator() First FitFlux =============" << endl;
#endif
    zeFit.SetWhatToFit(FitFlux); // then flux
    zeFit.UseGalaxyModel(true);  
    zeFit.DoTheFit();
    //zeFit.write("sn_init1",dir,WriteLightCurve|WriteVignetsInfo|WriteMatrices);
    
#ifdef DEBUG
    cout << " ============= SimFitPhot::operator() Now FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky); // then everything    
     zeFit.DoTheFit();     
#ifdef DEBUG
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
     }
#ifdef DEBUG
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.DoTheFit();
  }else{
    zeFit.DoTheFit(); 
  }
  
  zeFit.write("sn",dir, WriteLightCurve|WriteGalaxy|WriteResid|WriteWeight|WriteData|WriteVignetsInfo|WriteMatrices);
  ofstream lstream((string(dir+"/lc.dat")).c_str());
  Lc.write_short((ostream&)lstream);
  lstream.close();
  Lc.write_xml((string(dir+"/lc.xml")).c_str()); 
  
}



