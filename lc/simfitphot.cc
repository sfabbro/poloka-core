#include <algorithm> // min_element

#include <fastfinder.h>
#include <reducedutils.h>
#include <fstream>
#include "simfit.h"
#include "simfitphot.h"


SimFitPhot::SimFitPhot(LightCurveList& Fiducials,bool usegal)
{
#ifdef FNAME
  cout << " > SimFitPhot::SimFitPhot(LightCurveList& Fiducials)" << endl;
#endif
  //init_phot(Fiducials);
#ifdef DEBUG
  cout << "  zeFit.reserve ... " << endl;
#endif
  dowrite=true;
  zeFit.reserve(Fiducials.Images.size());
  
  zeFit.VignetRef = new SimFitRefVignet(Fiducials.RefImage,usegal); //  no data will be read cause no star is defined
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
  
  cout << "### LCREFTYPE=" << Lc.Ref->type << endl;
  
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
  
  if(dowrite && !IsDirectory(dir))
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
     for(int i=0;i<3;i++) {
      zeFit.SetWhatToFit(FitFlux); // then flux
      zeFit.UseGalaxyModel(true);  
      zeFit.DoTheFit(2);
      if(i==-12) {
	ofstream lstream((string(dir+"/lc_init0.dat")).c_str());
	Lc.write_lc2fit((ostream&)lstream);
	lstream.close();
	zeFit.write("sn_init0",dir, WriteWeight|WriteData);
      }
      zeFit.SetWhatToFit(FitPos); // then pos
      zeFit.UseGalaxyModel(true);  
      zeFit.DoTheFit(2);
     }
    ofstream lstream((string(dir+"/lc_init.dat")).c_str());
    Lc.write_lc2fit((ostream&)lstream);
    lstream.close();
    //return;
#ifdef DEBUG
    cout << " ============= SimFitPhot::operator() Now FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky); // then everything    
     zeFit.DoTheFit(30);     
#ifdef DEBUG
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
       (*itVig)->CheckWeight();
     }
#ifdef DEBUG
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky);
     zeFit.DoTheFit(30);
  }else{
    zeFit.DoTheFit(); 
  }
  if(dowrite) {
  zeFit.write("sn",dir, WriteLightCurve|WriteGalaxy|WriteResid|WriteWeight|WriteData|WriteVignetsInfo|WriteMatrices);
  ofstream lstream((string(dir+"/lc2fit.dat")).c_str());
  Lc.write_lc2fit((ostream&)lstream);
  lstream.close();
  Lc.write_xml((string(dir+"/lc.xml")).c_str()); 
  }
}



