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

#define DEBUG0
void SimFitPhot::operator() (LightCurve& Lc)
{

#ifdef DEBUG0
  switch (Lc.Ref->type)
    {
    case -1:  cout << " ============= SimFitPhot::operator() Fitting a Sn + Galaxy with FIXED position \"" << Lc.Ref->name << "\" =============" << endl; break;
    case 0: cout << " ============= SimFitPhot::operator() Fitting a Sn + Galaxy \"" << Lc.Ref->name << "\" =============" << endl; break;
    case 1: cout << " ============= SimFitPhot::operator() Fitting a star (without galaxy) \"" << Lc.Ref->name << "\"  =============" << endl; break;
    case 2: cout << " ============= SimFitPhot::operator() Fitting a galaxy (without star) \"" <<  Lc.Ref->name << "\"  =============" << endl; break;
    default: cout << " SimFitPhot::operator() : Error : unknown star type :" << Lc.Ref->type << endl; return;
    }
#endif

#ifdef DEBUG0
  cout << " ============= SimFitPhot::operator() zeFit.Load =============" << endl;
#endif  
  zeFit.Load(Lc);
  
  switch (Lc.Ref->type)
    {
    case -1: // sn+galaxy fixed position 
      break; 
    case 0: // sn+galaxy
      break;
    case 1: //star 
      zeFit.SetWhatToFit(FitFlux);
      zeFit.UseGalaxyModel(false);
      zeFit.DoTheFit();
      zeFit.SetWhatToFit(FitFlux  | FitPos | FitSky );
      zeFit.UseGalaxyModel(false);
      zeFit.DoTheFit(50,0.05);
      // robustify to get rid of other stars in the vignet
      for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
	(*itVig)->KillOutliers();
	(*itVig)->CheckWeight();
      }
      zeFit.SetWhatToFit(FitFlux  | FitPos | FitSky );
      zeFit.UseGalaxyModel(false);
      //zeFit.DoTheFit(50,0.05);
      break;
    case 2: //galaxy 
      zeFit.SetWhatToFit(                   FitGal); 
      break;
    }
  
  //string dir = Lc.Ref->name;
  string dir = ".";
  
  if(dowrite && !IsDirectory(dir))
    MKDir(dir.c_str());
  
  if(Lc.Ref->type == 0) {
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() FitInitialGalaxy =============" << endl;
#endif
    zeFit.FitInitialGalaxy(); // first fit inital galaxy
    zeFit.write("sn_init",dir,WriteGalaxy);
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() Fit Flux then Position 3 times =============" << endl;
#endif
     for(int i=0;i<3;i++) {
      zeFit.SetWhatToFit(FitFlux); // then flux
      zeFit.UseGalaxyModel(true);  
      zeFit.DoTheFit(2,0.05);
      zeFit.SetWhatToFit(FitPos); // then pos
      zeFit.UseGalaxyModel(true);  
      zeFit.DoTheFit(2,0.05);
     }
    ofstream lstream((string(dir+"/lc_init.dat")).c_str());
    Lc.write_lc2fit((ostream&)lstream);
    lstream.close();
    //return;
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() Now FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky); // then everything    
     zeFit.DoTheFit(30,0.05);     
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
       (*itVig)->CheckWeight();
     }
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitPos | FitGal | FitSky =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky);
     zeFit.DoTheFit(30,0.005);
  }

  //============================================================
  //============================================================

  
  if(Lc.Ref->type == -1) {
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() FitInitialGalaxy =============" << endl;
#endif
    zeFit.FitInitialGalaxy(); // first fit inital galaxy
    zeFit.write("sn_init",dir,WriteGalaxy);
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() Fit Flux then Position 3 times =============" << endl;
#endif
    zeFit.SetWhatToFit(FitFlux); // then flux
    zeFit.UseGalaxyModel(true);  
    zeFit.DoTheFit(2,0.05);
    ofstream lstream((string(dir+"/lc_init.dat")).c_str());
    Lc.write_lc2fit((ostream&)lstream);
    lstream.close();
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() Now FitFlux | FitGal | FitSky =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitSky); // then everything    
     zeFit.DoTheFit(30,0.05);     
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
       (*itVig)->CheckWeight();
     }
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitPos | FitGal =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitSky);
     zeFit.DoTheFit(30,0.005);
  }
  
  
  if(Lc.Ref->type>0){
    zeFit.DoTheFit(50,0.005); 
  }

  if(dowrite) {
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() Write =============" << endl;
#endif
  zeFit.write("sn",dir, WriteLightCurve|WriteGalaxy|WriteResid|WriteWeight|WriteData|WriteVignetsInfo|WriteMatrices);
  ofstream lstream((string(dir+"/lc2fit.dat")).c_str());
  Lc.write_lc2fit((ostream&)lstream);
  lstream.close();
  Lc.write_xml((string(dir+"/lc.xml")).c_str()); 
  }
}



