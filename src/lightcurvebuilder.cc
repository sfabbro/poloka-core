#include <iomanip>
#include <fstream>
#include "lightcurvebuilder.h"
#include "simultaneousfit.h"
#include "fitsimage.h"
#include "kernelfit.h"
#include "vignet.h"
#include "vignetfit.h"
#include "reducedutils.h"
#include "imagesubtraction.h"
#include "sestar.h"
#include "jimstar.h"


static double sqr(const double &x) {return x*x;}

static void OptimizeVignetSizes(const double &BestSeeing, const double &WorstSeeing, 
				int &HrefX, int &HrefY)
{
  // default is 1 FWHM of the worst seeing
  double nSigPsf = 2.3548;
  if (getenv("NSIGPSF")) nSigPsf = atof(getenv("NSIGPSF"));
  int hsize = int(ceil(WorstSeeing*nSigPsf));
  OptParams kernelParams;
  cout << " Worst seeing image half size: " << hsize << endl;
  kernelParams.OptimizeSizes(BestSeeing, WorstSeeing);
  hsize = max(hsize, kernelParams.HKernelSize) + kernelParams.HKernelSize;
  cout << " Chosen max half vignet size is " << hsize << endl;
  HrefX = HrefY = hsize;
}

LightCurveBuilder::LightCurveBuilder(const ReducedImageList &Images)
  : Reference(NULL), Nights(Images), BandName(Images.front()->Band())
{
  cout << "**************************" << endl;
  cout << "   EN VOITURE SIMONE      " << endl;
  cout << "**************************" << endl;
}


LightCurveBuilder::LightCurveBuilder(const LightCurveBuilder &Other) 
  : vector<LightCurve> (), Reference(NULL)
{
  *this = Other;
}

LightCurveBuilder& LightCurveBuilder::operator = (const LightCurveBuilder &Right)
{  
  // a real pain to do a copy constructor again.
  // when people stop put pointers in classes, we'll be all bishops
  *((vector<LightCurve> *) this)  =  Right; 
  if (Reference) delete Reference;
  if (Right.Reference) Reference = new Night(Right.Reference->Name());
  Nights = Right.Nights; 
  for (unsigned i=0; i<Nights.size(); ++i)
    {
      const Night *night = Nights[i];
      for (LightCurveIterator fid = begin(); fid != end(); ++fid)
	(*fid)[i]->night = night;
    }
  BandName = Right.BandName;
  nsn = Right.nsn;
  return *this;
}

LightCurveBuilder::~LightCurveBuilder()
{
  for (NightCIterator it = Nights.begin(); it != Nights.end(); ++it)
    {
      const Night *night = *it;
      for (size_t i=1; i<size(); ++i)
	{	  
	  string toremove = night->Name()+"_"+(*this)[i].FidName;
	  remove((toremove+"_data.fits").c_str());
	  remove((toremove+"_kernel.fits").c_str());
	  remove((toremove+"_weight.fits").c_str());
	}
    }

  if (Reference) delete Reference;
  cout << "**************************" << endl;
  cout << "      TERMINUS            " << endl;
  cout << "**************************" << endl;
}

void LightCurveBuilder::coincid()
{
  // to make sure we assign the right nights in right order of the fiducials
  sort(Nights.begin(), Nights.end(), IncreasingNightSeeing);
  for (LightCurveIterator fid = begin(); fid != end(); ++fid)
    sort(fid->begin(), fid->end(), ByIncreasingSeeing);
}

void LightCurveBuilder::InitFiducials(const BaseStarList &FidList)
{
  cout << " Now Initialize Fiducials LightCurves" << endl;

  // remember Fiducials[nsn] are the SNe
  unsigned count = 0;
  for(BaseStarCIterator bsi = FidList.begin(); bsi != FidList.end(); ++bsi, ++count)
    {
      LightCurve fid(Nights, **bsi);
      char toto[12];
      if (count < nsn) sprintf(toto,"sn%d", count);
      else sprintf(toto,"fid%d", count-nsn);
      fid.FidName = toto;
      push_back(fid);
    }
  
  if (count < 1) {cerr << " Did not select any fiducials !! " << endl;} 
  else cout << " " << count << " fiducials were initialized"  << endl;
  StarListPhotometry();
}

void LightCurveBuilder::BuildKernelsAndSubs()
{
  cout << "**************************" << endl;
  cout << "   TCHOO TCHOO            " << endl;
  cout << "**************************" << endl;
  cout << " Building kernels and subtractions" << endl;

  coincid();
  int hrefx,hrefy;
  OptimizeVignetSizes(Reference->seeing, Nights.back()->seeing, hrefx, hrefy);  

  for (size_t i=0; i<Nights.size(); ++i)
    {
      cout << "--------------------------------------------" << endl;
      Night *night = Nights[i];
      string subName = SubtractedName(Reference->Name(), night->Name());
      ImageSubtraction sub(subName, *Reference, *night);
      
      // check if subtractions and kernel images  are produced
      bool is_produced = true;
      for (LightCurveCIterator fid = begin(); fid != end(); ++fid)
	if (!FileExists(night->Name()+"_"+fid->FidName+"_kernel.fits")
	    || !FileExists(night->Name()+"_"+fid->FidName+"_data.fits")) is_produced=false;
      if (!FileExists(sub.FitsName())) is_produced = false;
      if (is_produced) continue;

      // if not then produce them
      sub.FitKernel();
      sub.Execute(DoFits);
      cout << " Old photometric ratio = " << night->photomRatio << endl;
      double chi2= sub.Chi2();
      double ratio= sub.PhotomRatio();
      if (*night == *Reference)
	{
	  chi2 = 0;
	  ratio = 1;
	}
      night->photomRatio = ratio;
      night->kernelChi2 = chi2;
      cout << " New photometric ratio = " << night->photomRatio << endl;

      // build vignets but don't keep in memory local kernel of night i
      FitsImage subim(sub.FitsName());
      FitsImage im(night->FitsName(),RW);
      FitsImage weight(night->FitsWeightName());
      
      im.AddOrModKey("KERNREF", Reference->Name().c_str(), " name of the seeing reference image");
      im.AddOrModKey("KERNCHI2",chi2, " chi2 of the kernel fit");
      im.AddOrModKey("PHORATIO",ratio, " photometric ratio with KERNREF");
      im.AddOrModKey("ISZERO",int(night->IsZeroRef)," is image zero reference");

      // loop over fiducials within the current image
      int count = 0;
      for (LightCurveIterator fid = begin(); fid != end(); ++fid, ++count)
	{
	  FiducialStar *fs = (*fid)[i];
	  // write on disk kernel and data vignets to save time during simultaneous fit
    	  //if (!SameSeeing(*Reference,*night))
	  
	  
	  Vignet vig(fs->x, fs->y, im, hrefx, hrefy);
	  vig.writeFits(night->Name()+"_"+fid->FidName+"_data.fits");


	  if (Reference->Name() != night->Name()) {
	      Kernel kern;
	      sub.KernelToWorst(kern, fs->x, fs->y);
	      kern.writeFits(night->Name()+"_"+fid->FidName+"_kernel.fits");
	  }
	  
	  Vignet w(fs->x, fs->y, weight, hrefx, hrefy);
	  w.writeFits(night->Name()+"_"+fid->FidName+"_weight.fits");
	}
    }

  FitsImage refim(Reference->FitsName());
  for (LightCurveCIterator fid = begin(); fid != end(); ++fid)
    {
      Vignet refvig(fid->refstar.x, fid->refstar.y, refim, hrefx, hrefy);
      refvig.writeFits(Reference->Name()+"_"+fid->FidName+"_data.fits");
    }
}

void LightCurveBuilder::StarListPhotometry()
{
  cout << " Perform quick photometry on fiducials using existing starlist's" << endl;

  coincid();
  
  // loop over nights
  for (size_t i=0; i<Nights.size(); ++i)
    {
      Night *night = Nights[i];
      SEStarList curList(night->CatalogName());
      // loop over stars
      for (LightCurveIterator fid = begin(); fid != end(); ++fid)
	{
	  FiducialStar *fs = (*fid)[i];	  
	  SEStar *star = curList.FindClosest(*fs);
	  *((PhotStar*) fs)  = PhotStar(*star);
	  // approximate (since we don't have any) error on positioning for
	  // future ponderation
	  fs->varx = fs->varflux * sqr(night->sigmaX / fs->flux);
	  fs->vary = fs->varflux * sqr(night->sigmaY/fs->flux);
	  fs->covxy = fs->varflux * sqr(night->thetaXY/fs->flux);
	}
    }  
  // a hack: restore sn position because usually the starlist 
  // has worse initial estimate than the subtraction where SN was detected
  for (size_t i=0; i<Nights.size(); ++i)
    {
      for (LightCurveIterator fid = begin(); fid != end(); ++fid)
	{
	  FiducialStar *fs = (*fid)[i];	  
	  fs->x = fid->refstar.x;
	  fs->y = fid->refstar.y;
	}
    }
  for (LightCurveIterator fid = begin(); fid != end(); ++fid) fid->RelativeFluxes();
  RefStarListPhotometry();
}

void LightCurveBuilder::RefStarListPhotometry()
{
  cout << " Performing starlist photometry on reference " << endl;

  // loop over fiducials
  for (LightCurveIterator fid = begin(); fid != end(); ++fid)
    {
      fid->refstar.flux = 0;
      fid->refstar.sky = 0;
      fid->refstar.varflux = 0.1;
      fid->refstar.varsky = 0;
    } 
  // end loop over fiducials
}

void LightCurveBuilder::SubAperPhotometry(const double &Nfwhm)
{
  cout << " Performing aperture photometry on subtractions " << endl;

  // loop over nights
  for (unsigned i=0; i<Nights.size(); ++i)
    {
      Night *night = Nights[i];
      string subName = SubtractedName(Reference->Name(), night->Name());
      Night subnight(subName);
      FitsImage subIm(subnight.FitsName());
      double radius = subnight.seeing * 2.3548 * Nfwhm;
      int halfvigsize = int(radius*3);
      unsigned count = 0;
      double sigma2 = sqr(subnight.sigmaSky);
      // loop over fiducials
      for (LightCurveIterator fid = begin(); fid != end(); ++fid, ++count)
	{
	  FiducialStar *fs = (*fid)[i];
	  //Vignet subvig(fs->x, fs->y, subIm, halfvigsize, halfvigsize);
	  Kernel toto;
	  VignetFit subvig(fs, &subnight, subIm, toto, halfvigsize, halfvigsize);
	  // save the first two for posterity
	  if (count<2) subvig.writeFits(night->Name()+"_"+fid->FidName+"_sub.fits");
	  if (night->Name() != Reference->Name())
	    {
	      fs->flux = subvig.Aperture(fs->varflux, radius, sigma2, 0);
	      // fs->flux = subvig.WeightedAperture(fs->varflux);
	      // subtracted image is in unit of the reference!
	      fs->flux *= night->photomRatio;
	      fs->varflux *= sqr(night->photomRatio);
	    }
	  else 
	    {
	      fs->flux = 0;
	      fs->varflux = 0.1;
	    }
	} // end loop over fiducials
    } // end loop over nights

  for (LightCurveIterator fid = begin(); fid != end(); ++fid) fid->RelativeFluxes();

  RefAperPhotometry(Nfwhm);
  
}

void LightCurveBuilder::RefAperPhotometry(const double &Nfwhm)
{
  cout << " Performing aperture photometry on reference " << endl;

  FitsImage refIm(Reference->FitsName());
  // loop over fiducials
  double sigma2 = sqr(Reference->sigmaSky);

  for (LightCurveIterator fid = begin(); fid != end(); ++fid)
    {
      double radius = Reference->seeing * 2.3548 * Nfwhm;
      int halfvigsize = int(radius*5);
      Vignet vig(fid->refstar.x, fid->refstar.y, refIm, halfvigsize, halfvigsize);
      //fid->refstar.flux = subvig.WeightedAperture(fid->refstar.varflux);      
      fid->refstar.flux = vig.Aperture(fid->refstar.varflux, radius, sigma2, Reference->backLevel);      
    } 
  // end loop over fiducials
}

void LightCurveBuilder::SimFitPhotometry()
{
  cout << " **************************" << endl;
  cout << "    TATAC TATOUM           " << endl;
  cout << " **************************" << endl;
  cout << " Fit vignets simultaneously" << endl;

  coincid();

  // loop over lightcurves
  unsigned count = 0;
  for (LightCurveIterator it = begin(); it != end(); ++it, ++count)
    {
      cout <<"--------------------------------------------------" << endl;
      cout << " Starting fiducial #" << count  << endl;
      cout <<"--------------------------------------------------" << endl;

      SimultaneousFit zeFit;
      bool is_sn = (count < nsn);
      // loop over nights
      // zeFit.VignetRef = new VignetFit(&rstar, it->refstar, it->FidName);
      // zeFit.VignetRef->IsRefResolution = true;
      // zeFit.push_back(refVig);      
      for (size_t i=0; i<Nights.size(); ++i)
	{
	  FiducialStar *fs = (*it)[i];
	  VignetFit *vignet = new VignetFit(fs, fs->night, it->FidName);
	  vignet->flux /= fs->night->photomRatio;
	  vignet->IsStarHere = !(fs->night->IsZeroRef);
	  //vignet->IsRefResolution = SameSeeing(*Reference, *fs->night);
	  vignet->IsRefResolution = (Reference->Name() ==  fs->night->Name());
	  vignet->Weight = Kernel(fs->night->Name()+"_"+it->FidName+"_weight.fits");
	  zeFit.push_back(vignet);
	}
      zeFit.VignetRef = zeFit.front();
      // do not refit position if it's not a variable object: initial flux estimates
      // are derived from differential fluxes, which are supposedly close to zero.
      int fitsky = 0;
      if (getenv("SIMFITSKY")) fitsky = FitSky;
      if (is_sn) zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | fitsky);
      else zeFit.SetWhatToFit(FitFlux | FitGal | fitsky);
      zeFit.FindMinimumScale(Nights.back()->seeing);
      //zeFit.Minim = LevenbergMarquardt;
      zeFit.DoTheFit();
      it->refstar.flux = zeFit.GetGalaxyFlux(it->refstar.varflux);
      // write things down if it was the supernova, so we have something saved
      // before all the loops on the fiducials
      if (is_sn)
	{
	  string midname = it->FidName+"_"+BandName;
	  for(VignetFitCIterator vi = zeFit.begin();vi != zeFit.end(); ++vi) 
	    {
	      (*vi)->writeAllFits(it->FidName);
	    }
	  
	  zeFit.write(midname);
	  string outstr = "lightcurve_"+midname+".dat";
	  ofstream out(outstr.c_str());
	  it->dump(out);
	  it->dump();
	  //if (is_sn) break;
	}
      // delete zeFit.VignetRef; //to comment if you want to fit VignetRef as well
    }
}

void LightCurveBuilder::Calibrate(const double &Nfwhm) const
{
  for (unsigned i=0; i<Nights.size(); ++i)
    {
      const Night *night = Nights[i];
      //if (!night->IsPhotometric) continue;
      double radius = night->seeing * 2.3548 * Nfwhm;
      int hsize = int(ceil(radius+5));
      double sigma2 = sqr(night->sigmaSky);
      FitsImage im(night->FitsName());
      cout << " Full aperture photometry on standards of " << night->Name() << endl;
      int count = 0;
      ofstream out(string(night->Dir()+"/standard.list").c_str());
      out << "# n : star number" << endl;
      PhotStar().WriteHeader(out);
      for (LightCurveCIterator it = begin(); it != end(); ++it,++count)
	{
	  const FiducialStar *fs = (*it)[i];
	  PhotStar standard(*fs);
	  Vignet vig(standard.x, standard.y, im, hsize, hsize);
	  standard.flux = vig.Aperture(standard.varflux, radius, sigma2, standard.sky);
	  out << count << ' ';
	  standard.writen(out);
	  out << endl;
	}
      out.close();
    }
}

bool LightCurveBuilder::Init(const vector<Point> &SNe)
{
  // init the night list
  if (!Nights.Init()) return false;
  if (!Reference) 
    {
      sort(Nights.begin(), Nights.end(), IncreasingSeeing);
      Reference = new Night(*Nights.front());
    }
  Nights.InitPhotomRatio(*Reference);
	
  // fiducials[0:nsn-1] are the supernovae
  nsn = SNe.size();
  for (unsigned i=0;i<nsn;++i)
    {
      BaseStar snb(SNe[i],0);
      LightCurve snlc(Nights, snb);
      for (unsigned j=0; j<Nights.size(); ++j) 
	snlc[j]->sky = snlc[j]->night->backLevel;
      char toto[12];
      sprintf(toto,"sn%d", i);
      snlc.FidName = toto;
      push_back(snlc);
    }

  // get all the fiducials of that band to start with
  sort(Nights.begin(),Nights.end(),DecreasingNightRatio);
  SEStarList fidSEList(Reference->CatalogName());
  Frame refframe = Nights.CommonFrame()*Reference->UsablePart();
  
  BaseStarList fidList;
  SE2Base(&fidSEList)->CopyTo(fidList);
  Nights.FilterAllObjects(fidList, refframe);
  InitFiducials(fidList);

  cout << " LightCurveBuilder is initialized " << endl;

  return (size() != 0);
}

bool LightCurveBuilder::Init(const vector<Point> &SNe, const string &catalogue)
{
  // init the night list
  if (!Nights.Init()) return false;
  if (!Reference) 
    {
      sort(Nights.begin(), Nights.end(), IncreasingSeeing);
      Reference = new Night(*Nights.front());
    }
  Nights.InitPhotomRatio(*Reference);
	
  // fiducials[0:nsn-1] are the supernovae
  nsn = SNe.size();
  for (unsigned i=0;i<nsn;++i)
    {
      BaseStar snb(SNe[i],0);
      LightCurve snlc(Nights, snb);
      for (unsigned j=0; j<Nights.size(); ++j) 
	snlc[j]->sky = snlc[j]->night->backLevel;
      char toto[12];
      sprintf(toto,"sn%d", i);
      snlc.FidName = toto;
      push_back(snlc);
    }

  // get all the fiducials of that band to start with
  sort(Nights.begin(),Nights.end(),DecreasingNightRatio);
  
  // reference frame where fiducials are selected
  Frame refframe = Nights.CommonFrame()*Reference->UsablePart();
  
  // get header of the reference image
  FitsHeader header(Reference->FitsImageName(Calibrated));
  if(!header.IsValid()) {
    cerr << "Error in LightCurveBuilder::Init : Reference is not valid" << endl;
    return false;
  }

  // Here we select the catalogue
  BaseStarList *fidList =  (BaseStarList *) GetSelectedJimStarList(header,refframe,catalogue);
  
  // and we init the light curves 
  InitFiducials(*fidList);
  
  cout << " LightCurveBuilder is initialized " << endl;

  return (size() != 0);
}

void LightCurveBuilder::write(const string &MiddleName) const
{
  string fidname = "fiducials_" + MiddleName + "_" + BandName + ".list";

  ofstream fidout(fidname.c_str());
  fidout << "# n : fiducial tag number " << endl;
  fidout << "# jd : reduced julian date " << endl;
  (*this)[0].front()->WriteHeader(fidout);
  fidout << setiosflags(ios::fixed);
  for (size_t i=0; i<Nights.size() ; ++i)
    {
      int ind = 0;
      for (LightCurveCIterator fid = begin(); fid != end(); ++fid, ++ind)
	{
	  const FiducialStar *fs  = (*fid)[i]; // fiducial on x,y and on night i
	  fidout << setiosflags(ios::left);
	  fidout << setw(6) << ind
		 << setw(12) << setprecision(2) << Nights[i]->julianDate;
	  fidout << resetiosflags(ios::left);
	  fs->writen(fidout);
	  fidout << endl;
	}
    }
  fidout << resetiosflags(ios::fixed);
  fidout.close();

  string refname = "refstars_" + MiddleName + "_" + BandName + ".list";
  ofstream refout(refname.c_str());
  refout << setiosflags(ios::fixed);
  refout << "# n : fiducial tag number " << endl;
  refout << "# jd : reduced julian date " << endl;
  (*this)[0].front()->WriteHeader(refout);
  refout << "# stand : standard flag " << endl;
  int ind = 0;
  for (LightCurveCIterator fid = begin(); fid != end(); ++fid, ++ind)
    {
      refout << ind << " " << Reference->julianDate << " ";
      fid->refstar.writen(refout);
      refout << " " << fid->IsStandard << endl;
    }
  refout << resetiosflags(ios::fixed);
}
