#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <ctime>

#include "simfitvignet.h"
#include "simfit.h"
 

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  :::::::::::::::::: SimFit stuff   ::::::::::::::::::::::::::
  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

SimFit::SimFit()
{
  refill = true;
  fit_flux = fit_gal = true; fit_sky = fit_pos = false ;
  use_gal = true;
  fluxstart = galstart = skystart = xind = yind = 0;
  fluxend = galend = skyend = 0;
  scale = 1., minscale = 0.; chi2 = 1e29;
  nfx = nfy = hfx = hfy = nparams = ndata = 0;
  dont_use_vignets_with_star = false;
}

void SimFit::UseGalaxyModel(bool useit) {
  use_gal = useit;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    (*itVig)->UseGal = use_gal;
  }
}

void SimFit::SetWhatToFit(unsigned int ToFit)
{
  bool resize = false;
  if (fit_flux != (ToFit & FitFlux)) { fit_flux = (ToFit & FitFlux); resize=true; }
  if (fit_pos  != (ToFit & FitPos))  { fit_pos  = (ToFit & FitPos ); resize=true; }
  if (fit_gal  != (ToFit & FitGal))  { fit_gal  = (ToFit & FitGal ); resize=true; }
  if (fit_sky  != (ToFit & FitSky))  { fit_sky  = (ToFit & FitSky ); resize=true; }
  use_gal |= fit_gal; // use the galaxy if we fit it
  dont_use_vignets_with_star = false;

  int count=0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    (*itVig)->FitPos = fit_pos && (*itVig)->CanFitPos;
    (*itVig)->UseGal = use_gal;
    (*itVig)->FitSky = fit_sky && (*itVig)->CanFitSky;
    (*itVig)->FitFlux = fit_flux && (*itVig)->CanFitFlux;
    cout << count << " vignet whattofit: " 
	 << (*itVig)->FitPos << " "
	 << (*itVig)->UseGal << " "
	 << (*itVig)->FitSky << " "
	 << (*itVig)->FitFlux << endl;
    count++;
  }
  
  cout << " SimFit::SetWhatToFit(" << ToFit << ") : fit_flux " << fit_flux << " fit_pos " << fit_pos 
       << " fit_gal " << fit_gal << " fit_sky " << fit_sky << endl;
}

void SimFit::FindMinimumScale(double WorstSeeing)
{
#ifdef FNAME
  cout << " > SimFit::FindMinimumScale(const double &WorstSeeing) WorstSeeing=" << WorstSeeing  << endl;
#endif

  int hmin = max(int(ceil(WorstSeeing*2.3548)), 5);
  int hrefx = VignetRef->Data.HSizeX();
  int hrefy = VignetRef->Data.HSizeY();
  int hkx = 0;
  int hky = 0;

#ifdef DEBUG
  cout << "VignetRef->Data.HSizeX() = " << hrefx << endl;
  cout << "VignetRef->Data.HSizeY() = " << hrefy << endl;
#endif

  for (SimFitVignetIterator it = begin(); it != end(); ++it)
    {
      SimFitVignet *vi = *it;
      if ((vi->Kern.HSizeX() > hkx) || (vi->Kern.HSizeY()> hky))
	{
	  hkx = vi->Kern.HSizeX();
	  hky = vi->Kern.HSizeY();
	  hrefx = max(hmin, hkx) + hkx;
	  hrefy = max(hmin, hky) + hky;
	}
    }

  minscale = double(min(hrefx,hrefy)) / double(min(VignetRef->Data.HSizeX(),VignetRef->Data.HSizeY()));

#ifdef DEBUG
  cout << " SimFit::FindMinimumScale(" << WorstSeeing 
       << ") : Minimum scaling factor = " << minscale << endl;
#endif

}

double SimFit::GetWorstSeeing() {
#ifdef FNAME
  cout << " > SimFit::GetWorstSeeing(LightCurve& Lc)" << endl;
#endif 
  double worstSeeing = 0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    double currentseeing = (*itVig)->Image()->Seeing();
    if(currentseeing>worstSeeing)
      worstSeeing=currentseeing;
  }
  return worstSeeing;
}

void SimFit::Load(LightCurve& Lc, bool keepstar)
{
#ifdef FNAME
  cout << " > SimFit::Load(LightCurve& Lc)" << endl;
#endif
  
  if (Lc.size() != size()) 
    {
      cerr << " SimFit::Load() : Error : trying to load a LightCurve of different size \n";
      return;
    }
  
  // define the size of the reference vignet
  
  // get worst seeing to set the size of reference vignette
  double worst_seeing = GetWorstSeeing();
  
  
  // now build kernels of vignets to get the worst kernel
  int worst_kernel = 0;
  LightCurve::const_iterator itLc = Lc.begin();
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig, ++itLc)
    {      
      //(*itVig)->Load(*itLc); // this does not modify the kernel so I would better use SetStar which does nothing but set the star
      if(!keepstar) {
	(*itVig)->SetStar(*itLc); // just set the star
	(*itVig)->BuildKernel(); // rebuild kernel
      } 
      if((*itVig)->Kern.HSizeX()>worst_kernel)
	worst_kernel = (*itVig)->Kern.HSizeX();
      if((*itVig)->Kern.HSizeY()>worst_kernel)
	worst_kernel = (*itVig)->Kern.HSizeY();
    }
  
  // 2.3548*sigma = full-width at half-maximum [2.3548 = 2.*sqrt(2*log(2.))]
  // seeing (from sextractor SESEEING) is sigma in pixel units
  
  // radius is the size of the reference vignet
  int radius = int(ceil(2.3548*worst_seeing+worst_kernel)); 
  // minscale  = min_radius/radius (min_radius is used for fitting the position)
  minscale = (worst_seeing+worst_kernel)/radius;
  
#ifdef DEBUG
  cout << " in SimFit::Load, worst_seeing = " << worst_seeing << endl;
  cout << " in SimFit::Load, worst_kernel = " << worst_kernel << endl;
  cout << " in SimFit::Load, Reference vignet radius = " << radius << endl;
  cout << " in SimFit::Load, minscale = " << minscale << endl;
#endif
  
  // the VignetRef has already been build
  if(!keepstar)
    VignetRef->SetStar(Lc.Ref); // just set the star
  VignetRef->Resize(radius,radius); // now resize, this reloads data, update psf, and makeInitialGalaxy if usegal
#ifdef DEBUG
  SimFitRefVignet *toto = VignetRef;
  cout << " in SimFit::Load,  VignetRef     = " << toto << endl;
  cout << " in SimFit::Load,  VignetRef->Hx() = " <<  VignetRef->Hx() << endl;
  cout << " in SimFit::Load,  VignetRef->Hy() = " <<  VignetRef->Hy() << endl;
#endif 

  // now resize all vignets and initialize residus
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig)
    {           
      // we will fit the flux according to Lc.Ref and date
      (*itVig)->CanFitFlux      = Lc.Ref->IsVariable((*itVig)->Image()->ModifiedModifiedJulianDate());
      (*itVig)->CanFitPos       = (*itVig)->CanFitFlux;
      (*itVig)->CanFitSky       = true;
      (*itVig)->DontConvolve = Lc.Ref->Image()->Name() == (*itVig)->Image()->Name(); 
      // this resizes the vignet and update everything (kernel, psf, residus) 
      (*itVig)->ModifiedResid();
      (*itVig)->AutoResize();
      // check whether there are weights>0 on the position of the star,
      // if not, FitFlux, FitSky =false 
      (*itVig)->CheckWeight();       
    }
}

void SimFit::Resize(const double& ScaleFactor)
{

#ifdef FNAME
  cout << " > SimFit::Resize(" << ScaleFactor << ");" << endl;
#endif
  
  int hrefx = VignetRef->Hx();
  int hrefy = VignetRef->Hy();
  
  //if(true) {
  if(fabs(ScaleFactor-1)>0.01) {
    if ((!fit_flux) && (!fit_pos) && (!fit_gal) && (!fit_sky) || (size()==0)) 
      {
	cerr << " SimFit::Resize(" << ScaleFactor 
	     << ") : Warning: nothing to fit, not resizing. " << endl;
	return;
      }
    
    cout << " SimFit::Resize(" << ScaleFactor << "): old scale = " << scale; 
    scale = max(ScaleFactor, minscale);
    cout << " new scale = " << scale << endl;

  
    // resize vignets
    
    VignetRef->Resize(int(ceil(scale*double(VignetRef->Hx()))),
		      int(ceil(scale*double(VignetRef->Hy()))));
  }
  hrefx = VignetRef->Hx();
  hrefy = VignetRef->Hy();
  
  // anyway, resize and update everything for all vignets
  for (SimFitVignetIterator it = begin(); it != end(); ++it) {
    (*it)->ModifiedResid();
    (*it)->AutoResize();
  }
  

#ifdef DEBUG
  cout <<   "   SimFit::Resize whattofit = " << fit_flux << "," << fit_pos << "," << fit_gal << "," << fit_sky << endl;
#endif

  // recompute matrix indices
  hfx = hfy = nfx = nfy = 0;
  fluxstart = fluxend = 0;
  
  // fluxes
  if (fit_flux)
    {
      for (SimFitVignetCIterator it=begin(); it != end(); ++it)
	if ((*it)->FitFlux) fluxend++;
    }  
  fluxend--; // we want fluxend = nfluxes-1
    
  
  // position
  xind = yind = fluxend;
  if (fit_pos) 
    {
      xind += 1;
      yind += 2;
    }
  
  // galaxy pixels
  galstart = galend = yind;
  if (fit_gal)
    {
      if(fluxend)
	galstart += 1;
      hfx = hrefx;
      hfy = hrefy;
      nfx = 2*hfx+1;
      nfy = 2*hfy+1;
      galend = galstart + nfx*nfy-1;
      refill = true;
    }
  
  // skies
  skystart = skyend = galend;
  if (fit_sky)
    {
      skystart += 1;
      for (SimFitVignetCIterator it=begin(); it != end(); ++it)
	if ((*it)->FitSky) skyend++;
    }
  
  
  nparams = (fluxend-fluxstart+1) + (yind-fluxend) + (galend-yind) + (skyend-galend);

  //#ifdef DEBUG
  cout << " SimFit::Resize(" << ScaleFactor << ") : indices: \n" 
       << "   fluxstart " << fluxstart  << "   fluxend " << fluxend << endl
       << "   xind " << xind << "   yind " << yind << endl
       << "   galstart " << galstart  << " galend " << galend << endl
       << "   skystart " << skystart << " skyend " << skyend << endl
       << "   nparams = " << nparams << endl;
  //#endif

  ndata = 0;
  
  for (SimFitVignetCIterator it = begin(); it != end(); ++it) 
    ndata += (*it)->NValidPixels();
  //ndata += (2*(*it)->Hx()+1) * (2*(*it)->Hy()+1);
  
  Vec.allocate(nparams);
  PMat.allocate(nparams, nparams);
  MatGal.allocate(nfx*nfy,nfx*nfy);
  
#ifdef DEBUG
  cout << "   nparams = " << nparams << endl;
  cout << "   nfx     = " << nfx << endl;
  cout << "   nfy     = " << nfy << endl;
  cout << "   fluxstart,fluxend  = " << fluxstart << "," << fluxend << endl;
  cout << "   xind,yind          = " << xind << "," << yind << endl;
  cout << "   galstart,galend    = " << galstart << "," << galend << endl;
  cout << "   skystart,skyend    = " << skystart << "," << skyend << endl;  
#endif

}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::Matrix filling routines::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IMPORTANT: right now, coordinates are all centered on (0,0), center of 
  the vignet. That makes it easier to read. If you change the origin for 
  the data, weights or psf kernels, you should recode everything. Do the 
  merguez barbecue instead. 
*/


void SimFit::FillMatAndVec()
{

#ifdef FNAME
  cout << " > SimFit::FillMatAndVec() : Initialize " << endl;  
#endif
  
  Vec.Zero();
  PMat.Zero();
  
#ifdef DEBUG
  cout << "   in SimFit::FillMatAndVec() : Compute matrix and vectors " << endl;  
#endif

  if (fit_flux)            fillFluxFlux();
  if (fit_flux && fit_pos) fillFluxPos();
  if (fit_flux && fit_gal) fillFluxGal();
  if (fit_flux && fit_sky) fillFluxSky();

  if (fit_pos)             fillPosPos();
  if (fit_pos && fit_gal)  fillPosGal();
  if (fit_pos && fit_sky)  fillPosSky();

  if (fit_gal)             fillGalGal();
  if (fit_gal && fit_sky)  fillGalSky();

  if (fit_sky)             fillSkySky();

  // no need to symmetrize the matrix

  if(PMat.SizeX()<20 && PMat.SizeY()<20) {
    cout << "===== PMat =====" << endl;
    cout << PMat << endl;
    cout << "===== PMat =====" << endl;
  }
}

int SimFit::galind(const int i, const int j) const
{
  return galstart + (i+hfx)*nfy + (j+hfy);
}

void SimFit::fillFluxFlux()
{
  //*********************************************
  // flux-flux matrix terms and flux vector terms
  //*********************************************

#ifdef FNAME
  cout << " > SimFit::fillFluxFlux()" << endl;
#endif

  // loop over vignets
  int fluxind = 0;
  DPixel *pw, *pres, *ppsf;
  double sumvec, summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;
      
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      sumvec = 0.;
      summat = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
#ifdef DEBUG_TOTO
	  cout << "pw,pres,ppsf " << *pw << "," << *pres << "," << *ppsf << endl;
#endif

	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvec += (*pres) * (*ppsf) * (*pw);
	      summat += (*ppsf) * (*ppsf) * (*pw);
	      ++ppsf; ++pres; ++pw;
	    }
	}

      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;

      Vec(ind) = sumvec;
      PMat(ind,ind) = summat;
      
      ++fluxind;
    }
}

void SimFit::fillFluxPos()
{
  //**********************
  // flux-pos matrix terms
  //**********************

#ifdef FNAME
  cout << " > SimFit::fillFluxPos()" << endl;
#endif

  // loop over vignets
  int fluxind = 0;
  DPixel *pw, *ppsf, *ppdx, *ppdy;
  double summatx, summaty;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;

      if (!vi->FitFlux || !vi->FitPos) continue;
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      summatx = 0.;
      summaty = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
	  ppdx = &(vi->Psf.Dx)(-hx,j);
	  ppdy = &(vi->Psf.Dy)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillFluxPos pw,ppsf,ppdx,ppdy " 
	       << *pw << "," 
	       << *ppsf << ","
	       << *ppdx << ","
	       << *ppdy << endl;
	    
#endif
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summatx += (*ppsf) * (*ppdx) * (*pw);
	      summaty += (*ppsf) * (*ppdy) * (*pw);
	      ++pw; ++ppsf; ++ppdx; ++ppdy;
	    }
	}
      
      summatx *= vi->Star->flux; //JG
      summaty *= vi->Star->flux; //JG
      
      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;


      PMat(xind,ind) = summatx;
      PMat(yind,ind) = summaty;

      fluxind++;
    }  
}

#define KERNIND(HK,H,IND,A,B)\
int A = (HK-IND < H) ? IND-HK : -H;\
int B = (HK+IND < H) ? IND+HK : H;\
A = A<B ? A : B;\
B = A>B ? A : B;\

void SimFit::fillFluxGal()
{ 
  //*********************
  // flux-gal matrix part
  //*********************

#ifdef FNAME
  cout << " > SimFit::fillFluxGal()" << endl;
#endif

  int fluxind = 0;
  DPixel *pw, *ppsf, *pkern;
  double summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;
      
      int hx = vi->Hx();
      int hy = vi->Hy();
  
      int ind = fluxstart+fluxind;

      // the dirac case : sum[ psf(x) * dirac(y) ] = psf(y)	  
      if (vi->DontConvolve)
	{
	  // loop over fitting coordinates 
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      ppsf = &(vi->Psf)   (-hx,j);
	      pw   = &(vi->OptWeight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  PMat(galind(i,j),ind) = (*ppsf) * (*pw);
		  ++ppsf; ++pw;
		}
	    }
	  ++fluxind;
	  continue;
	}
      
      // the normal cases
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      // loop over fitting coordinates
      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      summat = 0.;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  ppsf  = &(vi->Psf)   (ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*ppsf) * (*pkern) * (*pw);
		      ++pkern; ++ppsf; ++pw;
		    }
		}
	      PMat(galind(is,js),ind) = summat;
	    }
	}
      ++fluxind;

    }
}

void SimFit::fillFluxSky()
{
  //***********************
  // flux-sky matrix terms
  //***********************

#ifdef FNAME
  cout << " > SimFit::fillFluxSky()" << endl;
#endif

  int skyind  = 0;
  int fluxind = 0;
  double summat;
  DPixel *pw, *ppsf;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (vi->FitFlux && vi->FitSky) {

      // sum over smallest between psf and data coordinates
      int hx = vi->Hx();
      int hy = vi->Hy();
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summat += (*ppsf) * (*pw);
	      ++ppsf; ++pw;
	    }  
	}
      // now fill in matrix part

      PMat(fluxstart+fluxind,skystart+skyind) = summat;
      }
      if(vi->FitFlux)
	++fluxind; 
      if(vi->FitSky)
	++skyind; 
    }
}

void SimFit::fillPosPos()
{
  //*********************************************
  // pos-pos matrix terms and pos vector terms
  //*********************************************

#ifdef FNAME
  cout << " > SimFit::fillPosPos()" << endl;
#endif

  double sumvecx = 0.;
  double sumvecy = 0.;
  double summatx = 0.;
  double summaty = 0.;
  double summatxy = 0.;
  DPixel *pw, *pres, *ppdx, *ppdy;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitPos) continue;

      // sum over the smallest area between psf derivative and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();

      double flux = vi->Star->flux;
      double flux2 = flux*flux;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  ppdx = &(vi->Psf.Dx)(-hx,j);
	  ppdy = &(vi->Psf.Dy)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillPosPos pw,pres,ppdx,ppdy " 
	       << *pw << "," 
	       << *pres << ","
	       << *ppdx << ","
	       << *ppdy << endl;
#endif
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvecx  += (*pres) * (*ppdx) * (*pw) * flux;
	      sumvecy  += (*pres) * (*ppdy) * (*pw) * flux;
	      summatx  += (*ppdx) * (*ppdx) * (*pw) * flux2;
	      summaty  += (*ppdy) * (*ppdy) * (*pw) * flux2;
	      summatxy += (*ppdx) * (*ppdy) * (*pw) * flux2;
	      ++pw; ++pres; ++ppdx; ++ppdy;
	    }
	}    
    }
  
  Vec(xind) = sumvecx;
  Vec(yind) = sumvecy;
  PMat(xind,xind) = summatx;
  PMat(yind,yind) = summaty;
  PMat(yind,xind) = summatxy;
}

void SimFit::fillPosGal()
{
  //*********************
  // pos-gal matrix part
  //*********************

#ifdef FNAME
  cout << " > SimFit::fillPosGal()" << endl;
#endif

  DPixel *ppdx, *ppdy, *pw, *pkern;
  double summatx, summaty;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitPos) continue;

      int hx = vi->Hx();
      int hy = vi->Hy();
  
      // the dirac case
      if (vi->DontConvolve)
	{
	  // loop over psf deriv:sum[ dpdx(x) * dirac(y) ] = dpdx(y)	  
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      ppdx = &(vi->Psf.Dx)(-hx,j);
	      ppdy = &(vi->Psf.Dy)(-hx,j);
	      pw   = &(vi->OptWeight)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillPosGal pw,ppdx,ppdy " 
	       << *pw << "," 
	       << *ppdx << ","
	       << *ppdy << endl;
#endif
	      for (int i=-hx; i<=hx; ++i)
		{
		  PMat(galind(i,j),xind) += (*ppdx) * (*pw);
		  PMat(galind(i,j),yind) += (*ppdy) * (*pw);
		  ++ppdx; ++ppdy; ++pw;
		}
	    }
	  continue;
	}

      // loop over psf deriv in fitting coordinates
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      summatx = 0.;
	      summaty = 0.;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  ppdx = &(vi->Psf.Dx)(ikstart,jk);
		  ppdy = &(vi->Psf.Dy)(ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summatx += (*ppdx) * (*pkern) * (*pw);
		      summaty += (*ppdy) * (*pkern) * (*pw);
		      ++pkern; ++ppdx; ++ppdy; ++pw;
		    }
		}
	      
	      summatx *= vi->Star->flux; // JG
	      summaty *= vi->Star->flux; // JG
	      
	      PMat(galind(is,js),xind) += summatx;
	      PMat(galind(is,js),yind) += summaty;      

	    }
	} // end of loop on pixels

    } //end of loop on vignets
}

void SimFit::fillPosSky()
{
  //***********************
  // pos-sky matrix terms
  //***********************

#ifdef FNAME
  cout << " > SimFit::fillPosSky()" << endl;
#endif

  int skyind = 0;
  DPixel *pw, *ppdx, *ppdy;
  double summatx, summaty;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (vi->FitSky && vi->FitPos) {

	// sum over the smallest area between psf and vignet
	int hx = vi->Hx();
	int hy = vi->Hy();
	summatx = 0.;
	summaty = 0.;
	
	for (int j=-hy; j<=hy; ++j)
	  {
	    pw   = &(vi->OptWeight)(-hx,j);
	    ppdx = &(vi->Psf.Dx)(-hx,j);
	    ppdy = &(vi->Psf.Dy)(-hx,j);
	    for (int i=-hx; i<=hx; ++i) 
	      {
		summatx += (*ppdx) * (*pw);
		summaty += (*ppdy) * (*pw);
		++ppdx; ++ppdy; ++pw;
	      }
	  }
	
	summatx *= vi->Star->flux; // JG
	summaty *= vi->Star->flux; // JG
	
	// now fill matrix part
	PMat(skystart+skyind,xind) = summatx;
	PMat(skystart+skyind,yind) = summaty;
	
	++skyind;
      }
    }
}

void SimFit::fillGalGal()
{

  //******************************************
  // gal-gal matrix terms and gal vector terms
  //******************************************

#ifdef FNAME
  cout << " > SimFit::fillGalGal()" << endl;
#endif

  double sumvec;
  DPixel *pkern, *pres, *pw;

  // loop over vignets
#ifdef DEBUG_FILLMAT
  cout << "  Loop over vignets ..." << endl;
#endif
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      if((*it)->FitFlux && dont_use_vignets_with_star)
	continue;

      const SimFitVignet *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();

      // the dirac case
      if (vi->DontConvolve)
	{
	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		Vec(galind(i,j)) += (vi->Resid)(i,j) * (vi->OptWeight)(i,j);
	      }
	  continue;
	}


      /* galaxy contribution to vector terms
	 basically a convolution. Can't just simply use the Convolve routine: 
	 gotta choose some rules for the borders, depending on sizes user has chosen */

      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      sumvec = 0.;
	      KERNIND(hky,hy,js,jkstart,jkend);
	      
	      // sum over kernel, stay in fitting coordinates
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  pres  = &(vi->Resid) (ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);

#ifdef DEBUG_FILLMAT
		  //cout << "   fillGalGal pkern,pres,pw " << *pkern << "," << *pres << "," << *pw << endl;
#endif
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      sumvec += (*pkern) * (*pres) * (*pw);
		      ++pkern; ++pres; ++pw;
		    }
		}
	      Vec(galind(is,js)) += sumvec;
	    }
	}
    }

  // now if weights do not change (that is not robustify), 
  // we do not need to refill this part at each iteration
  if (false && !refill) 
    {

      int ngal = galend-galstart+1; 
#ifdef DEBUG
      cout << "     case notrefile " << endl;
#endif
      for (int j=0; j<ngal; ++j) 
	for (int i=j+1; i<ngal; ++i) {
	  PMat(galstart+i,galstart+j) = MatGal(i,j);
	}
      return;
    }
  
  /* galaxy-galaxy matrix terms: longest loop.
     use the simplified relation, which in one dimension can be written as:
     matrix(m,n) = sum_ik ker[ik-(in-im)] * ker(ik) * weight(ik+im). 
     This is also a convolution, but again, borders are messy. */
  
#ifdef DEBUG
  cout << " Computing galaxy-galaxy matrix terms (longest loop) ... " << endl;
  int count = 0;
  int zesize = size();
#endif

  double summat;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {

#ifdef DEBUG
      count++;
      cout << "  " << count << "/" << zesize << endl;
#endif    
       if((*it)->FitFlux && dont_use_vignets_with_star)
	continue;

      const SimFitVignet *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      // the dirac case
      if (vi->DontConvolve)
	{

	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		int ind = galind(i,j);
		PMat(ind,ind) += (vi->OptWeight)(i,j);
	      }
	  continue;
	}

   
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int min_m = -2*(hkx*nfy + hky);
      int npix = nfx*nfy;
      
      int m,n; // index of the galgal matrix
      int im,jm; // index of pixels in the galaxy image, corresponding to matrix index m, 
      int in,jn; // index of pixels in the galaxy image, corresponding to matrix index n
      // m = (im+hfx)*nfy + (jm+hfy)  same as galind
      int imn,jmn;
      int ik,jk;
      
      int ikmin,ikmax,jkmin,jkmax;
      DPixel *pkern1,*pkern2,*pw;
      
      // loop over fitting coordinates
      for ( m=0; m<npix; ++m)
	{
	  im = (m / nfy) - hfx;
	  jm = (m % nfy) - hfy;
	  for ( n=((min_m+m > 0) ? min_m+m : 0); n<=m; ++n) // ???
	    {
	      in = (n / nfy) - hfx;
	      jn = (n % nfy) - hfy;
	      imn = im-in;
	      jmn = jm-jn;
	      summat = 0.;
	      jkmin = max(max(-hky,-hy-jm),-hky-jmn);
	      jkmax = min(min(hky,hy-jm), hky-jmn);
	      
	      for (jk=jkmin; jk<=jkmax; ++jk) {
		ikmin = max(max(-hkx,-hx-im),-hkx-imn);
		ikmax = min(min(hkx,hx-im),hkx-imn);
		pkern1 = &( (vi->Kern)  (ikmin,jk));
		pkern2 = &( (vi->Kern)  (ikmin+imn,jk+jmn));
		pw     = &( (vi->OptWeight) (ikmin+im,jk+jm));
		for (ik=ikmin; ik<=ikmax; ++ik, ++pkern1, ++pkern2, ++pw)  
		  summat += (*pkern1)*(*pkern2)*(*pw);
	      }
	      PMat(galstart+m,galstart+n) += summat; 	
	    }
	}
    }
  
#ifdef DEBUG
  cout << " SimFit::fillGalGal is ending ... " << endl;
  cout << " galstart =" << galstart << endl;
  cout << " galend   =" << galend << endl;
  
#endif
  int ngal = galend-galstart+1;
  for (int j=0; j<ngal; ++j) 
    for (int i=j+1; i<ngal; ++i) { 
      MatGal(i,j) = PMat(galstart+i,galstart+j);
    }
  refill = false;
#ifdef DEBUG
  cout << " SimFit::fillGalGal done " << endl;
#endif

}



void SimFit::fillGalSky()
{
  //**********************
  // gal-sky matrix terms
  //*********************

#ifdef DEBUG
  cout << " SimFit::fillGalSky()" << endl;
#endif

  // loop over vignets
  DPixel *pw, *pkern;
  int skyind = 0;
  double summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if(vi->FitSky) {
      
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      int ind = skystart+skyind;

      // dirac case
      if (vi->DontConvolve)
	{
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      pw = &(vi->OptWeight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  PMat(galind(i,j),ind) = *pw++;
		}
	    }
	  ++skyind;
	  continue;
	}

      // loop over smallest area among fitting size and weight
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      summat = 0.;
	      // sum over kernel
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  pw    = &(vi->OptWeight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*pkern) * (*pw);
		      ++pkern; ++pw;
		    }
		}
	      // now fill matrix elements
	      PMat(galind(is,js),ind) = summat;
	    }
	}
      ++skyind;
      }
    }
}

void SimFit::fillSkySky()
{
  //*********************************************
  // sky-sky matrix terms and flux vector terms
  //*********************************************

#ifdef DEBUG
  cout << " SimFit::fillSkySky()" << endl;
#endif

  // loop over vignets
  int skyind = 0;
  double sumvec, summat;
  DPixel *pres, *pw;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if(vi->FitSky) {

      // sum over vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      sumvec = 0.;
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvec += (*pw) * (*pres);
	      summat += (*pw);
	      ++pres; ++pw;
	    }
	}

      // now fill out the matrix and vector
      int ind = skystart+skyind;
      Vec(ind) = sumvec;
      PMat(ind,ind) = summat;
      ++skyind;
      }
    }
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::    Solving routines   ::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/

double SimFit::computeChi2() const
{
#ifdef FNAME
  cout << " > SimFit::computeChi2()" << endl;
#endif
  double c = 0.;
  int count = 0;
  double ic = 0;
  for (SimFitVignetCIterator it=begin(); it != end(); ++it) {
    if((*it)->FitFlux && dont_use_vignets_with_star)
      continue;
    ic = (*it)->Chi2();
    c += ic;
    count ++;
#ifdef DEBUG
    cout << "   in SimFit::computeChi2() vignet " << count << " chi2 = " << ic << endl;
#endif   
  }
  return c;
}

bool SimFit::Update(double Factor)
{
  if(fit_pos) {
    float maxoffset=1;  // we don't want to get offsets larger than one pixel
    if(fabs(Vec(xind)*Factor)>maxoffset) {
      return Update(0.9*maxoffset/fabs(Vec(xind)));
    }
    if(fabs(Vec(yind)*Factor)>maxoffset) {
      return Update(0.9*maxoffset/fabs(Vec(yind)));
    }
  }

  
#ifdef DEBUG
  cout << " SimFit::Update(" << Factor << "): updating parameters \n";
#endif
  

  // update the galaxy
  if (fit_gal)
    for (int j=-hfy; j<=hfy; ++j) 
      for (int i=-hfx; i<=hfx; ++i)
	VignetRef->Galaxy(i,j) += Vec(galind(i,j))*Factor;

  // update the reference position and psf
  if (fit_pos) 
    {
      //#ifdef DEBUG
      cout << "   in SimFit::Update shifting of dx dy " << Vec(xind)*Factor << " " <<  Vec(yind)*Factor << " (xind,yind)=" << xind << "," << yind <<  endl;
      //#endif
      if (!(VignetRef->ShiftCenter(Point(Vec(xind)*Factor, Vec(yind)*Factor)))) {
	abort();
      }
      VignetRef->UpdatePsfResid();
    }

  // update all vignets
  int fluxind = fluxstart;
  int skyind  = skystart;
  int start = 1;
  
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      SimFitVignet *vi = *it;
      
      if(start && (vi->CanFitFlux)) {
	cout << "      in SimFit::Update DUMP : " << vi->Star->flux << " " << vi->Star->x << " " << vi->Star->y << endl;
	start=0;
      }
      
      // update flux
      if ((fit_flux) && (vi->FitFlux)) {
	vi->Star->flux += Vec(fluxind++)*Factor;
      }
 
      // update pos (not really required if not vignetref)
      if ((fit_pos) && !(vi->ShiftCenter(Point(Vec(xind)*Factor, Vec(yind)*Factor)))) {
	cout << "ShiftCenter failure" << endl;
	abort();
      }
      // update sky
      if (fit_sky && (vi->FitSky)) vi->Star->sky += Vec(skyind++)*Factor;
      
      // update residuals and convolved psf
      vi->ModifiedResid();
      vi->Update();
    }

  return true;
}

double SimFit::oneNRIteration(double OldChi2)
{
#ifdef FNAME
  cout << " > SimFit::oneNRIteration(double OldChi2)" << endl;
#endif
  FillMatAndVec();
  
  if (cholesky_solve(PMat,Vec,"L") != 0) {
    PMat.writeFits("DEBUG_pmat.fits");
    cout << "Vec:" << endl;
    cout << Vec << endl;
    cout << "Quit ..." << endl;
    abort();
  }

  Update(1.);  
  double curChi2 = computeChi2();
  const double minFact = 0.001;
  if (curChi2 > OldChi2) 
    {
#ifdef DEBUG      
      cout << " SimFit::oneNRIteration(" 
	   << OldChi2 << "):  chi2=" << curChi2 
	   << " increased, reducing corrections" << endl;
#endif

      double fact = 1.;
      while ((curChi2 > OldChi2) && (fact > minFact))
	{
	  Update(-fact);
	  fact *= 0.9;
	  Update(fact);
	  curChi2 = computeChi2();
#ifdef DEBUG
	  cout << " SimFit::oneNRIteration(" 
	       << OldChi2 << "):  curChi2=" << curChi2 << " fact=" << fact << endl;
#endif
	}

      // reducing corrections had no effect
      if (curChi2 > OldChi2)
	{
	  Update(-fact);
	  curChi2 = computeChi2();
#ifdef DEBUG
	  cout << " SimFit::oneNRIteration(" 
	       << OldChi2 << "):  chi2=" << curChi2 << " return to beginning \n";
#endif

	  curChi2 = OldChi2;
	}
    }

  return curChi2;
}

bool SimFit::IterateAndSolve(const int MaxIter,  double Eps)
{
#ifdef FNAME
  cout << " > SimFit::IterateAndSolve(const int MaxIter, double Eps): " << MaxIter << " " << Eps << "  Start : \n"
       << (*this) << endl;
#endif
  double oldchi2;
  chi2 = computeChi2();
  int iter = 0;
  double diff = 1000;
  do
    {
      cout << "   in SimFit::IterateAndSolve(" << MaxIter << ") : Iteration # " 
	   << iter << " chi2/dof = " << setprecision(4) << chi2/(ndata - nparams) << endl;
      oldchi2 = chi2;
      chi2 = oneNRIteration(oldchi2);
      diff = fabs(chi2-oldchi2);
      cout << "   diff = " << diff << endl;
    }while ((iter++ < MaxIter) && (diff>Eps));
  
  cout << "   in SimFit::IterateAndSolve() : Iteration # " 
       << iter << " chi2/dof = " << setprecision(4) << chi2/(ndata - nparams) << endl;
 

  cout << "   SimFit::IterateAndSolve(): Finish : \n"
       << (*this) << endl;

  return true;
}

bool SimFit::GetCovariance()
{
#ifdef FNAME
  cout << " > SimFit::GetCovariance()" << endl;
#endif
  
  // someday we should recompute the covariance ala MINOS
  
  int status = cholesky_invert(PMat,"L");
  
  if (status != 0) 
    {      
      cerr << " SimFit::GetCovariance() : Error: inverting failed. Lapack status=: " 
	   << status << endl;
      return false;
    }

  // rescale covariance matrix with estimated global sigma scale factor
  // it corrects for initially under-estimated (ex: correlated) weights if chi2/dof < 1 
  // or for error in our model if chi2/dof > 1

  //double sigscale = chi2 / double(ndata - nparams);
  double sigscale = 1;
  
  int fluxind = 0;
  int skyind  = 0;
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      if ((fit_flux) && (*it)->FitFlux) {
	// (*it)->Star->varflux = sigscale * PMat(fluxind,fluxind++); // DO NOT USE THIS
	(*it)->Star->varflux = sigscale * PMat(fluxind,fluxind);
	fluxind++;
      }
      if (fit_sky) {  
	// (*it)->Star->varsky  = sigscale * PMat(skyind,skyind++); // DO NOT USE THIS
	(*it)->Star->varsky  = sigscale * PMat(skyind,skyind);
	skyind++;
      }
    }

  if (fit_pos)
    {
      VignetRef->Star->varx  = sigscale * PMat(xind,xind);
      VignetRef->Star->vary  = sigscale * PMat(yind,yind);
      VignetRef->Star->covxy = sigscale * PMat(xind,yind);
    }
  
  return true;
}

void SimFit::FitInitialGalaxy() {
#ifdef FNAME
  cout << " > SimFit::FitInitialGalaxy()" << endl;  
#endif
  int nzero=0;
  for (SimFitVignetIterator it=begin(); it != end(); ++it) {
    if(!((*it)->FitFlux))
      nzero ++;
  }
  cout << "using " << nzero << " exposures for reference (in FitInitialGalaxy)" << endl;
  if(nzero==0) {
    abort();
  }
  SetWhatToFit(FitGal);
  dont_use_vignets_with_star = true;
  Resize(1);
  IterateAndSolve(2);
}


void SimFit::DoTheFit(int MaxIter)
{
#ifdef FNAME
  cout << " > SimFit::DoTheFit() : Starting" << endl;  
#endif
  if (!(fit_pos || fit_sky || fit_flux || fit_gal)) 
    {
      cerr << " SimFit::DoTheFit() : Error : nothing to fit" << endl;
      return;
    }

#ifdef DEBUG
  clock_t tstart = clock(); 
#endif 

  Resize(1);
  if (!IterateAndSolve(MaxIter)) return;
  if (!GetCovariance())    return;

#ifdef DEBUG
  clock_t tend = clock();
  cout << " SimFit::DoTheFit() : CPU comsumed : " 
       <<  float(tend-tstart) / float(CLOCKS_PER_SEC) << endl;
#endif
}

void SimFit::write(const string& StarName,const string &DirName, unsigned int whattowrite) 
{
#ifdef FNAME
  cout << " SimFit::write(" << StarName << ")" << endl;
#endif
  // write lc
  if(whattowrite & WriteLightCurve) {
    ofstream lstream(string(DirName+"/lightcurve_"+StarName+".dat").c_str());
    lstream << setiosflags(ios::fixed);
    front()->Star->WriteHeader(lstream);
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      lstream << *((*it)->Star) << endl;
    }
    lstream.close();
  }
  
  // write vignets info
  if(whattowrite & WriteVignetsInfo) {
    ofstream vignetstream(string(DirName+"/vignets_"+StarName+".dat").c_str());
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      SimFitVignet *vi = *it;
      vignetstream << *vi << endl;
    }
    vignetstream.close();
  }
  
  // write galaxy fits
  if(whattowrite & WriteGalaxy)
    VignetRef->Galaxy.writeFits(DirName+"/galaxy_"+StarName+".fits");
  

  // write residuals
  if(whattowrite & WriteResid) {
    for (SimFitVignetIterator it=begin(); it != end() ; ++it){      
      (*it)->ClearResidZeroWeight();
      (*it)->Resid.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_resid.fits");
    }
  }
  
  // write Data
  if(whattowrite & WriteData)
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->Data.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_data.fits");
  
  // write Psf
  if(whattowrite & WritePsf)
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->Psf.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_psf.fits");

  // write weights
  if(whattowrite & WriteWeight)
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->OptWeight.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_weight.fits");
  
  // write kernels
  if(whattowrite & WriteKernel)
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->Kern.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_kern.fits");

  //write matrices
  if(whattowrite & WriteMatrices) {
    PMat.writeFits(DirName+"/pmat_"+StarName+".fits");
    
    int i=0;
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      if((*it)->FitFlux)
	i++;
    Mat vm(1,i);
    i=0;
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      if((*it)->FitFlux) {
	vm(0,i)= (*it)->Star->flux;
	i++;
      }
    vm.writeFits(DirName+"/vec_"+StarName+".fits");
    
    // create matrix of Nights and Images
    fillNightMat();
    NightMat.writeFits(DirName+"/nightmat_"+StarName+".fits");
  }
}

void SimFit::fillNightMat(LightCurve& Lc) {
#ifdef FNAME
  cout << " > SimFit::fillNightMat" << endl;
#endif
  
  
      
  LightCurve::const_iterator itLc = Lc.begin();
  LightCurve::const_iterator endLc = Lc.end();
  double jd; // julian day
  vector<double> nightdates;
  bool isinnight;
  double timediff = 10./24.; // 10 hours
  int nimages = 0;
  for (;itLc != endLc; ++itLc) {
    jd = (*itLc)->Image()->ModifiedModifiedJulianDate();
    if(Lc.Ref->IsVariable(jd)) {
      nimages++;
      isinnight=false;
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  isinnight = true;
	  break;
	}
      }
      if(!isinnight) {
	nightdates.push_back(jd);
      }
    }
  } 
  int nnights = nightdates.size(); // number of nights for this lightcurve
  NightMat.allocate(nnights,nimages); // szie of matrix
  // now we fill this matrix
    //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  itLc = Lc.begin();
  int im = 0;
  for (;itLc != endLc; ++itLc) {
    jd = (*itLc)->Image()->ModifiedModifiedJulianDate();
    if(Lc.Ref->IsVariable(jd)) {
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  NightMat(day,im)=1;
	  break;
	}
      }
      im++;
    }
  }
  
  cout << "========= NightMat ==========" << endl;
  cout << NightMat << endl;
  return;
}


void SimFit::fillNightMat() { // using simfitvignet
#ifdef FNAME
  cout << " > SimFit::fillNightMat" << endl;
#endif
  
  double jd; // julian day
  vector<double> nightdates;
  bool isinnight;
  double timediff = 10./24.; // 10 hours
  int nimages = 0;
  
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    if((*itVig)->FitFlux) {
      jd = (*itVig)->Image()->ModifiedModifiedJulianDate();
      nimages++;
      isinnight=false;
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  isinnight = true;
	  break;
	}
      }
      if(!isinnight) {
	nightdates.push_back(jd);
      }
    }
  } 
  int nnights = nightdates.size(); // number of nights for this lightcurve
  NightMat.allocate(nnights,nimages); // szie of matrix
  // now we fill this matrix
    //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  
  int im = 0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    if((*itVig)->FitFlux) {
      jd = (*itVig)->Image()->ModifiedModifiedJulianDate();
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  NightMat(day,im)=1;
	  break;
	}
      }
      im++;
    }
  }
  //cout << "========= NightMat ==========" << endl;
  //cout << NightMat << endl;
  return;
}


ostream& operator << (ostream& Stream, const CountedRef<SimFitVignet> &Vig)
{
  Stream << *Vig;
  return Stream;
}


ostream& operator << (ostream& Stream, const SimFit &MyFit)
{

  Stream << " SimFit with " << MyFit.size() << " vignettes: " << endl
	 << "   Params : npar" << endl
	 << " -----------------------" << endl;
  if(MyFit.fluxend>0)
    Stream << "     flux : " << MyFit.fluxend - MyFit.fluxstart + 1 << endl;
  else
    Stream << "     flux : " << 0 << endl;

  Stream << " position : " << MyFit.yind - MyFit.fluxend << endl
	 << "   galaxy : " << MyFit.nfx   << "X"  << MyFit.nfy << endl
	 << "      sky : " << MyFit.skyend -  MyFit.galend  << endl
	 << " -----------------------" << endl
	 << "    Total     " << MyFit.nparams << " parameters "<< MyFit.ndata << " data "<< endl
	 << "  scale = " << MyFit.scale << " minscale = " << MyFit.minscale << " refill = "  << MyFit.refill << endl
	 << "  chi2 = "  << MyFit.chi2  << " dof = "      << MyFit.ndata-MyFit.nparams << endl;
#ifdef DEBUG  
//   Stream << "   Reference: " << endl
// 	 << MyFit.VignetRef << endl
// 	 << "   Vignets: " << endl
// 	 << " --------------------------" << endl;

//   copy(MyFit.begin(), MyFit.end(), ostream_iterator<CountedRef<SimFitVignet> >(Stream, "\n"));
#endif
  return Stream;
}


ostream& SimFit::DumpMatrices(ostream& Stream) const {
  Stream << "====  SimFit::DumpMatrices ====" << endl;
  Stream  << endl;
  Stream  << "Vec" << endl;
  Stream  << Vec << endl;
  Stream  << "PMat" << endl;
  Stream  << PMat << endl;
  return Stream;

}

// could not make it work. was designed to replace the ugly test in the fillGalGal loop above
#if 0
DPixel *pkerm, *pkern;
#define MAX3(A,B,C) (A>B ? (A>C ? A : C) : (B>C ? B : C))
#define MIN3(A,B,C) (A<B ? (A<C ? A : C) : (B<C ? B : C))

	  int ikstart = MAX3(-hkx, -hkx-imn, -hsx-im);
	  int jkstart = MAX3(-hky, -hky-jmn, -hsy-jm);
	  int ikend = MIN3(hkx, hkx-imn, hsx-im);
	  int jkend = MIN3(hky, hky-jmn, hsy-jm);
          ikstart = (ikstart < ikend) ? ikstart : ikend;  
          jkstart = (jkstart < jkend) ? jkstart : jkend;  
	  summat = 0.;
	  for (int jk=jkstart; jk<=jkend; ++jk)
	    {
	      pkerm = &(vi->Kern)  (ikstart,jk);
	      pkern = &(vi->Kern)  (ikstart+imn,jk+jmn);
	      pw    = &(vi->OptWeight)(ikstart+im,jk+jm);
	      for (int ik=ikstart; ik<=ikend; ++ik)
		{
		  summat += (*pkerm) * (*pkern) * (*pw); 
		  ++pkerm; ++pkern; ++pw;
		}
	    }
	  Mat[(galstart+m)*nparams+galstart+n] = summat; 	
	}
    }
}
#endif
