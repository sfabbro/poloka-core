#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <ctime>

#include "simfitvignet.h"
#include "simfit.h"

#define FNAME // name of functions ar dumped
#define DEBUG // some output
//#define DEBUG_FILLMAT // lots of debug from matrices filling

static void resize_vec(double* &vec, const int n)
{
#ifdef FNAME
  cout << " > resize_vec(double* &vec, const int n), vec= " << vec << " n= " << n << endl; 
#endif
  if (vec) delete [] vec;
  vec = new double[n];
}

static void resize_mat(double* &mat, const int nx, const int ny)
{
#ifdef FNAME
  cout << " > resize_mat(double* &mat, const int nx, const int ny), mat = " << mat << " nx,ny= " << nx << "," << ny << endl; 
#endif
  if (mat) delete [] mat;
  mat = new double[nx*ny];
}

extern "C" {
  void dposv_(char *, int *, int *, double *, int *, double *, int *, int *);
  void dpotri_(char *, int *, double *, int *, int *);
};

static int cholesky_solve(double *a, double *b, int n)
{
#ifdef FNAME
  cout << " >  cholesky_solve" << endl;
#endif  
  int nhrs = 1, info = 0;
  dposv_("L", &n, &nhrs, a, &n, b, &n, &info);

  if (info != 0) 
    cerr << " cholesky_solve(" << a << "," << b << "," << n
	 << ") : Warning: Cholesky factorization failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

static int cholesky_invert(double *a, int n)
{  
  int info = 0;

  dpotri_("L", &n, a, &n, &info);

  if (info != 0) 
    cerr << " cholesky_invert(" << a << "," << n
	 << ") : Warning: inversion failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  :::::::::::::::::: SimFit stuff   ::::::::::::::::::::::::::
  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

SimFit::SimFit()
{
  refill = true;
  solved = false;
  fit_flux = fit_gal = true; fit_sky = fit_pos = false ;
  use_gal = true;
  fluxstart = galstart = skystart = xind = yind = 0;
  fluxend = galend = skyend = 0;
  scale = 1., minscale = 0.; chi2 = 1e29;
  nfx = nfy = hfx = hfy = nparams = ndata = 0;
  Vec = 0;
  Mat = 0;
  MatGal = 0;
  dont_use_vignets_with_star = false;
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

  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    (*itVig)->FitPos = fit_pos;
    (*itVig)->UseGal = use_gal;
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
  VignetRef->Resize(radius,radius); // now resize, this reloads data, update psf, and makeInitialGalaxy
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
      (*itVig)->FitFlux      = Lc.Ref->IsVariable((*itVig)->Image()->JulianDate());
      (*itVig)->DontConvolve = Lc.Ref->Image()->Name() == (*itVig)->Image()->Name(); 
      // this resizes the vignet and update everything (kernel, psf, residus) 
      (*itVig)->ModifiedResid();
      (*itVig)->AutoResize();
    }
}

void SimFit::Resize(const double& ScaleFactor)
{

#ifdef DEBUG
  cout << " SimFit::Resize(" << ScaleFactor << ");" << endl;
#endif
  
  int hrefx = VignetRef->Hx();
  int hrefy = VignetRef->Hy();
  
  if(fabs(ScaleFactor-1)<0.001) {
#ifdef DEBUG
    cout << " actually do not resize " << endl;
#endif    
  }else{
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
    
    
    hrefx = VignetRef->Hx();
    hrefy = VignetRef->Hy();
    
    //for (SimFitVignetIterator it = begin(); it != end(); ++it)
    //(*it)->AutoResize();
  }
  
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
      
      fluxend--; // we want fluxend = nfluxes-1
    }
  
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
      skyend    = skystart + size()-1;
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
  
  resize_vec(Vec, nparams);
  resize_mat(Mat, nparams, nparams);
  resize_mat(MatGal, nfx*nfy,nfx*nfy);

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

  memset(Vec, 0, nparams*sizeof(double)); // big bug !!
  memset(Mat, 0, nparams*nparams*sizeof(double));

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

#ifdef DEBUG
  cout << "   in SimFit::FillMatAndVec() : Symmetrizing matrix" << endl;
#endif

  for (int i=0; i<nparams; ++i) 
    for (int j=i+1; j<nparams; ++j) 
      Mat[i*nparams+j] = Mat[j*nparams+i];
  //Mat[j*nparams+i] = Mat[i*nparams+j];
#ifdef DEBUG
  //DumpMatrices();
#endif

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
	  pw   = &(vi->Weight)(-hx,j);
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

#ifdef CHECK_MAT_BOUNDS
      cout << " in fillfluxflux ind = " << ind << endl;
      fillVec(ind) = sumvec;
      fillMat(ind*nparams+ind) = summat;
#else
      Vec[ind] = sumvec;
      Mat[ind*nparams+ind] = summat;
#endif      
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

      if (!vi->FitFlux) continue;
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      summatx = 0.;
      summaty = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
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

#ifdef CHECK_MAT_BOUNDS
      cout << " in fillfluxpos ind = " << ind << endl;
      cout << "   SimFit::fillFluxPos fillMat(" << xind*nparams+ind << ")=" <<summatx << endl;
      cout << "   SimFit::fillFluxPos fillMat(" << yind*nparams+ind << ")=" <<summaty << endl;
      fillMat(xind*nparams+ind) = summatx;
      fillMat(yind*nparams+ind) = summaty;
#else
      Mat[xind*nparams+ind] = summatx;
      Mat[yind*nparams+ind] = summaty;
#endif
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
	      pw   = &(vi->Weight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
#ifdef CHECK_MAT_BOUNDS
		  fillMat(galind(i,j)*nparams+ind) = (*ppsf) * (*pw);
#else
		  Mat[galind(i,j)*nparams+ind] = (*ppsf) * (*pw);
#endif
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
		  pw    = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*ppsf) * (*pkern) * (*pw);
		      ++pkern; ++ppsf; ++pw;
		    }
		}
#ifdef CHECK_MAT_BOUNDS
	      fillMat(galind(is,js)*nparams+ind) = summat;
#else
	      Mat[galind(is,js)*nparams+ind] = summat;
#endif
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

  for (SimFitVignetCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;

      // sum over smallest between psf and data coordinates
      int hx = vi->Hx();
      int hy = vi->Hy();
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summat += (*ppsf) * (*pw);
	      ++ppsf; ++pw;
	    }  
	}
      // now fill in matrix part
#ifdef CHECK_MAT_BOUNDS
      fillMat((fluxstart+fluxind)*nparams + skystart+skyind) = summat;
#else
      Mat[(fluxstart+fluxind)*nparams + skystart+skyind] = summat;
#endif
      ++fluxind;
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
      if (!vi->FitFlux) continue;

      // sum over the smallest area between psf derivative and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();

      double flux = vi->Star->flux;
      double flux2 = flux*flux;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
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
// 	      sumvecx  += (*pres) * (*ppdx) * (*pw);
// 	      sumvecy  += (*pres) * (*ppdy) * (*pw);
// 	      summatx  += (*ppdx) * (*ppdx) * (*pw);
// 	      summaty  += (*ppdy) * (*ppdy) * (*pw);
// 	      summatxy += (*ppdx) * (*ppdy) * (*pw);
	      
	      sumvecx  += (*pres) * (*ppdx) * (*pw) * flux;
	      sumvecy  += (*pres) * (*ppdy) * (*pw) * flux;
	      summatx  += (*ppdx) * (*ppdx) * (*pw) * flux2;
	      summaty  += (*ppdy) * (*ppdy) * (*pw) * flux2;
	      summatxy += (*ppdx) * (*ppdy) * (*pw) * flux2;
	      ++pw; ++pres; ++ppdx; ++ppdy;
	    }
	}    
    }


  // now fill in the matrix and vector
#ifdef CHECK_MAT_BOUNDS
  cout << "fillVec("<< xind<<") = " << sumvecx << endl;
  cout << "fillVec("<<yind<<") = " << sumvecy << endl;
  cout << "fillMat("<<xind*nparams + xind<<") = " << summatx << endl;
  cout << "fillMat("<<yind*nparams + yind<<") = " << summaty << endl;
  cout << "fillMat("<<yind*nparams + xind<<") = " << summatxy << endl;
  
  fillVec(xind) = sumvecx;
  fillVec(yind) = sumvecy;
  fillMat(xind*nparams + xind) = summatx;
  fillMat(yind*nparams + yind) = summaty;
  fillMat(yind*nparams + xind) = summatxy;
#else
  Vec[xind] = sumvecx;
  Vec[yind] = sumvecy;
  Mat[xind*nparams + xind] = summatx;
  Mat[yind*nparams + yind] = summaty;
  Mat[yind*nparams + xind] = summatxy;
#endif
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
      if (!vi->FitFlux) continue;

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
	      pw   = &(vi->Weight)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillPosGal pw,ppdx,ppdy " 
	       << *pw << "," 
	       << *ppdx << ","
	       << *ppdy << endl;
#endif
	      for (int i=-hx; i<=hx; ++i)
		{
#ifdef CHECK_MAT_BOUNDS
		  fillMat(galind(i,j)*nparams + xind) += (*ppdx) * (*pw);
		  fillMat(galind(i,j)*nparams + yind) += (*ppdy) * (*pw);
#else
		  Mat[galind(i,j)*nparams + xind] += (*ppdx) * (*pw);
		  Mat[galind(i,j)*nparams + yind] += (*ppdy) * (*pw);
#endif
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
		  pw    = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summatx += (*ppdx) * (*pkern) * (*pw);
		      summaty += (*ppdy) * (*pkern) * (*pw);
		      ++pkern; ++ppdx; ++ppdy; ++pw;
		    }
		}
	      
	      summatx *= vi->Star->flux; // JG
	      summaty *= vi->Star->flux; // JG
	      

#ifdef CHECK_MAT_BOUNDS
	      fillMat(galind(is,js)*nparams+xind) += summatx;
	      fillMat(galind(is,js)*nparams+yind) += summaty;      
#else
	      Mat[galind(is,js)*nparams+xind] += summatx;
	      Mat[galind(is,js)*nparams+yind] += summaty;
#endif
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
  for (SimFitVignetCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;

      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      summatx = 0.;
      summaty = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
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
#ifdef CHECK_MAT_BOUNDS
      cout << "   SimFit::fillPosSky() fillMat(" << (skystart+skyind)*nparams+xind << ")=" <<summatx << endl;
      cout << "   SimFit::fillPosSky() fillMat(" << (skystart+skyind)*nparams+yind << ")=" <<summaty << endl;
      
      fillMat((skystart+skyind)*nparams+xind) = summatx;
      fillMat((skystart+skyind)*nparams+yind) = summaty;
#else
      Mat[(skystart+skyind)*nparams+xind] = summatx;
      Mat[(skystart+skyind)*nparams+yind] = summaty;
#endif
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
#ifdef CHECK_MAT_BOUNDS		
		fillVec(galind(i,j)) += (vi->Resid)(i,j) * (vi->Weight)(i,j);
#else
		Vec[galind(i,j)] += (vi->Resid)(i,j) * (vi->Weight)(i,j);
#endif
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
		  pw    = &(vi->Weight)(ikstart,jk);

#ifdef DEBUG_FILLMAT
		  //cout << "   fillGalGal pkern,pres,pw " << *pkern << "," << *pres << "," << *pw << endl;
#endif

		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      sumvec += (*pkern) * (*pres) * (*pw);
		      ++pkern; ++pres; ++pw;
		    }
		}
#ifdef DEBUG_FILLMAT
	      //cout << " Vec(" << galind(is,js) << ") += " << sumvec << endl;
#endif
#ifdef CHECK_MAT_BOUNDS
	      fillVec(galind(is,js)) += sumvec;
#else
	      Vec[galind(is,js)] += sumvec;
#endif

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
#ifdef CHECK_MAT_BOUNDS
	  fillMat((galstart+i)*nparams+j+galstart) = fillMatGal(i*ngal+j);
#else
	  Mat[(galstart+i)*nparams+j+galstart] = MatGal[i*ngal+j];
#endif
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
#ifdef CHECK_MAT_BOUNDS		
		fillMat(ind*nparams+ind) += (vi->Weight)(i,j);
#else
		Mat[ind*nparams+ind] += (vi->Weight)(i,j);
#endif
	      }
	  continue;
	}

      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int min_m = -2*(hkx*nfy + hky);
      int npix = nfx*nfy;

      // loop over fitting coordinates
      for (int m=0; m<npix; ++m)
	{
	  int im = (m / nfy) - hfx;
	  int jm = (m % nfy) - hfy;
	  for (int n=((min_m+m > 0) ? min_m+m : 0); n<=m; ++n)
	    {
	      int in = (n / nfy) - hfx;
	      int jn = (n % nfy) - hfy;
	      int imn = im-in;
	      int jmn = jm-jn;
	      summat = 0.;
	      for (int jk=-hky; jk<=hky; ++jk)
		for (int ik=-hkx; ik<=hkx; ++ik)
		  {
		    // the ugly test (see for comments below)
		    if ((ik+im >= -hx) && (ik+im <= hx) && 
			(jk+jm >= -hy) && (jk+jm <= hy) && 			
			(ik+imn >= -hkx) && (ik+imn <= hkx) &&
			(jk+jmn >= -hky) && (jk+jmn <= hky))
		      {
			summat += (vi->Kern)  (ik,jk)
			        * (vi->Kern)  (ik+imn,jk+jmn)
			        * (vi->Weight)(ik+im,jk+jm);
		      }	    
		  }
#ifdef CHECK_MAT_BOUNDS
	      fillMat((galstart+m)*nparams+galstart+n) += summat; 	
#else
	      Mat[(galstart+m)*nparams+galstart+n] += summat; 
#endif
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
#ifdef CHECK_MAT_BOUNDS      
      fillMatGal(i*ngal+j) = fillMat((galstart+i)*nparams+j+galstart);
#else
      MatGal[i*ngal+j] = Mat[(galstart+i)*nparams+j+galstart];
#endif
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

  for (SimFitVignetCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const SimFitVignet *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      int ind = skystart+skyind;

      // dirac case
      if (vi->DontConvolve)
	{
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      pw = &(vi->Weight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
#ifdef CHECK_MAT_BOUNDS		  
		  fillMat(galind(i,j)*nparams+ind) = *pw++;
#else
		  Mat[galind(i,j)*nparams+ind] = *pw++;
#endif
		}
	    }
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
		  pw    = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*pkern) * (*pw);
		      ++pkern; ++pw;
		    }
		}
	      // now fill matrix elements
#ifdef CHECK_MAT_BOUNDS
	      fillMat(galind(is,js)*nparams+ind) = summat;
#else
	      Mat[galind(is,js)*nparams+ind] = summat;
#endif
	    }
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
  for (SimFitVignetCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const SimFitVignet *vi = *it;
      
      // sum over vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      sumvec = 0.;
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
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
#ifdef CHECK_MAT_BOUNDS
      fillVec(ind) = sumvec;
      fillMat(ind*nparams+ind) = summat;
#else
      Vec[ind] = sumvec;
      Mat[ind*nparams+ind] = summat;
#endif
      
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
  // update only if is fitted
  if (!solved) return false;

#ifdef DEBUG
  cout << " SimFit::Update(" << Factor << "): updating parameters \n";
#endif

  // update the galaxy
  if (fit_gal)
    for (int j=-hfy; j<=hfy; ++j) 
      for (int i=-hfx; i<=hfx; ++i)
	VignetRef->Galaxy(i,j) += Vec[galind(i,j)]*Factor;

  // update the reference position and psf
  if (fit_pos) 
    {
#ifdef DEBUG
      cout << "   in SimFit::Update shifting of x  y  " << VignetRef->Star->x << " " << VignetRef->Star-> y << endl
	   << " dx dy " << Vec[xind]*Factor << " " <<  Vec[yind]*Factor << " (xind,yind)=" << xind << "," << yind <<  endl;
#endif
      if (!(VignetRef->ShiftCenter(Point(Vec[xind]*Factor, Vec[yind]*Factor)))) {
	abort();
      }
      VignetRef->UpdatePsfResid();
    }

  // update all vignets
  int fluxind = fluxstart;
  int skyind  = skystart;
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      SimFitVignet *vi = *it;

      // update flux
      if ((fit_flux) && (vi->FitFlux)) vi->Star->flux += Vec[fluxind++]*Factor;

      // update pos (not really required if not vignetref)
      if ((fit_pos) && !(vi->ShiftCenter(Point(Vec[xind]*Factor, Vec[yind]*Factor)))) {
	cout << "ShiftCenter failure" << endl;
	return false;
      }
      // update sky
      if (fit_sky) vi->Star->sky += Vec[skyind++]*Factor;
      
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

  ofstream mat("mat.dat");
  DumpMatrices(mat);
  mat.close();

  solved = (cholesky_solve(Mat,Vec,nparams) == 0);
  //#include <vutils.h>
  //solved = MatSolve(Mat,nparams,Vec);
  if (!solved) return 1e29;
  //DumpMatrices();

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

  do
    {
      cout << "   in SimFit::IterateAndSolve() : Iteration # " 
	   << iter << " chi2/dof = " << setprecision(4) << chi2/(ndata - nparams) << endl;
      oldchi2 = chi2;
      chi2 = oneNRIteration(oldchi2);
#ifdef DEBUG
      cout << "   in SimFit::IterateAndSolve fabs(chi2-oldchi2)/((chi2+oldchi2)/2} = " << (fabs(chi2-oldchi2)/((chi2+oldchi2)/2)) << endl;
#endif 
    }while ((iter++ < MaxIter) && (fabs(chi2-oldchi2) > ((chi2+oldchi2)/2.*Eps)));
  
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
  if (!solved) 
    {
      cerr << " SimFit::GetCovariance() : system not solved yet" << endl;
      return false;
    }

  cout << " SimFit::GetCovariance() : starting" << endl;

  // someday we should recompute the covariance ala MINOS

  int status = cholesky_invert(Mat,nparams);

  if (status != 0) 
    {      
      cerr << " SimFit::GetCovariance() : Error: inverting failed. Lapack status=: " 
	   << status << endl;
      return false;
    }

  // rescale covariance matrix with estimated global sigma scale factor
  // it corrects for initially under-estimated (ex: correlated) weights if chi2/dof < 1
  // or for error in our model if chi2/dof > 1

  double sigscale = chi2 / double(ndata - nparams);
  int fluxind = 0;
  int skyind  = 0;
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      if ((fit_flux) && (*it)->FitFlux) 
	            (*it)->Star->varflux = sigscale * Mat[fluxind*nparams+fluxind++];
      if (fit_sky)  (*it)->Star->varsky  = sigscale * Mat[skyind*nparams+skyind++];
    }

  if (fit_pos)
    {
      VignetRef->Star->varx  = sigscale * Mat[xind*nparams+xind];
      VignetRef->Star->vary  = sigscale * Mat[yind*nparams+yind];
      VignetRef->Star->covxy = sigscale * Mat[xind*nparams+yind];
    }
  
  return true;
}

void SimFit::FitInitialGalaxy() {
#ifdef FNAME
  cout << " > SimFit::FitInitialGalaxy()" << endl;  
#endif
  SetWhatToFit(FitGal);
  dont_use_vignets_with_star = true;
  DoTheFit();
}


void SimFit::DoTheFit()
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
  
  // If position is to be fitted, iterate more (non linear), 
  // fit data to its mininimum vignet size, and do not fit sky

  const double oldscale = scale;

  if (fit_pos)
    {
#ifdef DEBUG
      cout << " SimFit::DoTheFit(): fit position and resize vignets at min scale" << endl;
#endif
      bool oldfit_sky = fit_sky;
      bool oldfit_gal = fit_gal;
      fit_sky = false;
      fit_gal = false;	
      if (fit_gal) 
	Resize(minscale);
      else
	Resize(1);
      if (!IterateAndSolve(30)) return;
      if (!GetCovariance())     return;
      fit_sky = oldfit_sky;
      fit_gal = oldfit_gal;
    }

  //return;

  // Freeze position, resize the vignets and fit linearly the rest
  bool oldfit_pos = fit_pos;
  fit_pos = false;
  // one has to tell all vignets we don't fit the pos anymore
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig)
    (*itVig)->FitPos = fit_pos;
  
  // only do those iteration if you are not at maximum scale, we want to fit the sky
  // or we did not fit anything yet.
  if (minscale != oldscale || fit_sky || oldfit_pos)
    {
#ifdef DEBUG
      cout << " SimFit::DoTheFit() : Freeze position and resize vignets at " << endl;
#endif
      Resize(oldscale);
      if (!IterateAndSolve(7)) return;
      if (!GetCovariance())    return;
    }
  
  fit_pos = oldfit_pos;


#ifdef DEBUG
  clock_t tend = clock();
  cout << " SimFit::DoTheFit() : CPU comsumed : " 
       <<  float(tend-tstart) / float(CLOCKS_PER_SEC) << endl;
#endif

}

void SimFit::write(const string& StarName,const string &DirName) 
{
  cout << " SimFit::write(" << StarName << ")" << endl;
  // dump light curve and vignets
  ofstream lstream(string(DirName+"/lightcurve_"+StarName+".dat").c_str());
  ofstream vignetstream(string(DirName+"/vignets_"+StarName+".dat").c_str());
  lstream << setiosflags(ios::fixed);
  front()->Star->WriteHeader(lstream);
  const string name = DirName+"/simfit_vignets_"+StarName+".fits";
  VignetRef->Galaxy.writeFits(DirName+"/galaxy_"+StarName+".fits");
  for (SimFitVignetIterator it=begin(); it != end() ; ++it)
    {      
      SimFitVignet *vi = *it;
      vi->ClearResidZeroWeight();
      lstream << *vi->Star << endl;
      vignetstream << *vi << endl;
      vi->Kern.writeFits(DirName+"/"+vi->Image()->Name()+"_"+StarName+"_kern.fits");
      vi->Resid.writeFits(DirName+"/"+vi->Image()->Name()+"_"+StarName+"_resid.fits");
      vi->Data.writeFits(DirName+"/"+vi->Image()->Name()+"_"+StarName+"_data.fits");
      vi->Psf.writeFits(DirName+"/"+vi->Image()->Name()+"_"+StarName+"_psf.fits");
      vi->Weight.writeFits(DirName+"/"+vi->Image()->Name()+"_"+StarName+"_weight.fits");
    }
  lstream.close();
  vignetstream.close();

  // dump covariance matrix
  if (Mat != 0 && Mat != 0) 
    {
      ofstream cstream(string(DirName+"/cov_"+StarName+".dat").c_str());
      cstream << setiosflags(ios::fixed);
      for(int j=fluxstart; j<=fluxend; ++j)
	{
	  for(int i=fluxstart; i<=fluxend; ++i)
	    cstream << Mat[i*nparams+j] << " ";
	  cstream << endl;
	}
    }
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
  Stream << "   Reference: " << endl
	 << MyFit.VignetRef << endl
	 << "   Vignets: " << endl
	 << " --------------------------" << endl;

  copy(MyFit.begin(), MyFit.end(), ostream_iterator<CountedRef<SimFitVignet> >(Stream, "\n"));
#endif
  return Stream;
}


ostream& SimFit::DumpMatrices(ostream& Stream) const {
  Stream << "====  SimFit::DumpMatrices ====" << endl;
  Stream  << endl;
  ios::fmtflags oldflags = Stream.flags();
  //Stream << setiosflags(ios::fixed);
  for(int i=0;i<nparams;i++) {
    //if(i%10==0)
    Stream << "Vec[" << i << "]\t= " << float(Vec[i]) << endl;
  }
  for(int j=0;j<nparams;j++) {
    Stream << "Mat[i*nparams+" << j << "]\t= ";
    for(int i=0;i<nparams;i++) {
      Stream << Mat[i*nparams+j] << " ";
    }
    Stream << endl;
  }
//   for(int j=0;j<nfy;j++) {
//     Stream << "MatGal[i+" << j << "*nfx]\t= ";
//     for(int i=0;i<nfx;i++) {
//       Stream << MatGal[i+j*nfx] << " ";
//     }
//     Stream << endl;
//   } 
  Stream.flags(oldflags);
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
	      pw    = &(vi->Weight)(ikstart+im,jk+jm);
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
