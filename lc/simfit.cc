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
}

void SimFit::SetWhatToFit(unsigned int ToFit)
{
  fit_flux = (ToFit & FitFlux);
  fit_pos  = (ToFit & FitPos);
  fit_gal  = (ToFit & FitGal);
  fit_sky  = (ToFit & FitSky);
  use_gal |= fit_gal; // use the galaxy if we fit it
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

void SimFit::Load(LightCurve& Lc)
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
      (*itVig)->SetStar(*itLc); // just set the star
      (*itVig)->BuildKernel(); // rebuild kernel 
      if((*itVig)->Kern.HSizeX()>worst_kernel)
	worst_kernel = (*itVig)->Kern.HSizeX();
      if((*itVig)->Kern.HSizeY()>worst_kernel)
	worst_kernel = (*itVig)->Kern.HSizeY();
    }
  
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
  VignetRef->SetStar(Lc.Ref); // just set the star
  VignetRef->Resize(radius,radius); // now resize, this reloads data, update psf, and makeInitialGalaxy
#ifdef DEBUG
  SimFitRefVignet *toto = VignetRef;
  cout << " in SimFit::Load,  VignetRef     = " << toto << endl;
  cout << " in SimFit::Load,  VignetRef->Hx() = " <<  VignetRef->Hx() << endl;
  cout << " in SimFit::Load,  VignetRef->Hy() = " <<  VignetRef->Hy() << endl;
#endif 

  // now resize all vignets and initialize residus
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig, ++itLc)
    {           
      // we will fit the flux according to Lc.Ref and date
      (*itVig)->FitFlux      = Lc.Ref->IsVariable((*itVig)->Image()->JulianDate()); 
      (*itVig)->DontConvolve = Lc.Ref->Image()->Name() == (*itVig)->Image()->Name();
      (*itVig)->FitPos = fit_pos;
      (*itVig)->UseGal = use_gal;
      
      // this resizes the vignet and update everything (kernel, psf, residus) 
      (*itVig)->AutoResize();
    }
}

void SimFit::Resize(double ScaleFactor)
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
    
    scale = max(ScaleFactor, minscale);
    cout << " SimFit::Resize(" << ScaleFactor 
	 << "): with chosen scale factor = " << scale << endl;
    
    // resize vignets
    
    VignetRef->Resize(int(ceil(scale*double(VignetRef->Hx()))),
		     int(ceil(scale*double(VignetRef->Hy()))));
    
    
    hrefx = VignetRef->Hx();
    hrefy = VignetRef->Hy();
    
    // resize and update everything for all vignets
    for (SimFitVignetIterator it = begin(); it != end(); ++it)
      (*it)->AutoResize();
  }
  
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
  ndata = 0;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it) 
    ndata += (2*(*it)->Hx()+1) * (2*(*it)->Hy()+1);

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

  memset(Vec, 0, nparams*sizeof(double));
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
    for (int j=i+1; j<nparams; ++j) Mat[i*nparams+j] = Mat[j*nparams+i];

#ifdef DEBUG
  DumpMatrices();
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
#ifdef DEBUG
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
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summatx += (*ppsf) * (*ppdx) * (*pw);
	      summaty += (*ppsf) * (*ppdy) * (*pw);
	      ++pw; ++ppsf; ++ppdx; ++ppdy;
	    }
	}

      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;

#ifdef CHECK_MAT_BOUNDS
      fillMat(xind*nparams+ind) = summatx;
      fillMat(yind*nparams+ind) = summaty;
#else
      Mat[xind*nparams+ind] = summatx;
      Mat[yind*nparams+ind] = summaty;
#endif
      ++fluxind;
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

      // sum over the smallest area between psf derivative and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->Weight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  ppdx = &(vi->Psf.Dx)(-hx,j);
	  ppdy = &(vi->Psf.Dy)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvecx  += (*pres) * (*ppdx) * (*pw);
	      sumvecy  += (*pres) * (*ppdy) * (*pw);
	      summatx  += (*ppdx) * (*ppdx) * (*pw);
	      summaty  += (*ppdy) * (*ppdy) * (*pw);
	      summatxy += (*ppdx) * (*ppdy) * (*pw);
	      ++pw; ++pres; ++ppdx; ++ppdy;
	    }
	}    
    }

  // now fill in the matrix and vector
#ifdef CHECK_MAT_BOUNDS
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
      // now fill matrix part
#ifdef CHECK_MAT_BOUNDS
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
		  cout << " pkern pres pw " << *pkern << " " << *pres << " " << *pw << endl;
#endif

		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      sumvec += (*pkern) * (*pres) * (*pw);
		      ++pkern; ++pres; ++pw;
		    }
		}
#ifdef DEBUG_FILLMAT
	      cout << " Vec(" << galind(is,js) << ") += " << sumvec << endl;
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
  if (!refill) 
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
      if (!(VignetRef->ShiftCenter(Point(Vec[xind]*Factor, Vec[yind]*Factor)))) return false;
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
      cout << "HELLO !!!!" << endl;
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
  solved = (cholesky_solve(Mat,Vec,nparams) == 0);
  if (!solved) return 1e29;

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

  // rescale covariance matrix with estimated global scale factor
  // it corrects for initially under-estimated (ex: correlated) weights if chi2/dof < 1
  // or for error in our model if chi2/dof > 1

  const double scale = chi2 / double(ndata - nparams);
  int fluxind = 0;
  int skyind  = 0;
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      if ((fit_flux) && (*it)->FitFlux) (*it)->Star->varflux = scale * Mat[fluxind*nparams+fluxind++];
      if (fit_sky)  (*it)->Star->varsky  = scale * Mat[skyind*nparams+skyind++];
    }

  if (fit_pos)
    {
      VignetRef->Star->varx  = scale * Mat[xind*nparams+xind];
      VignetRef->Star->vary  = scale * Mat[yind*nparams+yind];
      VignetRef->Star->covxy = scale * Mat[xind*nparams+yind];
    }
  
  return true;
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
      fit_sky = false;
      if (fit_gal) Resize(minscale);
      if (!IterateAndSolve(30)) return;
      if (!GetCovariance())     return;
      fit_sky = oldfit_sky;
    }

  // Freeze position, resize the vignets and fit linearly the rest

  bool oldfit_pos = fit_pos;
  fit_pos = false;
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

void SimFit::write(const string& StarName) const
{
  // dump light curve and vignets
  ofstream lstream(string("lightcurve_"+StarName+".dat").c_str());
  lstream << setiosflags(ios::fixed);
  front()->Star->WriteHeader(lstream);
  const string name = "simfit_vignets_"+StarName+".fits";
  
  for (SimFitVignetCIterator it=begin(); it != end() ; ++it)
    {
      const SimFitVignet *vi = *it;
      lstream << *vi->Star << endl;
      //vi->writeInImage(name);
    }
  lstream.close();

  // dump covariance matrix
  if (Mat != 0 && Mat != 0) 
    {
      ofstream cstream(string("cov_"+StarName+".dat").c_str());
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


ostream& operator << (ostream& Stream, const SimFit &SimFit)
{

  Stream << " SimFit with " << SimFit.size() << " vignettes: " << endl
	 << "   Params : npar" << endl
	 << " -----------------------" << endl
	 << "     flux : " << SimFit.fluxend - SimFit.fluxstart + 1 << endl
	 << " position : " << SimFit.yind - SimFit.fluxend << endl
	 << "   galaxy : " << SimFit.nfx   << "X"  << SimFit.nfy << endl
	 << "      sky : " << SimFit.skyend -  SimFit.galend  << endl
	 << " -----------------------" << endl
	 << "    Total     " << SimFit.nparams << " parameters "<< SimFit.ndata << " data "<< endl
	 << "  scale = " << SimFit.scale << " minscale = " << SimFit.minscale << " refill = "  << SimFit.refill << endl
	 << "  chi2 = "  << SimFit.chi2  << " dof = "      << SimFit.ndata-SimFit.nparams << endl;
  
  Stream << "   Reference: " << endl
	 << SimFit.VignetRef << endl
	 << "   Vignets: " << endl
	 << " --------------------------" << endl;

  copy(SimFit.begin(), SimFit.end(), ostream_iterator<CountedRef<SimFitVignet> >(Stream, "\n"));

  return Stream;
}


void SimFit::DumpMatrices() {
  cout << "====  SimFit::DumpMatrices ====" << endl;
  cout  << endl;
  for(int i=0;i<nparams;i++) {
    //if(i%10==0)
    cout << "Vec[" << i << "]\t= " << float(Vec[i]) << endl;
  }
  return;
  for(int j=0;j<nparams;j++) {
    cout << "Mat[i+" << j << "*nparams]\t= ";
    for(int i=0;i<nparams;i++) {
      cout << Mat[i+j*nparams] << " ";
    }
    cout << endl;
  }
  for(int j=0;j<nfy;j++) {
    cout << "MatGal[i+" << j << "*nfx]\t= ";
    for(int i=0;i<nfx;i++) {
      cout << MatGal[i+j*nfx] << " ";
    }
    cout << endl;
  } 
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
