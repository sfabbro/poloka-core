#include <ctime>
#include <fstream>
#include <iomanip>

#include "vignetfit.h"
#include "simultaneousfit.h"
#include "vutils.h"
#include "lapackutils.h"
#include "photstar.h"
#include "dimage.h"
#include <blas++.h>

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  :::::::::::::::::: SimultaneousFit stuff   ::::::::::::::::::::::::::
  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#define KERNIND(HK,H,IND,A,B)\
int A = (HK-IND < H) ? IND-HK : -H;\
int B = (HK+IND < H) ? IND+HK : H;\
A = A<B ? A : B;\
B = A>B ? A : B;\

SimultaneousFit::SimultaneousFit() 
  : refill(true),fit_flux(true),fit_gal(true),fit_sky(false),fit_pos(false),
    nfx(0), nfy(0), hfx(0), hfy(0),
    fluxstart(0),galstart(0),skystart(0), xind(0), yind(0), 
    fluxend(0),galend(0),skyend(0),scale(1), minscale(0), 
    nparams(0),ndata(0),iter(0),chi2(0),printlevel(false),
    Minim(NewtonRaphson),VignetRef(NULL), galaxy(NULL)
{
}

SimultaneousFit::~SimultaneousFit()
{
  if (galaxy) delete galaxy;
}

void SimultaneousFit::SetWhatToFit(const int ToFit)
{
  fit_flux = (ToFit & FitFlux);
  fit_pos = (ToFit & FitPos);
  fit_gal = (ToFit & FitGal);
  fit_sky = (ToFit & FitSky);
}

void SimultaneousFit::Resize(const double &ScaleFactor)
{

  if ((!fit_flux) && (!fit_pos) && (!fit_gal) && (!fit_sky) || (size()==0)) 
    {
      cerr << " SimultaneousFit::Resize detected nothing to fit " << endl;
      return;
    }
  // resize vignets
  scale = (ScaleFactor > minscale) ? ScaleFactor : minscale;
  cout << " Resizing vignets with scale factor = " << scale << endl;

  VignetRef->Resize(scale);
  int hrefx = VignetRef->Hx();
  int hrefy = VignetRef->Hy();
  for (VignetFitIterator it = begin(); it != end(); ++it) 
    {
      VignetFit *vi = *it;
      int hnewx = hrefx - vi->Kern.HSizeX();
      int hnewy = hrefy - vi->Kern.HSizeY();
      vi->Resize(hnewx, hnewy);
    }

  // redo indices
  // fluxes
  fluxstart = 0;
  fluxend  = 0;
  int nflux = 0;
  if (fit_flux)
    {
      for (VignetFitCIterator it=begin(); it != end(); ++it)
	if ((*it)->IsStarHere) ++nflux;
      fluxend = nflux-1;
    }

  // position
  xind = fluxend;
  yind = fluxend;
  int npos = 0;
  if (fit_pos) 
    {
      xind += 1;
      yind += 2;
      npos = 2;
    }
  
  // galaxy
  galstart = yind;
  galend = yind;
  nfx = 0;
  nfy = 0;
  int holdx = hfx > hrefx ? hrefx : hfx;
  int holdy = hfy > hrefy ? hrefy : hfy;
  hfx = 0;
  hfy = 0;
  if (galaxy) 
    {
      Vignet temp(*galaxy);      
      delete galaxy;
      galaxy = new Vignet(VignetRef->ic+VignetRef->dxc, VignetRef->jc+VignetRef->dyc, hrefx, hrefy);
      for (int j=-holdy; j<=holdy; ++j)
	for (int i=-holdx; i<=holdx; ++i)
	  {
	    (*galaxy)(i,j) = temp(i,j);
	  }
    }
  else galaxy = new Vignet(VignetRef->ic+VignetRef->dxc, VignetRef->jc+VignetRef->dyc, hrefx, hrefy);

  if (fit_gal)
    {
      galstart += 1;
      hfx = galaxy->HSizeX();
      hfy = galaxy->HSizeY();
      nfx = 2*hfx+1;
      nfy = 2*hfy+1;
      galend = galstart + nfx*nfy-1;
      refill = true;
    }
  
  // sky
  skystart = galend;
  skyend = galend;
  int nsky = 0;
  if (fit_sky)
    {
      skystart += 1;
      nsky = size();
      skyend = skystart+nsky-1;
    }

  nparams = nflux + npos + nfx*nfy + nsky;
  ndata = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it) 
    ndata += (2*(*it)->Hx()+1) * (2*(*it)->Hy()+1);

  cout << "   To fit : ?  npar" << endl
       << " -----------------------" << endl
       << "    flux  : " << fit_flux << "  " << nflux << endl
       << " position : " << fit_pos << "  " << npos << endl
       << "   galaxy : " << fit_gal << "  " << nfx << "X" <<nfy << endl
       << "     sky  : " << fit_sky << "  " << nsky << endl
       << " -----------------------" << endl
       << "    Total     " << nparams << " parameters "<< ndata << " data "<< endl;

  Vec.resize(nparams,1);
  Mat.resize(nparams,nparams);
  MatGal.resize(nfx*nfy,nfx*nfy);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::Matrix filling routines::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IMPORTANT: right now, coordinates are all centered on (0,0), center of 
  the vignet. That makes it easier to read. If you change the origin for 
  the data, weights or psf kernels, you should recode everything. Do the 
  merguez barbecue instead. 
*/


void SimultaneousFit::FillMatAndVec()
{
  if (printlevel)  cout << " Initialize " << endl;  
  Vec = 0;
  Mat = 0;

  if (printlevel) cout << " Compute matrix and vectors " << endl;  
  if (fit_flux) fillFluxFlux();
  if (fit_flux && fit_pos) fillFluxPos();
  if (fit_flux && fit_gal) fillFluxGal();
  if (fit_flux && fit_sky) fillFluxSky();

  if (fit_pos)  fillPosPos();
  if (fit_pos && fit_gal)  fillPosGal();
  if (fit_pos && fit_sky)  fillPosSky();

  if (fit_gal) fillGalGal();
  if (fit_gal && fit_sky)  fillGalSky();

  if (fit_sky) fillSkySky();

  if (printlevel) cout << " Symmetrizing matrix" << endl;
  for (int i=0; i<nparams; ++i) 
    for (int j=i+1; j<nparams; ++j) Mat(i,j) = Mat(j,i);

  if (printlevel) cout << " Finished computation of matrix and vector" << endl;
}

int SimultaneousFit::galind(const int i, const int j) const
{
  return galstart + (i+hfx)*nfy + (j+hfy);
}

#define OPTIMIZED

void SimultaneousFit::fillFluxFlux()
{
  //*********************************************
  // flux-flux matrix terms and flux vector terms
  //*********************************************
  if (printlevel) cout << " Filling flux-flux..." << endl;

  // loop over vignets
  int fluxind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      double sumvec = 0;
      double summat = 0;

#ifndef OPTIMIZED

      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  {
	    double dum = vi->Psf(i,j) * vi->Weight(i,j);
	    sumvec += vi->Resid(i,j) * dum;
	    summat += vi->Psf(i,j) * dum;
	  }

#else

      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *presid = &(vi->Resid)(-hx,j);
	  DPixel *ppsf = &(vi->Psf)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      double dum = (*ppsf) * (*pw);
	      sumvec += (*presid) * dum;
	      summat += (*ppsf) * dum;
	      ++ppsf; ++presid; ++pw;
	    }
	}

#endif

      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;
      Vec(ind) = sumvec;
      Mat(ind,ind) = summat;
      ++fluxind;
    }
}

void SimultaneousFit::fillFluxPos()
{
  //**********************
  // flux-pos matrix terms
  //**********************
  if (printlevel) cout << " Filling flux-pos..." << endl;

  // loop over vignets
  int fluxind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      double summatx = 0;
      double summaty = 0;

#ifndef OPTIMIZED

      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  {
	    double dum = vi->Psf(i,j) * vi->Weight(i,j);
	    summatx += vi->DpDx(i,j) * dum;
	    summaty += vi->DpDy(i,j) * dum;
	  }

#else

      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *ppsf = &(vi->Psf)(-hx,j);
	  DPixel *pdpdx = &(vi->DpDx)(-hx,j);
	  DPixel *pdpdy = &(vi->DpDy)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summatx += (*ppsf) * (*pdpdx) * (*pw);
	      summaty += (*ppsf) * (*pdpdy) * (*pw);
	      ++pw; ++ppsf; ++pdpdx; ++pdpdy;
	    }
	}

#endif

      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;
      Mat(xind,ind) = summatx;
      Mat(yind,ind) = summaty;
      ++fluxind;
    }
  
}

void SimultaneousFit::fillFluxGal()
{ 
  //*********************
  // flux-gal matrix part
  //*********************

  if (printlevel) cout << " Filling flux-gal..." << endl;

  // loop over vignets
  int fluxind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      
      int hx = vi->Hx();
      int hy = vi->Hy();
  
      int ind = fluxstart+fluxind;

      // the dirac case : sum[ psf(x) * dirac(y) ] = psf(y)	  
      if (vi->IsRefResolution)
	{
	  // loop over fitting coordinates 

#ifndef OPTIMIZED

	  for (int j=-hy;  j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      Mat(galind(im,j),ind) = vi->Psf(i,j) * vi->Weight(i,j);

#else

	  for (int j=-hy;  j<=hy; ++j)
	    {
	      DPixel *ppsf = &(vi->Psf)(-hx,j);
	      DPixel *pw = &(vi->Weight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  Mat(galind(i,j),ind) = (*ppsf) * (*pw);
		  ++ppsf; ++pw;
		}
	    }
#endif
	  ++fluxind;
	  continue;
	}

      // loop over fitting coordinates
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

#ifndef OPTIMIZED

      for (int is=-hsx;  is<=hsx; ++is)
	for (int js=-hsy;  js<=hsy; ++js)
	  {
	    double summat = 0;
	    // sum over kernel centered on (i,j)
	    for (int jk=-hky; jk<=hky; ++jk)
	      for (int ik=-hkx; ik<=hkx; ++ik)
		summat += vi->Psf(is+hkx-ik, js+hky-jk) 
		  * vi->Kern(ik,jk) 
		  * vi->Weight(is+hkx-ik, js+hky-jk);
	    Mat(galind(is,js), ind) = summat;
	  }

#else

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      double summat = 0;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  DPixel *pkern = &(vi->Kern)(ikstartis,jk-js);
		  DPixel *ppsf = &(vi->Psf)(ikstart,jk);
		  DPixel *pw = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*ppsf) * (*pkern) * (*pw);
		      ++pkern; ++ppsf; ++pw;
		    }
		}
	      Mat(galind(is,js), ind) = summat;
	    }
	}
#endif

      ++fluxind;

    }
}

void SimultaneousFit::fillFluxSky()
{
  //***********************
  // flux-sky matrix terms
  //***********************

  if (printlevel) cout << " Filling flux-sky..." << endl;
  int skyind = 0;
  int fluxind = 0;
  // loop over vignets
  for (VignetFitCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      // sum over smallest between psf and data coordinates
      int hx = vi->Hx();
      int hy = vi->Hy();
      double summat = 0;
#ifndef OPTIMIZED
      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  sumat += vi->Psf(i,j) * vi->Weight(i,j);
#else
      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *ppsf = &(vi->Psf)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summat += (*ppsf) * (*pw);
	      ++ppsf; ++pw;
	    }  
	}
#endif
      // now fill in matrix part
      Mat(fluxstart+fluxind, skystart+skyind) = summat;
      ++fluxind;
    }
}

void SimultaneousFit::fillPosPos()
{
  //*********************************************
  // pos-pos matrix terms and pos vector terms
  //*********************************************
  if (printlevel) cout << " Filling pos-pos..." << endl;

  double sumvecx = 0;
  double sumvecy = 0;
  double summatx = 0;
  double summaty = 0;
  double summatxy = 0;

  // loop over vignets
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      
      // sum over the smallest area between psf derivative and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();

#ifndef OPTIMIZED

      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  {
	    double res = vi->Resid(i,j);
	    double dx = vi->DpDx(i,j);
	    double dy = vi->DpDy(i,j);
	    double w = vi->Weight(i,j);
	    sumvecx += res * dx * w;
	    sumvecy += res * dy * w;
	    summatx += dx * dx * w;
	    summaty += dy * dy * w;
	    summatxy += dx * dy * w;
	  }
#else

      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *presid = &(vi->Resid)(-hx,j);
	  DPixel *pdpdx = &(vi->DpDx)(-hx,j);
	  DPixel *pdpdy = &(vi->DpDy)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvecx += (*presid) * (*pdpdx) * (*pw);
	      sumvecy += (*presid) * (*pdpdy) * (*pw);
	      summatx += (*pdpdx) * (*pdpdx) * (*pw);
	      summaty += (*pdpdy) * (*pdpdy) * (*pw);
	      summatxy += (*pdpdx) * (*pdpdy) * (*pw);
	      ++pw; ++presid; ++pdpdx; ++pdpdy;
	    }
	}    
#endif
    }
  // now fill in the matrix and vector
  Vec(xind) = sumvecx;
  Vec(yind) = sumvecy;
  Mat(xind,xind) = summatx;
  Mat(yind,yind) = summaty;
  Mat(yind,xind) = summatxy;

}

void SimultaneousFit::fillPosGal()
{
  //*********************
  // pos-gal matrix part
  //*********************

  if (printlevel) cout << " Filling pos-gal..." << endl;

  // loop over vignets
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      
      int hx = vi->Hx();
      int hy = vi->Hy();
  
      // the dirac case
      if (vi->IsRefResolution)
	{

	  // loop over psf deriv:sum[ dpdx(x) * dirac(y) ] = dpdx(y)	  
#ifndef OPTIMIZED
	  for (int j=-hy;  j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		Mat(galind(i,j), xind) += vi->DpDx(i,j) * vi->Weight(i,j);
		Mat(galind(i,j), yind) += vi->DpDy(i,j) * vi->Weight(i,j);
	      }
#else
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      DPixel *pdpdx = &(vi->DpDx)(-hx,j);
	      DPixel *pdpdy = &(vi->DpDy)(-hx,j);
	      DPixel *pw = &(vi->Weight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  Mat(galind(i,j), xind) += (*pdpdx) * (*pw);
		  Mat(galind(i,j), yind) += (*pdpdy) * (*pw);
		  ++pdpdx; ++pdpdy; ++pw;
		}
	    }
#endif
	  continue;
	}

      // loop over psf deriv in fitting coordinates
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

#ifndef OPTIMIZED

      for (int is=-hsx;  is<=hsx; ++is)
	for (int js=-hsy;  js<=hsy; ++js)
	  {
	    double summatx = 0;
	    double summaty = 0;
	    // sum over kernel centered on (i,j)
	    for (int jk=-hky; jk<=hky; ++jk)
	      for (int ik=-hkx; ik<=hkx; ++ik)
		{
		  double dum = vi->Kern(ik,jk) *vi->Weight(is+hkx-ik,js+hky-jk);
		  summatx += vi->DpDx(is+hkx-ik,js+hky-jk) * dum;
		  summaty += vi->DpDy(is+hkx-ik,js+hky-jk) * dum;
		}
	    Mat(galind(is,js), xind) += summatx;
	    Mat(galind(is,js), yind) += summaty;	      
	  }

#else

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      double summatx = 0;
	      double summaty = 0;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  DPixel *pkern = &(vi->Kern)(ikstartis,jk-js);
		  DPixel *pdpdx = &(vi->DpDx)(ikstart,jk);
		  DPixel *pdpdy = &(vi->DpDy)(ikstart,jk);
		  DPixel *pw = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summatx += (*pdpdx) * (*pkern) * (*pw);
		      summaty += (*pdpdy) * (*pkern) * (*pw);
		      ++pkern; ++pdpdx; ++pdpdy; ++pw;
		    }
		}
	      Mat(galind(is,js), xind) += summatx;
	      Mat(galind(is,js), yind) += summaty;	      
	    }
	} // end of loop on pixels
#endif
    } //end of loop on vignets
}

void SimultaneousFit::fillPosSky()
{
  //***********************
  // pos-sky matrix terms
  //***********************

  if (printlevel) cout << " Filling pos-sky..." << endl;
  // loop over vignets
  int skyind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const VignetFit *vi = *it;
      if (!vi->IsStarHere) continue;
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      double summatx = 0;
      double summaty = 0;

#ifndef OPTIMIZED
      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  {
	    summatx += vi->DpDx(i,j) * vi->Weight(i,j);
	    summaty += vi->DpDy(i,j) * vi->Weight(i,j);
	  }  
#else
      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *pdpdx = &(vi->DpDx)(-hx,j);
	  DPixel *pdpdy = &(vi->DpDy)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summatx += (*pdpdx) * (*pw);
	      summaty += (*pdpdy) * (*pw);
	      ++pdpdx; ++pdpdy; ++pw;
	    }  
	}
#endif
      // now fill matrix part
      Mat(skystart+skyind, xind) = summatx;
      Mat(skystart+skyind, yind) = summaty;
    }
}

void SimultaneousFit::fillGalGal()
{

  //******************************************
  // gal-gal matrix terms and gal vector terms
  //******************************************

  if (printlevel) cout << " Filling gal-gal..." << endl;
  
  // loop over vignets
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();

      // the dirac case
      if (vi->IsRefResolution)
	{
	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		Vec(galind(i,j)) += (vi->Resid)(i,j) * (vi->Weight)(i,j);
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

#ifndef OPTIMIZED
      for (int is=-hsx;  is<=hsx; ++is)
	for (int js=-hsy;  js<=hsy; ++js)
	  {
	    double sumvec = 0;
	    // sum over kernel, stay in fitting coordinates
	    for (int jk=-hky; jk<=hky; ++jk)
	      for (int ik=-hkx; ik<=hkx; ++ik)
		{
		  sumvec += vi->Kern(ik,jk) 
		    * vi->Resid(is+hkx-ik, js+hky-jk) 
		    * vi->Weight(is+hkx-ik, js+hky-jk);
		}
	    Vec(galind(is,js)) += sumvec;
	  }	
#else
      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      double sumvec = 0;
	      KERNIND(hky,hy,js,jkstart,jkend);
	      
	      // sum over kernel, stay in fitting coordinates
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  DPixel *pkern = &(vi->Kern)(ikstartis,jk-js);
		  DPixel *presid = &(vi->Resid)(ikstart,jk);
		  DPixel *pw = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      sumvec += (*pkern) * (*presid) * (*pw);
		      ++pkern; ++presid; ++pw;
		    }
		}
	      Vec(galind(is,js)) += sumvec;
	    }
	}
#endif
    }

  if (!refill) 
    {
      LaIndex galinds(galstart,galend);
      Mat(galinds, galinds) = MatGal;
      return;
    }
  
      /* galaxy-galaxy matrix terms: longest loop.
	 use the simplified relation, which in one dimension can be written as:
	 matrix(m,n) = sum_ik ker[ik-(in-im)] * ker(ik) * weight(ik+im). 
	 This is also a convolution, but again, borders are messy. */
  for (VignetFitCIterator it = begin(); it != end(); ++it)
    {
      const VignetFit *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      // the dirac case
      if (vi->IsRefResolution)
	{

	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		int ind = galind(i,j);
		Mat(ind,ind) += (vi->Weight)(i,j);
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
	      double summat = 0;
	      for (int jk=-hky; jk<=hky; ++jk)
		for (int ik=-hkx; ik<=hkx; ++ik)
		  {
		    // the ugly test (see for comments below)
		    if ((ik+im >= -hx) && (ik+im <= hx) && 
			(jk+jm >= -hy) && (jk+jm <= hy) && 			
			(ik+imn >= -hkx) && (ik+imn <= hkx) &&
			(jk+jmn >= -hky) && (jk+jmn <= hky))
		      {
			summat += (vi->Kern)(ik,jk)
			  * (vi->Kern)(ik+imn,jk+jmn)
			  * (vi->Weight)(ik+im,jk+jm);
		      }	    
		  }
	      Mat(galstart+m,galstart+n) += summat; 	
	    }
	}
    }

  LaIndex galinds(galstart,galend);
  MatGal = Mat(galinds,galinds);
  refill = false;
}


void SimultaneousFit::fillGalSky()
{
  //**********************
  // gal-sky matrix terms
  //*********************

  if (printlevel) cout << " Filling gal-sky..." << endl;
  // loop over vignets
  int skyind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const VignetFit *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      int ind = skystart+skyind;

      // dirac case
      if (vi->IsRefResolution)
	{
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      DPixel *pw = &(vi->Weight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  Mat(galind(i,j), ind) = *pw;
		  ++pw;
		}
	    }
	  continue;
	}

      // loop over smallest area among fitting size and weight
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

#ifndef OPTIMIZED
      for (int is=-hsx;  is<=hsx; ++is)
	for (int js=-hsy;  js<=hsy; ++js)
	  {
	    double summat = 0;
	    for (int jk=-hkx; jk<=hky; ++jk)
	      for (int ik=-hky; ik<=hkx; ++ik)
		{
		  summat += vi->Kern(ik,jk) * vi->Weight(is+hkx-ik,js+hky-jk);
		}
	    Mat(galind(is,js), ind) = summat;
	  }
#else
      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      double summat = 0;
	      // sum over kernel
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  DPixel *pkern = &(vi->Kern)(ikstartis,jk-js);
		  DPixel *pw = &(vi->Weight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*pkern) * (*pw);
		      ++pkern; ++pw;
		    }
		}
	      // now fill matrix elements
	      Mat(galind(is,js), ind) = summat;
	    }
	}
#endif
    }
}

void SimultaneousFit::fillSkySky()
{
  //*********************************************
  // sky-sky matrix terms and flux vector terms
  //*********************************************
  if (printlevel)   cout << " Filling sky-sky..." << endl;

  // loop over vignets
  int skyind = 0;
  for (VignetFitCIterator it = begin(); it != end(); ++it, ++skyind)
    {
      const VignetFit *vi = *it;
      
      // sum over vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      double sumvec = 0;
      double summat = 0;
#ifndef OPTIMIZED
      for (int j=-hy; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i) 
	  {
	    sumvec += vi->Resid(i,j) * vi->Weight(i,j);
	    summat += vi->Weight(i,j);
	  }
#else

      for (int j=-hy; j<=hy; ++j)
	{
	  DPixel *pw = &(vi->Weight)(-hx,j);
	  DPixel *presid = &(vi->Resid)(-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvec += (*presid)* (*pw);
	      summat += (*pw);
	      ++presid; ++pw;
	    }
	}

#endif
      // now fill out the matrix and vector
      int ind = skystart+skyind;
      Vec(ind) = sumvec;
      Mat(ind,ind) = summat;
    }
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::    Solving routines   ::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/

bool SimultaneousFit::Solve(const double &lambda, const bool invert)
{
  if (printlevel) 
    {
      cout << " Solving normal equations..." << endl;
      cout << " lambda = " << lambda << " invert = " << invert << endl;
    }

  if (!invert || lambda > 0) for (int i=0; i<nparams; ++i) Mat(i,i) *= 1.0+lambda;

  int inversion = SpdSolveAndInvert(Mat, Vec, invert);

  if (inversion > 0) 
    {      
      cerr << " SimultaneousFit::Solve Matrix does not look pos.def " << endl;
      return false;
    }

  if (inversion < 0) 
    {
      cerr << " SimultaneousFit:: Solve something went wrong during solving" << endl;
      return false;
    }

  return true;
}

void SimultaneousFit::MakePsfs()
{
  VignetRef->MakeNormalizedPsf();
  // one round to build the convolved psf

  for (VignetFitIterator it=begin(); it != end(); ++it)
    {
      VignetFit *vi = *it;
      if (vi == VignetRef) continue;
      if (vi->IsRefResolution) vi->MakeNormalizedPsf();
      else vi->MakeConvolvedPsf(VignetRef->Psf, VignetRef->DpDx, VignetRef->DpDy);
      //vi->MakeNormalizedPsf();
      vi->DpDx *= vi->flux;
      vi->DpDy *= vi->flux;
    }

  VignetRef->DpDx *= VignetRef->flux;
  VignetRef->DpDy *= VignetRef->flux;
}

void SimultaneousFit::MakeResidsAndWeights()
{ 
  if (printlevel) cout << " Updating weights " << endl;
 for (VignetFitIterator it = begin(); it != end(); ++it) 
    {
      VignetFit *vi = *it;
      vi->MakeModel(galaxy);
      vi->MakeResid();
      // if more than 3 iterations, recompute weights with model
      /*
	if (iter > 3)
	{
	//	  vi->UpdateWeight(iter > 4);
	vi->UpdateWeight(false); // we should use weight maps now.
	refill = true;
	}
	//      else vi->MakeInitialWeight();
	*/
    }
}

double SimultaneousFit::ComputeChi2()
{
  double totchi2 = 0;
  double totresid = 0;

  for (VignetFitIterator it=begin(); it != end(); ++it)
    {
      VignetFit *vi = *it;
      vi->MakeChi2();
      totresid += vi->residsigma;
      totchi2 += vi->chi2;
    }

  if (printlevel) cout << " Total mean residual " << totresid/ndata << endl;
  return totchi2;
}

bool SimultaneousFit::ApplyCorrections(const double Factor)
{
  if (printlevel) cout << " Applying corrections " << endl;

  double precision = 0.01;
  bool status = false;
  if (fit_flux)
    {
      int fluxind = fluxstart;
      for (VignetFitIterator it=begin(); it != end(); ++it)
	{
	  VignetFit *vi = *it;
	  if (vi->IsStarHere) 
	    {
	      double delta_flux = Vec(fluxind)*Factor;
	      vi->flux += delta_flux; 
	      ++fluxind;
	      if (printlevel) cout << " delta_flux " << delta_flux << endl;
	    }
	}
    }

  // corrections in x and y
  // update psf's at the end with new values of fluxes
  if (fit_pos)
    {
      // compute correction in x
      double delta_xc = Vec(xind)*Factor;
      VignetRef->dxc += delta_xc;
      status = fabs(delta_xc) < precision;


      // compute correction in y
      double delta_yc = Vec(yind)*Factor;
      VignetRef->dyc += delta_yc;
      status = fabs(delta_yc) < precision;

      if (printlevel) cout << " delta_xc delta_yc " << delta_xc << " " << delta_yc << endl;

      // apply same corrections as for VignetRef 
      // these positions are not used during the fit, 
      // except if the PSFs on each image are produced without convolving by the kernel
      for (VignetFitIterator it=begin(); it != end(); ++it)
	{
	  VignetFit *vi = *it;
	  if (vi != VignetRef)
	    {
	      vi->dxc += delta_xc; 
	      vi->dyc += delta_yc;	
	    }
	}      
      MakePsfs();
    }

  // correction in sky similar to flux  
  if (fit_sky)
    {
      int skyind = skystart;
      for (VignetFitIterator it=begin(); it != end(); ++it, ++skyind) 
	{
	  VignetFit *vi = *it;
	  double delta_sky = Vec(skyind)*Factor;
	  vi->sky += delta_sky; 
	  status = fabs(delta_sky) < vi->sky*precision;
	  if (printlevel) cout << " delta_sky " << delta_sky << endl;
	}
    }

  // fitting galaxy
  if (fit_gal)
    {
      double totalgal = 0;
      for (int j=-hfy; j<=hfy; ++j) 
	for (int i=-hfx; i<=hfx; ++i)
	  {	    
	    double delta_gal = Vec(galind(i,j))*Factor;
	    (*galaxy)(i,j) += delta_gal;
	    totalgal += delta_gal;
	  }
      status =  (totalgal < galaxy->sum()*precision);
      if (printlevel) cout << " delta_gal " << totalgal << endl;
    }

  MakeResidsAndWeights();

  return status;
}

inline void operator *= (const LaGenMatDouble &A, const double &alpha)
{
  int M = A.size(0);
  int N = A.size(1);
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j)  A(i,j) *= alpha;
}

void SimultaneousFit::GetChi2(double &newchi2, const double &oldchi2, double &lambda)
{
  ApplyCorrections(1);
  
  newchi2 = ComputeChi2();
  if (printlevel) cout << " new chi2 = " << newchi2 << endl;

  switch (Minim)
    {

    case LevenbergMarquardt:	
      {
	double factor = 10;
	if (newchi2 > oldchi2) 
	  {
	    if (printlevel) cout << " Chi2 increased. Multiply LM parameter by " << factor << endl;
	    lambda *= factor;
	    ApplyCorrections(-1);
	  }
	else lambda /= factor;
	break;
      }

    case NewtonRaphson: default:	
      {
	if (newchi2 > oldchi2) 
	  {
	    if (printlevel) cout << " Chi2 increased. Apply dichotomy " << endl;
	    lambda = 0;
	    double fact = 1;
	    while ((newchi2 > oldchi2) && (fact > 0.001))
	      {
		ApplyCorrections(-fact);
		fact /= 2;
		ApplyCorrections(fact);
		newchi2 = ComputeChi2();
		if (printlevel) cout << " dichotomy applied. new chi2 = " << newchi2 << endl;
	      }

	    if ((newchi2 > oldchi2) &&  (fact < 0.002))
	      {
		ApplyCorrections(-fact);
		newchi2 = ComputeChi2();
		if (printlevel) cout << " dichotomy unapplied. new chi2 = " << newchi2 << endl;
	      }
	    // cout << " After dichotomy new chi2 = " << newchi2 << endl;
	  }
      }
    }
}

bool SimultaneousFit::IterateAndSolve(const int MaxIterations, const double Epsilon)
{
  if (printlevel) 
    {
      cout << setiosflags(ios::fixed);
      cout << " Starting iterations "
	   << "\n   Initial SN fluxes and sky =  ";
      for (VignetFitCIterator it = begin(); it != end(); ++it) 
	cout << "\n    " << (*it)->Name() << ' ' << (*it)->flux << ' ' << (*it)->sky << endl;
      cout <<   "   Initial SN xc yc  =  " << VignetRef->ic+VignetRef->dxc 
	   << ' ' << VignetRef->jc+VignetRef->dyc
  	   << "\n   Initial Gal flux  =  " << galaxy->sum();
	cout << endl;
    }

  iter = 0;
  double oldchi2 = 1e30;
  double newchi2 = ComputeChi2();
  double lambda = 0.001; // dumping factor for levenberg marquardt method
  double dof = ndata-nparams;
  bool inverted = false;

  do
    {
      cout << " Iteration # " << iter << " chi2/dof = " << setprecision(4) << newchi2/dof << endl;
      // check if we are finished
      /*      if (((newchi2 < oldchi2) && (oldchi2 - newchi2 < Epsilon*dof) || 
	      (iter == MaxIterations-1)) && (iter>=3))*/
      if (( (fabs(oldchi2/newchi2)- 1 < Epsilon*dof) || 
	    (iter == MaxIterations-1) ) && (iter>=3))
	{
	  cout << " Fit has converged. One last iteration to invert the matrix " << endl;
	  lambda = 0;
	  inverted = true;
	}

      // fill in the matrices
      if (printlevel) cout << " Fill and solve " << endl;
      FillMatAndVec();
      
      if (!Solve(lambda,inverted)) return false;

      // get the proper chi2 and compute the corrections
      if (oldchi2 > newchi2) oldchi2 = newchi2;
      GetChi2(newchi2, oldchi2, lambda);

      ++iter;
    }
  while (!inverted);

  // print out mean pixel residuals per vignets
  if (printlevel) 
    for (VignetFitCIterator it = begin(); it != end(); ++it) 
      cout << **it << endl;

  // rescale covariance matrix with estimated global scale factor
  chi2 = newchi2;
  double scalor = chi2/dof;
  Mat *= scalor;
  
  if (printlevel) 
    {
      cout << " IterateAndSolve successul \n"
	   << "\n   Final SN fluxes and sky =  ";
      for (VignetFitCIterator it = begin(); it != end(); ++it) 
	cout << "\n    " << (*it)->Name() << ' ' << (*it)->flux << ' ' << (*it)->sky << endl;
      cout <<   "   Final SN xc yc  =  " << VignetRef->ic+VignetRef->dxc 
	   << ' ' << VignetRef->jc+VignetRef->dyc
  	   << "\n   Final Gal flux  =  " << galaxy->sum();
      cout << endl;
    }
  return true;
}

void SimultaneousFit::AssignStar()
{
  if (fit_pos)
    {
      VignetRef->star->x = VignetRef->dxc + VignetRef->ic;
      VignetRef->star->y = VignetRef->dyc + VignetRef->jc;
      VignetRef->star->varx = Mat(xind,xind);
      VignetRef->star->vary = Mat(yind,yind);
      VignetRef->star->covxy = Mat(xind,yind);
      // ignore positional error on other vignets??
    }
  
  int skyind = skystart;
  int fluxind = fluxstart;
  
  for (VignetFitIterator it = begin(); it != end(); ++it, ++skyind)
    {
      VignetFit *vi = *it;
      if (fit_flux && vi->IsStarHere)
	{
	  vi->star->flux = vi->flux;
	  vi->star->varflux = Mat(fluxind,fluxind);
	  ++fluxind;
	}
      else if (!vi->IsStarHere) 
	{
	  vi->star->flux = 0;
	  vi->star->varflux = 0;
	}
      if (fit_sky)
	{
	  vi->star->sky = vi->sky;
	  vi->star->varsky = Mat(skyind,skyind);
	}
    }
}

void SimultaneousFit::DoTheFit(const double &MaxScale)
{
  printlevel = (getenv("FITPRINT") != NULL);
  if (printlevel) cout << " Starting the fit" << endl;  
  clock_t tstart = clock(); 
  
  cout << " Initialize the fit " << endl;
  MakeInitialModel();
  //  for (VignetFitIterator it = begin(); it != end(); ++it) (*it)->writeAllFits("toto");
  //  galaxy->writeFits("galaxy_toto.fits");

  if (!(fit_pos || fit_sky || fit_flux || fit_gal)) 
    {
      cerr << " SimultaneousFit::DoTheFit has nothing to fit" << endl;
      return;
    }

  /* If position is to be fitted, iterate more (non linear) 
     and set minimum size of galaxy to be fitted.
     In that case, sky is worthless to fit, initial estimate 
     is usually good enough */
  chi2 = 1e30;
  if (fit_pos)
    {
      bool oldfit_sky = fit_sky;
      fit_sky = false;
      if (fit_gal) Resize(minscale);
      if (!IterateAndSolve(30)) 
	{
	  cerr << " SimultaneousFit::DoTheFit failed" << endl;
	  return;
	}
      AssignStar();
      fit_sky = oldfit_sky;
    }

  /* Freeze position, properly resize the vignet to at least one full fwhm and
     fit linearly the rest. 5 iterations to update the weights 
     and remove cosmics */

  bool oldfit_pos = fit_pos;
  fit_pos = false;

  if ((scale != MaxScale || chi2 > 1e29))
    {
      if (printlevel) cout << " Freeze position and resize vignets at maximum size " << endl;
      Resize(MaxScale);
      MakePsfs();
      MakeResidsAndWeights();
      if (!IterateAndSolve(7))
	{
	  cerr << " SimultaneousFit::DoTheFit failed" << endl;
	  return;
	} 
      AssignStar();
    }

  fit_pos = oldfit_pos;
  clock_t tend = clock();
  if (printlevel) cout << " CPU comsumed for the SimultaneousFit " 
		       <<  float(tend-tstart) / float(CLOCKS_PER_SEC) << endl;
}

void SimultaneousFit::FindMinimumScale(const double &WorstSeeing)
{
  int hmin = max(int(ceil(WorstSeeing*2.3548)),5);
  int hrefx = VignetRef->HSizeX();
  int hrefy = VignetRef->HSizeY();
  int hkx = 0;
  int hky = 0;

  for (VignetFitIterator it = begin(); it != end(); ++it)
    {
      VignetFit *vi = *it;
      if ((vi->Kern.HSizeX() > hkx) || (vi->Kern.HSizeY()> hky))
	{
	  hkx = vi->Kern.HSizeX();
	  hky = vi->Kern.HSizeY();
	  hrefx = max(hmin,hkx) + hkx;
	  hrefy = max(hmin,hky) + hky;
	}
    }

  minscale = double(min(hrefx,hrefy))/double(min(VignetRef->HSizeX(),VignetRef->HSizeY()));
  if (printlevel) cout << " Minimum scaling factor = " << minscale << endl;
}

void SimultaneousFit::MakeInitialModel()
{
  Resize(1);
  MakePsfs();
  // if (fit_gal && (!galaxy));

  // produce initial galaxy: best resolution image, sky and psf subtracted
  if (!VignetRef->IsStarHere)
    {
      int hpx = VignetRef->Psf.HSizeX();
      int hpy = VignetRef->Psf.HSizeY();
      if ((hpx==0) || (hpy==0)) VignetRef->MakeNormalizedPsf();
      for (int j=-hfy; j<=hfy; ++j) 
	for (int i=-hfx; i<=hfx; ++i) 
	  {
	    double star = 0;
	    if (i<=hpx && i>=-hpx && j<=hpy && j>=-hpy) 
	      star = VignetRef->flux * VignetRef->Psf(i,j);
	    (*galaxy)(i,j) = (*VignetRef)(i,j) - VignetRef->sky - star;
	  }
    }
  else 
    {
      for (int j=-hfy; j<=hfy; ++j) 
	for (int i=-hfx; i<=hfx; ++i) 
	  (*galaxy)(i,j) = (*VignetRef)(i,j) - VignetRef->sky;
    }

  MakeResidsAndWeights();
}

void SimultaneousFit::write(const string &MiddleName) const
{
  if (Mat.size(0) != 0 && Mat.size(1) != 0) 
    {
      string file = "cov_"+MiddleName+".dat";
      ofstream stream(file.c_str());
      stream << setiosflags(ios::fixed);
      stream << Mat(LaIndex(fluxstart,fluxend),LaIndex(fluxstart,fluxend));
    }
  galaxy->writeFits("galaxy_"+MiddleName+".fits");
}

double SimultaneousFit::GetGalaxyFlux(double &VarGalFlux) const
{
  if (!galaxy) 
    { 
      cerr << "SimultaneousFit::GetGalaxyFlux no galaxy available." << endl;
      return 0;
    }

  double galflux = 0;
  double sum = 0;

  for (int j=-hfy; j<=hfy; ++j)
    for (int i=-hfx; i<=hfx; ++i)
      {
	int ind = galind(i,j);
	double pixelVar = Mat(ind,ind);
	double psf = VignetRef->Psf(i,j);
	double weight = psf / pixelVar;
	galflux += (*galaxy)(i,j) * weight;
	sum += weight * psf;
      }

  galflux /= sum;
  VarGalFlux = 1 / sum;

  return galflux;
}

#ifdef UN_JOUR_CA_MARCHERA

// could not make it work. was designed to replace the ugly test in the fillGalGal loop above

#define MAX3(A,B,C) (A>B ? (A>C ? A : C) : (B>C ? B : C))
#define MIN3(A,B,C) (A<B ? (A<C ? A : C) : (B<C ? B : C))

	  int ikstart = MAX3(-hkx, -hkx-imn, -hsx-im);
	  int jkstart = MAX3(-hky, -hky-jmn, -hsy-jm);
	  int ikend = MIN3(hkx, hkx-imn, hsx-im);
	  int jkend = MIN3(hky, hky-jmn, hsy-jm);
          ikstart = (ikstart < ikend) ? ikstart : ikend;  
          jkstart = (jkstart < jkend) ? jkstart : jkend;  
	  summat = 0;
	  for (int jk=jkstart; jk<=jkend; ++jk)
	    {
	      DPixel *pkerm = &(vi->Kern)(ikstart,jk);
	      DPixel *pkern = &(vi->Kern)(ikstart+imn,jk+jmn);
	      DPixel *pw = &(vi->Weight)(ikstart+im,jk+jm);
	      for (int ik=ikstart; ik<=ikend; ++ik)
		{
		  summat += (*pkerm) * (*pkern) * (*pw); 
		  ++pkerm; ++pkern; ++pw;
		}
	    }
	  Mat(galstart+m,galstart+n) = summat; 	
	}
    }
}
#endif
