#include "vignetfit.h"
#include "photstar.h"
#include "night.h"
#include "daophotpsf.h" // for DaoPsf

static double sqr(const double &x) {return x*x;}

VignetFit::VignetFit() 
  : IsStarHere(true),IsRefResolution(false),
    flux(0),sky(0),chi2(0),residsigma(0)
{
}

VignetFit::VignetFit(PhotStar *aStar, const Night *aNight,
		     const Image &Source, const Kernel &aKernel,
		     const int HRefX, const int HRefY)
  :Vignet(aStar->x, aStar->y, Source, HRefX-aKernel.HSizeX(), HRefY-aKernel.HSizeY()),
   night(aNight), star(aStar), IsStarHere(true), IsRefResolution(false), 
   flux(aStar->flux), sky(aStar->sky), chi2(0), residsigma(0), Kern(aKernel)
{
  if (night->backLevel == 0) sky = 0;
}

VignetFit::VignetFit(PhotStar *aStar, const Night *aNight, const string &Name)
  : Vignet(aNight->Name()+"_"+Name+"_data.fits"), night(aNight), star(aStar),
  IsStarHere(true), IsRefResolution(false), flux(aStar->flux), sky(aStar->sky), 
  chi2(0), residsigma(0), Kern(aNight->Name()+"_"+Name+"_kernel.fits")
{
  hx = hSizeX - Kern.HSizeX();
  hy = hSizeY - Kern.HSizeY();
  if (night->backLevel == 0) sky = 0;
}

const double ALPHA = 3;  // clipping factor in sigma residal unit
const double BETA  = 6;  // clipping exponent

// assume weight was already assigned externally or by MakeInitialWeight
void VignetFit::UpdateWeight(const bool robustify)
{

  if (Model.IsEmpty() || Resid.IsEmpty()) 
    {
      cerr << " VignetFit::UpdateWeight missing Model or Resid for " << night->Name() << endl;
      return;
    }

  if (hx != Weight.HSizeX() || hy != Weight.HSizeY()) MakeInitialWeight();

  // loop over residual map to see which pixels are badly modelized (cosmics,defects)

  int hmodx = Model.HSizeX();
  int hmody = Model.HSizeY();

  /* The variance on one pixel is the sum of:

      (1) Readout Noise
      (2) Poisson noise from sky and stars, given by the model
      (3) Flatfielding errors
      (4) Errors made in approximating the sky and the model (galaxy, kernel and psf).
      
    Here we assume that these quantities are given by:

      Var_init(i,j) = ReadoutNoise^2 + Sky(i,j)/Gain + (FlatError(i,j)*Sky(i,j))^2 
                    + (SkyError*Sky(i,j))^2

    which we update with:
      
      Var_update(i,j) = Var_init(i,j) + Model(i,j)/Gain + (FlatError(i,j)*Model(i,j))^2
                      + (ModelError*Model(i,j))^2

    Practically:
      - Var_init(i,j) should be built from the weight map, because it includes 
        cosmics, glitches deweighting. It is constructed by measuring the sky rms, 
	the flat, and masks, and includes the read out noise, sky and the flat error.
      - SkyError is neglected (mostly because never really measured it)
      - ModelError could be more accurate, but is just approximated as ProfileError
  */

  for (int j=-hy; j<=hy; ++j)
    for (int i=-hx; i<=hx; ++i)
      {
	double w = Weight(i,j);
	double varPixel = (w > 1e-10) ? 1./w : 1e20;
	if ((i<=hmodx && i>=-hmodx && j<=hmody && j>=-hmody) && (w > 1e-10))
	  {
	    // take abs because it could be that the model is erronous
	    // model = flux*psf+galaxy+local_sky
	    double vmodel = fabs(Model(i,j)-night->backLevel);

	    vmodel = vmodel/night->gain 
	      + sqr(vmodel * night->flatError)
	      + sqr(vmodel * night->profileError);
	    varPixel += vmodel;
	    Weight(i,j) = 1./varPixel;
	  }
	
	/* robustify weight by depondering high residuals
	   the empirical weighting function was chosen after Stetson chi2 lectures.
	   see http://nedwww.ipac.caltech.edu/level5/Stetson/frames.html */
	if (robustify) Weight(i,j) *= 1. / 
			 (1. + pow(fabs(Resid(i,j))/(sqrt(varPixel)*ALPHA), BETA));
      }
}

void VignetFit::MakeInitialWeight()
{
  if (hx == Weight.HSizeX()&& hy == Weight.HSizeY()) return;


  double sigma2 = sqr(night->sigmaSky/night->gain);

  for (int j=-hy; j<=hy; ++j)
    for (int i=-hx; i<=hx; ++i)
      Weight(i,j) = 1./sigma2;
}

void VignetFit::MakeConvolvedPsf(const Kernel &PsfRef, const Kernel &DpDxRef, const Kernel &DpDyRef)
{
  if (IsRefResolution)
    {
      MakeNormalizedPsf();
      return;
    }

  // get the psf at least of PsfRef convolved size
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  int hrx = PsfRef.HSizeX();
  int hry = PsfRef.HSizeY();

  Kernel bigPsf(PsfRef, hkx, hky);
  Kernel bigDpDx(DpDxRef, hkx, hky);
  Kernel bigDpDy(DpDyRef, hkx, hky);

  if ((Psf.HSizeX() < hrx)  && (Psf.HSizeY() < hry)) Psf = Kernel(hrx,hry);
  if ((DpDx.HSizeX() < hrx) && (DpDx.HSizeY() < hry)) DpDx = Kernel(hrx,hry);
  if ((DpDy.HSizeX() < hrx) && (DpDy.HSizeY() < hry)) DpDy = Kernel(hrx,hry);

  // convolve all of them at same time
  for (int j=-hry; j<=hry; ++j)
    for (int i=-hrx; i<=hrx; ++i)
      {
	double sump = 0;
	double sumx = 0;
	double sumy = 0;
	DPixel *pkern = Kern.begin();
	for (int jk =-hky; jk <= hky; ++jk)
	  {
	    DPixel *ppsf = &bigPsf(i+hkx, j-jk);
	    DPixel *pdx  = &bigDpDx(i+hkx, j-jk);
	    DPixel *pdy  = &bigDpDy(i+hkx, j-jk);
	    for (int ik = -hkx; ik <= hkx; ++ik) 
	      {
		sump += (*pkern) * (*ppsf); 
		sumx += (*pkern) * (*pdx); 
		sumy += (*pkern) * (*pdy); 
		++pkern; --ppsf; --pdx; --pdy;
	      }
	  }
	Psf(i,j)  = sump;
	DpDx(i,j) = sumx;
	DpDy(i,j) = sumy;
      }
}

void VignetFit::MakeResid()
{
  if (Model.IsEmpty())
    {
      cerr << " VignetFit::MakeResid missing Model for " 
	   << night->Name() << endl;
      return;
    }

  if (hx != Resid.HSizeX() || hy != Resid.HSizeY()) Resid = Kernel(hx,hy);

  for (int j=-hy; j<=hy; ++j) 
    for (int i=-hx; i<=hx; ++i) 
      Resid(i,j) = (*this)(i,j) - Model(i,j);
}

void VignetFit::MakeChi2()
{
  if (Resid.IsEmpty() || Weight.IsEmpty())
    {
      cerr << " VignetFit::MakeChi2 missing Resid or Weight for " 
	   << night->Name() << endl;
      return;
    }

  if (hx != Chi2.HSizeX() || hy != Chi2.HSizeY()) Chi2 = Kernel(hx,hy);

  chi2 = 0;
  residsigma = 0;
  for (int j=-hy; j<=hy; ++j) 
    for (int i=-hx; i<=hx; ++i) 
      {
	double resid = Resid(i,j);
	double weight = Weight(i,j);
	double pixelchi2 = weight * resid*resid;
	chi2 += pixelchi2;
	residsigma += resid*sqrt(weight);
	Chi2(i,j) = pixelchi2;
      }
}

void VignetFit::MakeModel(const Vignet *Galaxy)
{
  // combine them all into model kernel vignet

  if (Model.HSizeX() != hx || Model.HSizeY() != hy) Model = Kernel(hx, hy);

  /* 
     The Model Vignet is composed as:
     
     Model(i,j) = flux X [PSF_best*Kernel](i-xc,j-yc) + [Galaxy*Kernel](i,j) + Sky

     where Sky is the sky/pixel under the galaxy+sn:

       Sky = Global_Sky + Local_Sky

       Global_Sky = current sky in the full image (subtracted or nor from the original)
                  = night->backLevel

       Local_Sky = local sky difference from the global sky. 
                   Usually close to zero, but it can be fited.
		  = this->sky  
                  
  */

  if (Galaxy)
    {
      int hgx = Galaxy->HSizeX();
      int hgy = Galaxy->HSizeY();
      if (Kern.IsEmpty()) 
	{
	  for (int j=-hgy; j<=hgy; ++j)
	    for (int i=-hgx; i<=hgx; ++i) 
	      Model(i,j) = (*Galaxy)(i,j); 
	} 
      else ConvolveKernels(Model, *Galaxy, Kern);
    }

  if (IsStarHere)
    {
      for (int j=-hy; j<=hy; ++j) 
	for (int i=-hx; i<=hx; ++i) 
	  Model(i,j) += flux*Psf(i,j); 
    }

  for (int j=-hy; j<=hy; ++j)
    for (int i=-hx; i<=hx; ++i) 
      Model(i,j) += night->backLevel + sky; 
}

double VignetFit::WeightedAperture(double &VarAper)
{
  double sigma2 = night->sigmaSky*night->sigmaSky;
  double aper = fabs(Aperture(VarAper,night->seeing*2.3548,sigma2, sky));
  MakeNormalizedPsf();
  return Vignet::WeightedAperture(VarAper, Psf*aper, sigma2, sky);
}

void VignetFit::QuickPhotometry(const double &Nfwhm)
{
  WeightedRecentroid(star->flux, star->sky, night->seeing, night->seeing);
  star->x = ic + dxc;
  star->y = jc + dyc;

  // compute the covariance of estimated positions with that mean
  double vx = 0;
  double vy = 0;
  double cxy = 0;
  double sumw = 0;
  double alphax = -0.5/sqr(night->sigmaX);
  double alphay = -0.5/sqr(night->sigmaY);

  for (int j=-hy; j<=hy; ++j)
    for (int i=-hx; i<=hx; ++i)
      {
	double pixelVal = (*this)(i,j);
	double w = exp((i-dxc)*(i-dxc)*alphax + (j-dyc)*(j-dyc)*alphay);
	vx += sqr(i*pixelVal*w - dxc);
	vy += sqr(j*pixelVal*w - dyc);
	cxy += (i*pixelVal*w - dxc)*(j*pixelVal*w - dyc);
	sumw += w;
      }
  star->varx = vx / sumw;
  star->vary = vy / sumw;
  star->covxy = cxy / sumw;

  // flux shit
  double radius = max(Nfwhm*night->seeing*2.3548,min(hSizeX-1.0,hSizeY-1.0));
  star->flux = Aperture(star->varflux, radius, sqr(night->sigmaSky), star->sky);
}

void VignetFit::writeAllFits(const string &MiddleName) const 
{
  string fullname = night->Name()+"_"+MiddleName;
  writeFits(fullname+"_data.fits");  
  if (!Kern.IsEmpty()) Kern.writeFits(fullname+"_kernel.fits");
  if (!Weight.IsEmpty()) Weight.writeFits(fullname+"_weight.fits");
  if (!Resid.IsEmpty()) Resid.writeFits(fullname+"_resid.fits");
  if (!Model.IsEmpty()) Model.writeFits(fullname+"_model.fits");
  if (IsRefResolution) Psf.writeFits(fullname+"_psf.fits");
}
string VignetFit::Name() const {return night->Name();}
void VignetFit::readAllFits(const string &MiddleName)
{  
  string fullname = night->Name()+"_"+MiddleName;
  readFits(fullname+"_data.fits");
  Weight.readFits(fullname+"_weight.fits");
  Psf.readFits(fullname+"_psf.fits");
  Chi2.readFits(fullname+"_chi2.fits");
  Resid.readFits(fullname+"_resid.fits");
  Model.readFits(fullname+"_model.fits");
  Kern.readFits(fullname+"_kernel.fits");
  DpDx.readFits(fullname+"_dpdx.fits");
  DpDy.readFits(fullname+"_dpdy.fits");
}

void VignetFit::dump(ostream &Stream) const
{
  int npix = (2*hx+1)*(2*hy+1);
  Stream << " " << night->Name() 
	 << " hx " << hx 
	 << " chi2/pixel " << chi2/npix
	 << " resid/pixel " << residsigma/npix;
}


void VignetFit::MakeNormalizedPsf()
{
  if (Psf.HSizeX() != hx || Psf.HSizeY() != hy) Psf = Kernel(hx, hy);
  if (DpDx.HSizeX() != hx || DpDx.HSizeY() != hy) DpDx = Kernel(hx, hy);
  if (DpDy.HSizeX() != hx || DpDy.HSizeY() != hy) DpDy = Kernel(hx, hy);

  DaoPsf daopsf(night->ImagePsfName());
  double rad2 = sqr(daopsf.Radius());
  double xc = ic+dxc;
  double yc = jc+dyc;

  for (int j=-hy; j<=hy; ++j) 
    {
      double dy2 = sqr(j-dyc);
      for (int i=-hx; i<=hx; ++i) 
	{
	  if (sqr(i-dxc)+dy2 < rad2)
	    {
	      double dpdx = 0 ,dpdy = 0;
	      Psf(i,j) = daopsf.Value(ic+i,jc+j, xc, yc, dpdx, dpdy);
	      DpDx(i,j) = dpdx;
	      DpDy(i,j) = dpdy;
	    }
	}
    }


  /* should be obsolete by now.
  double norm = 1. / Psf.sum() * night->photomRatio;

  // normalize all of them with PSF norm
  DPixel *ppsf = Psf.begin();
  DPixel *pdx = DpDx.begin();
  DPixel *pdy = DpDy.begin();
  for (int i=Psf.Nx()*Psf.Ny(); i; --i)
    {
      (*ppsf) *= norm;
      (*pdx)  *= norm;
      (*pdy)  *= norm;
      ++ppsf; ++pdx; ++pdy;
    }
  */
}
