#include <psfmatch.h>

#include "simfitvignet.h"

// for kernel persistence
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"

static double sqr(const double& x) { return x*x; }

TabulatedPsf::TabulatedPsf(const Point& Pt, const DaoPsf& Dao, const int Radius)
  : Kernel(Radius), Dx(Radius), Dy(Radius)
{ 

  const int ic = int(Pt.x);
  const int jc = int(Pt.y);
  const double dxc = ic - Pt.x;
  const double dyc = jc - Pt.y;
  const double maxrad2 = sqr(Dao.Radius());
  double dy2;
  DPixel *ppsf=begin(), *ppdx=Dx.begin(), *ppdy=Dy.begin();

  for (int j=-hSizeY; j<=hSizeY; ++j) 
    {
      dy2  = sqr(j-dyc);
      for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
	{
	  if (sqr(i-dxc)+dy2 < maxrad2)
	    {
	      *ppsf = Dao.Value(ic+i,jc+j, Pt.x, Pt.y, *ppdx, *ppdy);
	    }
	  else { *ppsf = *ppdx = *ppdy = 0.; }
	}
    }
}

TabulatedPsf::TabulatedPsf(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
  : Kernel(), Dx(hSizeX, hSizeY), Dy(hSizeX, hSizeY)
{ 

  const double maxrad2 = sqr(Dao.Radius());
  double dy2;
  DPixel *ppsf=begin(), *ppdx=Dx.begin(), *ppdy=Dy.begin();
  for (int j=Rect.ystart; j<Rect.yend; ++j)
    {
      dy2  = sqr(j-Pt.y);
      for (int i=Rect.xstart; i<Rect.xend; ++i, ++ppsf, ++ppdx, ++ppdy) 
	{
	  if (sqr(i-Pt.x)+dy2 < maxrad2)
	    {
	      *ppsf = Dao.Value(i,j, Pt.x, Pt.y, *ppdx, *ppdy);
	    }
	  else { *ppsf = *ppdx = *ppdy = 0.; }
	}
    }
}

void TabulatedPsf::Scale(const double& s)
{
  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();
  for (int i=Nx()*Ny(); i; --i, ++ppsf, ++ppdx, ++ppdy)
    {
      *ppsf *= s;
      *ppdx *= s;
      *ppdy *= s;
    }
}


//=========================================================================================

SimFitRefVignet::SimFitRefVignet(const PhotStar *Star, const ReducedImage *Rim, const int Radius)
  : Vignet(Star, Rim, Radius), Psf(*Star, DaoPsf(*Rim), Radius)
{
}

void SimFitRefVignet::makeInitialGalaxy()
{  

  // allocate galaxy as biggest vignet possible, so resizing it is easy.
  Galaxy.Allocate(Data.Nx(), Data.Ny(), 1);

  // initial galaxy: best resolution image - sky - [flux*psf]  
  int hpx = Psf.HSizeX();
  int hpy = Psf.HSizeY();

  for (int j=-hy; j<=hy; ++j) 
    for (int i=-hx; i<=hx; ++i)
      {
	double pixStar=0.;
	if (i<=hpx && i>=-hpx && j<=hpy && j>=-hpy && Star->flux != 0.)
	  pixStar = Star->flux * Psf(i,j);
	Galaxy(i,j) = Data(i,j) - Star->sky - pixStar;
      }
}

void SimFitRefVignet::UpdatePsfResid()
{
  // re-allocate Resid if too small
  if (hx > Resid.HSizeX() || hy > Resid.HSizeY())
    {
      cout << " SimFitRefVignet::UpdatePsfResid() : resizing residuals" << endl;
      Resid = Kernel(hx,hy);
    }

  cout << " SimFitRefVignet::UpdatePsfResid() : updating residuals with galaxy" << endl;

  const double maxrad2 = sqr(psf->Radius());
  const int ic = int(Star->x);
  const int jc = int(Star->y);
  const double dxc = Star->x - ic;
  const double dyc = Star->y - jc;

  DPixel *pdat, *pres, *pgal;
  DPixel *ppsf, *ppdx, *ppdy;

  double dy2;

  for (int j=-hy; j<=hy; ++j) 
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      pgal = &Galaxy (-hx,j);

      dy2  = sqr(j-dyc);

      for (int i=-hx; i<=hx; ++i, ++ppdx, ++ppdy) 
	{
	  if (sqr(i-dxc)+dy2 < maxrad2)
	    {
	      *ppsf = psf->Value(i+ic,j+jc, Star->x, Star->y, *ppdx, *ppdy);
	    }
	  else { *ppsf = *ppdx = *ppdy = 0.; }
	  *pres = *pdat++ - Star->flux * *ppsf++ - *pgal++ - Star->sky; 
	}
    }
}

//=========================================================================================
SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim, const ReducedImage *Ref, 
			   const int Radius)
  : Vignet(Star, Rim, Radius), FitFlux(false)
{
  // load kernel
  PsfMatch psfmatch(*Ref,*Rim);
  {
    std::string kernelpath = Rim->Dir()+"/kernel_from_"+Ref->Name()+".xml";
    if(FileExists(kernelpath.c_str())) {
      KernelFit *kernel = new KernelFit();
      obj_input<xmlstream> oi(kernelpath);
      oi >> *kernel;
      oi.close();
      psfmatch.SetKernelFit(kernel);
    }else{
      std::cerr << "WARNING cannot find kernel " << kernelpath << ", so we do it" << std::endl;
      psfmatch.FitKernel(false);
      obj_output<xmlstream> oo(kernelpath);
      oo << *(psfmatch.GetKernelFit());
      oo.close();
    }
  }
  psfmatch.KernelToWorst(Kern, Star->x, Star->y);

  // make psf=refpsf*kernel
  TabulatedPsf refPsf(*Star, DaoPsf(*Ref), Radius + Kern.HSizeX()); 
  Psf = TabulatedPsf(refPsf.HSizeX()-Kern.HSizeX(), refPsf.HSizeY()-Kern.HSizeY());
  // resize properly vignets ready for convolution
  Resize(refPsf.HSizeX()-Kern.HSizeX(), refPsf.HSizeY()-Kern.HSizeY());
}

void SimFitVignet::UpdatePsfResid(const SimFitRefVignet& Ref)
{
  // resize data if too small
  // if (hx < Gal.HSizeX()-Kern.HSizeX() || hy <  Gal.HSizeY()-Kern.HSizeY()) 
  //  Resize(Gal.HSizeX()-Kern.HSizeX(), Gal.HSizeY()-Kern.HSizeY());
 
  // re-allocate Resid if too small
  if (hx > Resid.HSizeX() || hy > Resid.HSizeY()) 
    {
      cout << " SimFitVignet::UpdatePsfResid() : resizing residuals" << endl;
      Resid = Kernel(hx,hy);
    }

  // re-allocate Psf if too small
  if (Psf.HSizeX() < Ref.Psf.HSizeX()-Kern.HSizeX() || Psf.HSizeY() < Ref.Psf.HSizeY()-Kern.HSizeY()) 
    {
      cout << " SimFitVignet::UpdatePsfResid() : resizing Psf " << endl;
      Psf = TabulatedPsf(Ref.Psf.HSizeX()-Kern.HSizeX(), Ref.Psf.HSizeY()-Kern.HSizeY());
    }

  // convolve all of them at same time  
  double sump, sumx, sumy, sumg;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy, *prgal;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();

  cout << " SimFitVignet::UpdatePsfResid() : convolving Psf, Galaxy and updating residuals " << endl;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy)
	{
	  sumg = sump = sumx = sumy = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prpsf = &Ref.Psf   (i+hkx, j-jk);
	      prpdx = &Ref.Psf.Dx(i+hkx, j-jk);
	      prpdy = &Ref.Psf.Dy(i+hkx, j-jk);
	      prgal = &Ref.Galaxy(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik, ++pkern)
		{
		  sump += (*pkern) * (*prpsf--); 
		  sumx += (*pkern) * (*prpdx--);
		  sumy += (*pkern) * (*prpdy--);
		  sumg += (*pkern) * (*prgal--);
		}
	    }
	  
	  *pres = *pdat++ - Star->flux*sump - sumg - Star->sky;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
	}
    }
}
  
// update psf and residuals with a convolved psf
void SimFitVignet::UpdatePsfResid(const TabulatedPsf& RefPsf)
{

  // re-allocate Resid if too small
  if (hx > Resid.HSizeX() || hy > Resid.HSizeY()) 
    {
      cout << " SimFitVignet::UpdatePsfResid() : resizing residuals" << endl;
      Resid = Kernel(hx,hy);
    }

  // re-allocate Psf if too small
  if (Psf.HSizeX() < RefPsf.HSizeX()-Kern.HSizeX() || Psf.HSizeY() < RefPsf.HSizeY()-Kern.HSizeY()) 
    {
      cout << " SimFitVignet::UpdatePsfResid() : resizing Psf " << endl;
      Psf = TabulatedPsf(RefPsf.HSizeX()-Kern.HSizeX(),RefPsf.HSizeY()-Kern.HSizeY());
    }

  // convolve all of them at same time  
  double sump, sumx, sumy;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();

  cout << " SimFitVignet::UpdatePsfResid() : convolving Psf and updating residuals " << endl;
  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy)
	{
	  sump = sumx = sumy = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prpsf = &RefPsf   (i+hkx, j-jk);
	      prpdx = &RefPsf.Dx(i+hkx, j-jk);
	      prpdy = &RefPsf.Dy(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik, ++pkern) 
		{
		  sump += (*pkern) * (*prpsf--); 
		  sumx += (*pkern) * (*prpdx--); 
		  sumy += (*pkern) * (*prpdy--);
		}
	    }
	  
	  *pres = *pdat++ - Star->flux * sump - Star->sky;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
	}
    }
}


void SimFitVignet::UpdateResid(const Kernel &RefGal)
{
  // resize data if too small
  // if (hx < Gal.HSizeX()-Kern.HSizeX() || hy <  Gal.HSizeY()-Kern.HSizeY()) 
  //  Resize(Gal.HSizeX()-Kern.HSizeX(), Gal.HSizeY()-Kern.HSizeY());
 
  // re-allocate Resid if too small
  if (hx > Resid.HSizeX() || hy > Resid.HSizeY()) 
    {
      cout << " SimFitVignet::UpdateResid() : resizing residuals" << endl;
      Resid = Kernel(hx,hy);
    }

  double sumg;
  DPixel *pdat, *pres, *ppsf, *pkern, *prgal;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();

  cout << " SimFitVignet::UpdateResid() : convolving galaxy and updating residuals " << endl;
  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data  (-hx,j);
      pres = &Resid (-hx,j);
      ppsf = &Psf   (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres)
	{
	  sumg = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prgal = &RefGal(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik)
		{
		  sumg += (*pkern++) * (*prgal--);
		}
	    }	  
	  *pres = *pdat++ - Star->flux* *ppsf++ - sumg - Star->sky;
	}
    }
}


// update residuals with existing psf, no galaxy
void SimFitVignet::UpdateResid()
{
  // re-allocate Resid if too small
  if (hx > Resid.HSizeX() || hy > Resid.HSizeY()) 
    {
      cout << " SimFitVignet::UpdateResid() : resizing residuals" << endl;
      Resid = Kernel(hx,hy);
    }

  DPixel *pdat, *pres, *ppsf;

  cout << " SimFitVignet::UpdateResid() : updating residuals" << endl;
  for (int j=-hy; j<=hy; ++j) 
    {
      pdat = &Data  (-hx,j);
      pres = &Resid (-hx,j);
      ppsf = &Psf   (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres) 
	{	
	  *pres = *pdat++ - Star->flux * *ppsf++ - Star->sky;
	}
    }
}

ostream& operator << (ostream& Stream, const SimFitVignet &Vig)
{
  Stream << " SimFitVignet: " << endl
	 << "     Vignet: " << (Vignet) Vig << endl
	 << "     Psf   : " << Vig.Psf      << endl
	 << "     Kernel: " << Vig.Kern     << endl;
  return Stream;
}


//=============================================================
