#include <psfmatch.h>

#include "simfitvignet.h"

// for kernel persistence
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"

//#define DEBUG

TabulatedPsf::TabulatedPsf(const Point& Pt, const DaoPsf& Dao, const int Radius)
  : Kernel(Radius), Dx(Radius), Dy(Radius)
{ 
  Tabulate(Pt, Dao, Radius);
}

TabulatedPsf::TabulatedPsf(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
  : Kernel(Rect.Hx(), Rect.Hy()), Dx(hSizeX, hSizeY), Dy(hSizeX, hSizeY)
{ 
  Tabulate(Pt, Dao, Rect);
}

void TabulatedPsf::Resize(const int Hx, const int Hy)
{
  if (HSizeX() != Dx.HSizeX() ||
      HSizeY() != Dx.HSizeY() ||
      HSizeX() != Dy.HSizeX() ||
      HSizeY() != Dy.HSizeY() ||
      HSizeX() != Hx ||
      HSizeX() != Hy)
    {
      Allocate(2*Hx+1, 2*Hy+1);
      Dx.Allocate(Nx(), Ny());
      Dy.Allocate(Nx(), Ny());     
    }
}

void TabulatedPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius)
{
  Resize(Radius, Radius);

  const int ic = int(Pt.x);
  const int jc = int(Pt.y);

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(ic+i,jc+j, Pt, *ppdx, *ppdy);
      }
}

void TabulatedPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
{
  Resize(Rect.Hx(), Rect.Hy());

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  for (int j=Rect.ystart; j<Rect.yend; ++j)
    for (int i=Rect.xstart; i<Rect.xend; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(i,j, Pt, *ppdx, *ppdy);
      }
}

void TabulatedPsf::Scale(const double& s)
{  

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  for (int i=Nx()*Ny(); i; --i)
    {
      *ppsf++ *= s;
      *ppdx++ *= s;
      *ppdy++ *= s;
    }
}


//=========================================================================================

SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim, const int Radius)
  : Vignet(Rim, Radius), psf(new DaoPsf(*Rim))
{
}

SimFitRefVignet::SimFitRefVignet(const PhotStar *Star, const ReducedImage *Rim, const int Radius)
  : Vignet(Star, Rim, Radius), psf(new DaoPsf(*Rim)), Psf(*Star, *psf, Radius)
{
}

// make sure that psf, vignet and data have proper sizes.
// allocate if needed
void SimFitRefVignet::Resize(const int Hx, const int Hy)
{
  if  (Hx < 0 || Hy < 0) 
    {
      cerr << " SimFitRefVignet::Resize(" << Hx << "," << Hy << ") : impossible \n";
      return;
    } 

  // resize Data, Weight, Resid if necessary
  Vignet::Resize(Hx,Hy);
  Vignet::Allocate(Nx(), Ny());

  // resize Psf, Psf.Dx, Psf.Dy if necessary
  Psf.Resize(Hx,Hy);

  // resize Galaxy 
  makeInitialGalaxy();
}

void SimFitRefVignet::Load(const PhotStar *Star)
{
  if (!Star) return;
  Vignet::Load(Star);
  makeInitialGalaxy(); 
  if (!psf) psf = new DaoPsf(*rim);
  Psf.Tabulate(*Star, *psf, hx);
  UpdatePsfResid();
}

void SimFitRefVignet::makeInitialGalaxy()
{  
  // allocate galaxy as biggest vignet possible, so resizing it is easy.
  Galaxy.Allocate(Data.Nx(), Data.Ny(), 1);

  // initial galaxy: best resolution image - sky - [flux*psf]  
  DPixel *pdat = Data.begin(), *pgal = Galaxy.begin(), *ppsf = Psf.begin();

  for (int i=Nx()*Ny(); i ; --i)
    {
      *pgal++ = *pdat++ - Star->sky - Star->flux * *ppsf++;;
    }
}


void SimFitRefVignet::UpdatePsfResid()
{
  // re-allocate Resid if too small

  cout << " SimFitRefVignet::UpdatePsfResid() : updating residuals with galaxy" << endl;

  DPixel *pdat = Data.begin(), *pres = Resid.begin(), *pgal = Galaxy.begin();
  DPixel *ppsf = Psf.begin(), *ppdx = Psf.Dx.begin(), *ppdy = Psf.Dy.begin();

  for (int j=ystart; j<yend; ++j)
    for (int i=xstart; i<xend; ++i)
      {
	*ppsf = psf->Value(i, j, *Star, *ppdx, *ppdy);
	*pres++ = *pdat++ - Star->flux * *ppsf - *pgal++ - Star->sky; 
	++ppsf; ++ppdx; ++ppdy;
      }
    
}

//=========================================================================================
SimFitVignet::SimFitVignet(const ReducedImage *Rim)
  : Vignet(Rim), FitFlux(false)
{
}

SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim,  const SimFitRefVignet& Ref)
  : Vignet(Star, Rim, 0), FitFlux(false)
{  
  BuildKernel(Ref.Image());
  UpdatePsfResid(Ref);
}

// make sure that psf, vignet and data have proper sizes.
// allocate if needed, do only with a Star
void SimFitVignet::Resize(const int Hx, const int Hy)
{

  if (!Star) return;

  // resize data, weight, resid if necessary
  Vignet::Resize(Hx,Hy);

  // resize psf if necessary
  Psf.Resize(Hx,Hy);

}

void SimFitVignet::BuildKernel(const ReducedImage* Ref)
{

  if (!Star) return;

  // build kernel
  PsfMatch psfmatch(*Ref, *rim);
  {
    const string kernelpath = rim->Dir()+"/kernel_from_"+Ref->Name()+".xml";
    if(FileExists(kernelpath)) 
      {
	KernelFit *kernel = new KernelFit();
	obj_input<xmlstream> oi(kernelpath);
	oi >> *kernel;
	oi.close();
	psfmatch.SetKernelFit(kernel);
      }
    else
      {
	cerr << " SimFitVignet::BuildKernelPsf() : cannot find kernel " 
	     << kernelpath << ", so we do it" << endl;
	psfmatch.FitKernel(false);
	obj_output<xmlstream> oo(kernelpath);
	oo << *(psfmatch.GetKernelFit());
	oo.close();
      }
  }

  psfmatch.KernelToWorst(Kern, Star->x, Star->y);

}

void SimFitVignet::UpdatePsfResid(const SimFitRefVignet& Ref)
{

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

  DPixel *pdat = Data.begin(), *pres = Resid.begin(), *ppsf = Psf.begin();

  cout << " SimFitVignet::UpdateResid() : updating residuals, no galaxy" << endl;

  for (int i=Nx()*Ny(); i; --i)
    {	
      *pres++ = *pdat++ - Star->flux * *ppsf++ - Star->sky;
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
