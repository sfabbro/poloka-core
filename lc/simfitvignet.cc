#include <psfmatch.h>

#include "simfitvignet.h"

// for kernel persistence
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"

#define FNAME
#define DEBUG

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
#ifdef FNAME
  cout << " > TabulatedPsf::Resize(const int Hx, const int Hy)" << endl;
#endif
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
#ifdef FNAME
  cout << " > TabulatedPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius)" << endl;
#endif

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
SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim)
  : Vignet(Rim), psf(new DaoPsf(*Rim))
{
#ifdef FNAME
  cout << " > SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim)" << endl;
#endif
}

SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim, const int Radius)
  : Vignet(Rim, Radius), psf(new DaoPsf(*Rim))
{
#ifdef FNAME
  cout << " > SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim, const int Radius)" << endl;
#endif
}

SimFitRefVignet::SimFitRefVignet(const PhotStar *Star, const ReducedImage *Rim, const int Radius)
  : Vignet(Star, Rim, Radius), psf(new DaoPsf(*Rim)), Psf(*Star, *psf, Radius)
{
#ifdef FNAME
  cout << " > SimFitRefVignet(const PhotStar *Star, const ReducedImage *Rim, const int Radius)" << endl;
#endif
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
  
  // resize Psf, Psf.Dx, Psf.Dy if necessary
  //Psf.Resize(Hx,Hy);
  Psf.Tabulate(*Star,*psf,Hx);

  // resize Galaxy 
  makeInitialGalaxy();
}

void SimFitRefVignet::Load(const PhotStar *Star)
{
#ifdef FNAME
  cout << " > SimFitRefVignet::Load(const PhotStar *Star)" << endl;
#endif
  if (!Star) return;
  Vignet::Load(Star);
  
  if (!psf) psf = new DaoPsf(*rim);
  Psf.Tabulate(*Star, *psf, hx);
  makeInitialGalaxy(); 
  UpdatePsfResid();
}

void SimFitRefVignet::makeInitialGalaxy()
{  
#ifdef FNAME
  cout << " > SimFitRefVignet::makeInitialGalaxy()" << endl;
#endif

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
#ifdef FNAME
  cout << " > SimFitVignet::UpdatePsfResid() : updating residuals with galaxy" << endl;
#endif
  // re-allocate Resid if too small

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
SimFitVignet::SimFitVignet(const ReducedImage *Rim,  SimFitRefVignet* Ref)
  : Vignet(Rim)
{
  VignetRef = Ref;
  kernel_updated = false;
  psf_updated = false;
  resid_updated = false;
  FitFlux = false;
  FitPos = false;
  UseGal = false;
}

SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim,   SimFitRefVignet* Ref)
  : Vignet(Star, Rim, 0), FitFlux(false)
{  
#ifdef FNAME
  cout << " > SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim,  const SimFitRefVignet& Ref)" << endl;
#endif
  VignetRef = Ref;
  kernel_updated = false;
  psf_updated = false;
  resid_updated = false;
  FitFlux = false;
  FitPos = false;
  UseGal = false;
  //Update();
}

void SimFitVignet::AutoResize() {
#ifdef FNAME
  cout << " > SimFitVignet::AutoResize()" << endl;
#endif
  if(!Star) {
    cerr << "SimFitVignet::AutoResize ERROR you need a star to update this vignet, use SetStar for this" << endl;
  }
  int hx_ref = VignetRef->Hx();
  int hy_ref = VignetRef->Hy();
  if(!kernel_updated)
    BuildKernel();
  int hx_kernel = Kern.HSizeX();
  int hy_kernel = Kern.HSizeY();

#ifdef DEBUG
  SimFitRefVignet *toto = VignetRef;
  cout << "   in SimFitVignet::AutoResize VignetRef = " << toto << endl;
  cout << "   in SimFitVignet::AutoResize hx_ref    = " << hx_ref << endl;
  cout << "   in SimFitVignet::AutoResize hy_ref    = " << hy_ref << endl;
  cout << "   in SimFitVignet::AutoResize hx_kernel = " << hx_kernel << endl;
  cout << "   in SimFitVignet::AutoResize hy_kernel = " << hy_kernel << endl;  
#endif
  Resize(hx_ref-hx_kernel,hy_ref-hy_kernel); // resize vignet, this will actually read the data
}

void SimFitVignet::Update()
{
#ifdef FNAME
  cout << " > SimFitVignet::Update()" << endl;
#endif
  if(!Star) {
    cerr << "SimFitVignet::Update ERROR you need a star to update this vignet, use SetStar for this" << endl;
  }
  if(!kernel_updated)
    BuildKernel();
  if(!psf_updated)
    BuildPsf();
  if(!resid_updated) {
    if(FitPos && UseGal) UpdateResid_psf_gal();
    else if(FitPos && (!UseGal))  UpdateResid_psf();
    else if((!FitPos) && UseGal)  UpdateResid_gal();
    else  UpdateResid(); 
  }
}

// make sure that psf, vignet and data have proper sizes.
// allocate if needed, do only with a Star
void SimFitVignet::Resize( int Hx_new,  int Hy_new)
{
#ifdef FNAME
  cout << " > SimFitVignet::Resize( int Hx_new,  int Hy_new)" << endl;
#endif
  if (!Star) return;
  
  // resize data, weight, resid if necessary
  if(Hx()==Hx_new && Hy()==Hy_new) {
#ifdef DEBUG
    cout << "    actually do not resize" << endl;
#endif
  }else{
    Vignet::Resize(Hx_new,Hy_new);
    kernel_updated = false;
    psf_updated    = false;
    resid_updated  = false;
  }
  Update();
}

void  SimFitVignet::BuildPsf() {
  if(!VignetRef->psf) {
    cerr << " SimFitVignet::BuildPsf ERROR no VignetRef->psf to tabulate !!" << endl;
  }
  Psf.Tabulate(*Star,*(VignetRef->psf),*this);
  psf_updated = true;
}

void SimFitVignet::BuildKernel()
{

#ifdef FNAME
  cout << " > SimFitVignet::BuildKernel(const ReducedImage* Ref)" << endl;
#endif
  if (!Star) return;
  
  // build kernel
  PsfMatch psfmatch(*(VignetRef->Image()), *rim);
  {
    const string kernelpath = rim->Dir()+"/kernel_from_"+VignetRef->Image()->Name()+".xml";
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
  kernel_updated = true;
}

void SimFitVignet::UpdateResid_psf_gal()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_psf_gal() : convolving Psf, Galaxy and updating residuals " << endl;
#endif
  
  SimFitRefVignet& Ref = *VignetRef;

  // convolve all of them at same time  
  double sump, sumx, sumy, sumg;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy, *prgal;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  
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
  resid_updated = true;
}
  
// update psf and residuals with a convolved psf
void SimFitVignet::UpdateResid_psf()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_psf() : convolving Psf and updating residuals " << endl;
#endif
  // convolve all of them at same time  
  double sump, sumx, sumy;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();

  const TabulatedPsf& RefPsf = VignetRef->Psf;

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
   resid_updated = true;
}


void SimFitVignet::UpdateResid_gal()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_gal() : convolving galaxy and updating residuals " << endl;
#endif
#ifdef DEBUG
  cout << " in SimFitVignet::UpdateResid Kern  = " << Kern  << endl;
  cout << " in SimFitVignet::UpdateResid Data  = " << Data  << endl;
  cout << " in SimFitVignet::UpdateResid Resid = " << Resid << endl;
  cout << " in SimFitVignet::UpdateResid Psf   = " << Psf   << endl;
  cout << " in SimFitVignet::UpdateResid Star->flux = " << Star->flux << endl;
#endif

  double sumg;
  DPixel *pdat, *pres, *ppsf, *pkern, *prgal;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();

  const Kernel& RefGal = VignetRef->Galaxy;

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
  resid_updated = true;
}


// update residuals with existing psf, no galaxy
void SimFitVignet::UpdateResid()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid() : updating residuals, no galaxy" << endl;
#endif
  
  DPixel *pdat = Data.begin(), *pres = Resid.begin(), *ppsf = Psf.begin();

  

  for (int i=Nx()*Ny(); i; --i)
    {	
      *pres++ = *pdat++ - Star->flux * *ppsf++ - Star->sky;
    }
   resid_updated = true;
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
