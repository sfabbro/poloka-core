#include <psfmatch.h>
 
#include "simfitvignet.h"

// for kernel persistence
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"

//#define FNAME
//#define DEBUG

#define VALCUTOFF 0.

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

 cout << "before Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;

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
  cout << "after Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;
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
  

  double integrale = 0;

  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(ic+i,jc+j, Pt, *ppdx, *ppdy);
	integrale += *ppsf;
      }
  
  if(!(integrale>0.5)) {
    cout << "Very strange psf integral in TabulatedPsf::Tabulate (1) : " << integrale << endl;
    cout << "better quit (first write psf as fits)" << endl;
    writeFits("psf_DEBUG.fits");
    abort();
  }

  // normalization
  double norme = 1./integrale;
  ppsf = begin();
  ppdx = Dx.begin();
  ppdy = Dy.begin();
  
  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf *= norme;
	*ppdx *= norme;
	*ppdy *= norme;
      }
}

void TabulatedPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
{
  Resize(Rect.Hx(), Rect.Hy());

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  double integrale = 0;

  for (int j=Rect.ystart; j<Rect.yend; ++j)
    for (int i=Rect.xstart; i<Rect.xend; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(i,j, Pt, *ppdx, *ppdy);
	integrale += *ppsf;
      }

  if(!(integrale>0.5)) {
    cout << "Very strange psf integral in TabulatedPsf::Tabulate (2) : " << integrale << endl;
    cout << "better quit (first write psf as fits)" << endl;
    writeFits("psf_DEBUG.fits");
    abort();
  }

  // normalization
  double norme = 1./integrale;
  ppsf = begin();
  ppdx = Dx.begin();
  ppdy = Dy.begin();
  
  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf *= norme;
	*ppdx *= norme;
	*ppdy *= norme;
      }
}

void TabulatedPsf::Scale(const double& s)
{  

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  for (int i=Nx()*Ny(); i; --i)
    {
      *ppsf *= s;
      *ppdx *= s;
      *ppdy *= s;
      ++ppsf; ++ ppdx; ++ppdy;
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
#ifdef FNAME
  cout << " > SimFitRefVignet::Resize(const int Hx, const int Hy)" << endl;
#endif
  
  if  (Hx < 0 || Hy < 0) 
    {
      cerr << " SimFitRefVignet::Resize(" << Hx << "," << Hy << ") : impossible \n";
      return;
    } 
  
  // resize Data, Weight, Resid if necessary
  Vignet::Resize(Hx,Hy);
  
  // resize Psf, Psf.Dx, Psf.Dy if necessary
  Psf.Tabulate(*Star,*psf,Hx);

  // resize Galaxy 
  if(Galaxy.Nx()!=Data.Nx() || Galaxy.Ny()!=Data.Ny())
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
      //*pgal = *pdat - Star->sky - Star->flux * *ppsf;
      *pgal = 0.;
      ++pgal; ++pdat; ++ppsf;
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
	//*ppsf = psf->Value(i+ic, j+jc, *Star, *ppdx, *ppdy); 
	*pres = *pdat - Star->flux * *ppsf - *pgal - Star->sky; // JG 
	++ppsf; ++ppdx; ++ppdy; ++pres; ++pdat; ++pgal; // JG
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
  inverse_gain = 1./Rim->Gain();
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
  inverse_gain = 1./Rim->Gain();
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

#ifndef ONEPSFPERIMAGE
  if(!psf_updated)
    BuildPsf();
#endif
  
#ifdef DEBUG
  // check sizes
  //int hx_kernel = Kern.HSizeX();
  //int hy_kernel = Kern.HSizeY();
  int hx_data = Data.HSizeX();
  int hy_data = Data.HSizeY();
  int hx_weight = Weight.HSizeX();
  int hy_weight = Weight.HSizeY();
  int hx_resid = Resid.HSizeX();
  int hy_resid = Resid.HSizeY();
  int hx_psf = Psf.HSizeX();
  int hy_psf = Psf.HSizeY();
#ifndef ONEPSFPERIMAGE
  if (hx_data !=  hx_weight || hx_data != hx_resid || hx_data != hx_psf
      || hy_data !=  hy_weight || hy_data != hy_resid || hy_data != hy_psf ) {
#else
 if (hx_data !=  hx_weight || hx_data != hx_resid
      || hy_data !=  hy_weight || hy_data != hy_resid ) {
#endif
    cout << "   SimFitVignet::Update ERROR wrong size " << endl;
    cout << Kern << endl;
    cout << Data << endl;
    cout << Weight << endl;
    cout << Resid << endl;
    cout << Psf << endl;
    abort();
  }
#endif

  if(!resid_updated) { 
#ifdef ONEPSFPERIMAGE
    BuildPsf(); // this loads a daophot psf
    if(UseGal)
      UpdateResid_gal();
    else
      UpdateResid();
#else
    if(UseGal) 
      UpdateResid_psf_gal(); // compute psf as convolution of reference psf*kernel
    else
      UpdateResid_psf();
#endif
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
    OptWeight.Allocate(Nx(),Ny());
    kernel_updated = false;
    psf_updated    = false;
    resid_updated  = false;
  }
  Update();
}

void  SimFitVignet::BuildPsf() {
#ifdef FNAME
  cout << " > SimFitVignet::BuildPsf()" << endl;
#endif

#ifndef ONEPSFPERIMAGE
  
  if(!VignetRef->psf) {
    cout << "SimFitVignet::BuildPsf ERROR VignetRef has no psf" <<endl;
    abort();
  }
  Psf.Tabulate(*Star,*(VignetRef->psf),*this);
  // now scale this according to kernel integral
  if(!kernel_updated)
    BuildKernel();
  double photomratio = Kern.sum();
  Psf*=photomratio;
  
#else
  
  if(!psf) {
    psf = new DaoPsf(*rim);
  }
  Psf.Tabulate(*Star,*psf,*this);

#endif
  
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
	cout << "   Reading kernel " << kernelpath << " ..." << endl;
	obj_input<xmlstream> oi(kernelpath);
	oi >> *kernel;
	oi.close();
	cout << "   done" << endl;
	psfmatch.SetKernelFit(kernel);
	cout << "   kernel loaded" << endl;
      }
    else
      {
	cout << " SimFitVignet::BuildKernelPsf() : cannot find kernel " 
	     << kernelpath << ", so we do it" << endl;
	psfmatch.FitKernel(false);
	obj_output<xmlstream> oo(kernelpath);
	oo << *(psfmatch.GetKernelFit());
	oo.close();
      }
  }

  psfmatch.KernelToWorst(Kern, Star->x, Star->y);
  Star->photomratio = Kern.sum();
  kernel_updated = true;
}

void SimFitVignet::UpdateResid_psf_gal()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_psf_gal() : convolving Psf, Galaxy and updating residuals " << endl;
  cout << "Ref.Galaxy(0,0) = " << VignetRef->Galaxy(0,0) << endl;
#endif
  
  SimFitRefVignet& Ref = *VignetRef;

  // convolve all of them at same time  
  double sump, sumx, sumy, sumg;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy, *prgal;
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy, ++pw, ++pow)
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
		  sump += (*pkern) * (*prpsf); 
		  sumx += (*pkern) * (*prpdx);
		  sumy += (*pkern) * (*prpdy);
		  sumg += (*pkern) * (*prgal);
		  --prpsf; --prpdx; --prpdy; --prgal;
		}
	    }
	  
	  if( (!(sump>0)) && (!(sump<=0))) {
	    cout << "ERROR nan with sump in UpdateResid_psf_gal" << sump << endl;
	    DumpDebug();
	    abort();
	  }
	  val = Star->flux*sump+sumg+Star->sky;
	  *pres = *pdat - val;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
	  if(*pw == 0)
	    *pow = 0;
	  else {
	    if(val>VALCUTOFF)
	      *pow = 1./(1./(*pw)+val*inverse_gain);
	    else
	      *pow = *pw;
	  }
	  ++pdat;
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
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;
  
  const TabulatedPsf& RefPsf = VignetRef->Psf;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy, ++pw, ++pow)
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
		  sump += (*pkern) * (*prpsf); 
		  sumx += (*pkern) * (*prpdx); 
		  sumy += (*pkern) * (*prpdy);
		  --prpsf; --prpdx ; --prpdy;
		}
	    }
	  val =  Star->flux * sump + Star->sky;
	  *pres = *pdat - val;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
	   if(*pw == 0)
	    *pow = 0;
	   else {
	     if(val>VALCUTOFF)
	       *pow  = 1./(1./(*pw)+val*inverse_gain);
	     else
	       *pow = *pw;
	   }
	   ++pdat;
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
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;
  const Kernel& RefGal = VignetRef->Galaxy;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data  (-hx,j);
      pres = &Resid (-hx,j);
      ppsf = &Psf   (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++pw, ++pow)
	{
	  sumg = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prgal = &RefGal(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik)
		{
		  sumg += (*pkern) * (*prgal);
		  ++pkern; --prgal;
		}
	    }
	  val = Star->flux* *ppsf + sumg + Star->sky;
	  *pres = *pdat - val;
	  if(*pw == 0)
	    *pow = 0;
	  else {
	    if(val>VALCUTOFF)
	      *pow = 1./(1./(*pw)+val*inverse_gain);
	    else
	      *pow = *pw;
	  }
	  ++pdat;  ++ppsf; 
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
  DPixel *pw = Weight.begin();
  DPixel *pow = OptWeight.begin();
  double val;
  for (int i=Nx()*Ny(); i; --i)
    {
      val = Star->flux * *ppsf + Star->sky;
      *pres = *pdat - val;
       if(*pw == 0)
	 *pow = 0;
       else {
	 if(val>VALCUTOFF)
	   *pow  = 1./(1./(*pw)+val*inverse_gain);
	 else
	   *pow=*pw;
       }
       ++pres; ++pdat; ++ppsf; ++pw; ++pow;
    }
   resid_updated = true;
}



void SimFitVignet::DumpDebug() const {  
  cout << "############# SimFitVignet::DumpDebug #############" << endl;
  cout << " star flux = " << Star->flux << endl; 
  cout << " galaxy(0,0) = " << VignetRef->Galaxy(0,0) << endl;
  cout << " Data(0,0) = " << Data(0,0) << endl;
  cout << " Psf(0,0) = " << Psf(0,0) << endl;
  cout << " writing *_DEBUG.fits" << endl;
  Weight.writeFits("weight_DEBUG.fits");
  Resid.writeFits("resid_DEBUG.fits");
  Data.writeFits("data_DEBUG.fits");
  Psf.writeFits("psf_DEBUG.fits");
  VignetRef->Psf.writeFits("refpsf_DEBUG.fits");
  Kern.writeFits("kern_DEBUG.fits");
}



double SimFitVignet::Chi2() const {
  
  double chi2 = 0;
  DPixel *pow=OptWeight.begin(), *pres=Resid.begin();
  for (int i=Nx()*Ny(); i; --i, ++pres)
    chi2 += *pow++ * (*pres) * (*pres);
  
  if(chi2>=0)
    return chi2;
  
  // now debug 
  cout << "############# SimFitVignet::Chi2 ERROR chi2=" << chi2 << " #############" << endl;
  DumpDebug();
  cout << " now I quit ..." << endl;
  abort();
}
 



//=============================================================
