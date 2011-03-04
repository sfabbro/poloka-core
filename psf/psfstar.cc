#include "psfstar.h"
#include "image.h"
#include "imagepsf.h"


#include <cmath>

#define DEBUG 0


using namespace std;

PSFStar::PSFStar(const AperSEStar &A) : AperSEStar(A)
{
  psfX = x;
  psfY = y;
  psfFlux = flux;
  oldPsfFlux = 0;
  psfChi2 = 0;
}

PSFStar::PSFStar(const double X, const double Y, const double Flux)
{
  x = X;
  y =Y;
  flux = Flux;

  psfX = X;
  psfY = Y;
  psfFlux = Flux;
  oldPsfFlux = 0;
  psfChi2 = 0;
}

void PSFStar::SetPSFParams(const Vect &Params)
{
  psfParams = Params;
  psfParamsWeight.allocate(Params.Size(), Params.Size());
}


/* Auxiliary routine for FitPSFParams :compute chi2, and its first and
   second derivatives. */
static double ComputeAandB(const int starti, const int startj, 
			   const int endi, const int endj,
			   const Vect &Params,
			   const double FluxForWeight, const double Gain,
			   const Image &I, const Image &W,
			   const ImagePSF &PSF,
			   Mat *A=NULL, Vect *B=NULL)
{
  unsigned npar = PSF.NPar();
  unsigned nparTot = (B)? B->Size(): 1;
  Vect h(nparTot);
  Vect gradPar(npar);
  Vect gradPos(2);
  const double x = Params(npar);
  const double y = Params(npar+1);
  const double flux = Params(npar+2);

  double chi2 = 0;
  
  if (A && B) { A->Zero(); B->Zero();}

  for (int j=startj; j <endj; ++j)
    for (int i=starti ; i < endi; ++i)
      {
	double w = W(i,j);
	if (w<=0) continue;
	double psfVal = PSF.PSFValueWithParams(x,y, i, j, Params, 
					       &gradPos, &gradPar);
	double res = (double(I(i,j)) - flux*psfVal);
	/* w only accounts for the sky noise. add Poisson variance
	   of the star photons. There is a trap here: in the variance
	   associated to the pixel, there should not be any noise
	   correlated to the data: this leads to biassed estimates (in
	   particular of the flux). This is why we use the model value
	   rather than the pixel content to compute this Poisson contribution. */
	w = 1./(1./w+psfVal*FluxForWeight/Gain);
	chi2 += res*res*w;
	if (A && B)
	  {
	    for (unsigned k = 0; k < npar; ++k)
	      h(k) = gradPar(k) * flux;
	    h(npar) = gradPos(0) * flux;
	    h(npar+1) = gradPos(1) * flux;
	    h(npar+2) = psfVal;
	    
	    for (unsigned ki=0; ki<nparTot; ++ki)
	      {
		(*B)(ki) += h(ki)*res*w;
		for (unsigned kj=ki; kj<nparTot; ++kj)
		  {
		    (*A)(ki,kj) += h(ki)*h(kj)*w;
		  }
	      }
	  }
      }
  return chi2;
}



/*! Yet another Gauss-Newton minimizer. It fits both the analytical
  PSF parameters, but also the position and flux. Nothing very clever,
  and in particular, there is nothing to fit several (close) stars
  at a time, as DAOPHOT does. Not a problem for the CFHTLS deep fields. */
bool PSFStar::FitPSFParams(const Image &I, const Image &W, const ImagePSF &PSF)
{

  int starti,endi,startj,endj;
  PSF.StampLimits(psfX, psfY, starti, endi, startj, endj);

  unsigned npar = PSF.NPar();
  double gain = PSF.Gain();



  unsigned nparTot = npar+3;// +3 for x, y, flux
  
#if (DEBUG>=1)
  cout << " fitting star at " << x << ' ' << y << ' ' << " f = " << flux << endl;
#endif


  Vect b(nparTot);
  Mat a(nparTot, nparTot);
  Mat fitWeight;


  // initialize the vector of fitted parameters
  Vect params(nparTot);
  for (unsigned k=0; k < npar; ++k) params(k) += psfParams(k);
  params(npar) = psfX;
  params(npar+1) = psfY;
  params(npar+2) = psfFlux;


  int niter = 0;
  int maxiter = 20;

  /* we carry out 2 minimizations because the flux enters into the error computation (see
  ComputeAandB), and it should not be accounted for in the derivatives. So we do a first "rough" fit
  to get the flux right, and a final fit using the "correct" flux. In between we reinitialize
  in particular the chi2 "monitoring".
  */
  double diff[2] = {0.1,0.01}; // chi2 difference to stop.
  for (int loop=0; loop<=1; ++loop)
    {
      double minDiff = diff[loop];
      double oldChi2=1e30;
      psfChi2 = 1e20;
      double fluxForWeight = psfFlux;
      if (niter>=maxiter) break;
      niter = 0;
#if (DEBUG>=2)
      cout << " starting iteration loop with fluxForWeight = " << fluxForWeight << endl;
#endif
      do { // minimization loop 
	oldChi2 = psfChi2;
	
	// use fluxForWeight for computing variances rather than psfFlux to avoid rising chi2.
	double c0 = ComputeAandB(starti, startj, endi, endj, params, fluxForWeight, gain, I, W, PSF, &a, &b);
	// actually solve
	if (cholesky_solve(a,b,"U") != 0) 
	  {
	    psfChi2 = 1e30; return false;
	  }
	fitWeight = a; // to extract covariance once at minimum.
	
	Vect dir = b;
	double x0 = 0;
	double x2 = 2;
	double c2;
	do {
	  c2 = ComputeAandB(starti, startj, endi, endj, params+x2*dir, fluxForWeight, gain, I, W, PSF);
	  if (isnan(c2)) x2*=0.6;
	}while (isnan(c2));
	double x1 = 0.5*x2;
	double c1 = ComputeAandB(starti, startj, endi, endj, params+x1*dir, fluxForWeight, gain, I, W, PSF);
	double xmin=x1, cmin=c1;
	//    if (c1<c2 && c1<c0) { cmin=c1; xmin = x1;}
	if (c2<c1 && c2<c0) { cmin=c2; xmin = x2;}
	if (c0<c1 && c0<c2) { cmin=c0; xmin = x0;}
	// optimize along  dir
	for (int kk=0;kk<10;++kk)
	  {
	    double xnew = ( (c0-c1)*(x0*x0-x2*x2)-(c0-c2)*(x0*x0-x1*x1) )/ ( (c0-c1)*(x0-x2)-(c0-c2)*(x0-x1) ) / 2.;
	    double cnew = ComputeAandB(starti, startj, endi, endj, params+xnew*dir, fluxForWeight, gain, I, W, PSF);
#if (DEBUG>=2)
	    cout << " c0 c1 c2 cnew " << c0 << ' ' << c1 << ' ' << c2 << ' ' << cnew << endl;
	    cout << " x0 x1 x2 xnew " << x0 << ' ' << x1 << ' ' << x2 << ' ' << xnew << endl;
#endif

	    if (cnew<cmin) { cmin=cnew; xmin=xnew;}
	    
	    if (isnan(cmin))
	      {
		psfChi2 = 1e30; 
		return false;
	      }
	    x0=x1; c0=c1;
	    x1=x2; c1=c2;
	    x2=xnew; c2=cnew;
	    if (fabs(c1-c2)< minDiff) break;
	  }
	psfChi2 = cmin;



#if (DEBUG>=2)
	cout << " niter: " << niter << " chi2: " << psfChi2 
	     << " dchi2 " << oldChi2 - psfChi2 << " step " << xmin
	     << endl;
#endif
	
	// update parameters
	for (unsigned k=0; k < nparTot; ++k) params(k) += xmin*dir(k);
	
	// some sanity checks flux == nan is really bad
	for (unsigned k=0; k < nparTot; ++k) 
	  if (isnan(params(k)))
	    {
	      psfChi2 = 1e30;
	      return false;
	    }

	niter ++;
      } // end of minimization loop
      while (niter < maxiter && (fabs(oldChi2-psfChi2)>minDiff));
      psfFlux = params(npar+2);
    } // end of for (int loop=0 ; ...)

  // copy fitted parameters in the right places:
  for (unsigned k=0; k < npar; ++k) psfParams(k) = params(k);
  psfX = params(npar);
  psfY = params(npar+1);
  psfFlux = params(npar+2);

  
  /* We have to extract the weight matrix of the psf parameters (the first npar of params vector).
     This involves "marginalization" over the position and flux of the star. Compute 
     full covariance matrix (from chi2 second derivatives), extract the psf params sub block,
     and invert it back to get a weight matrix
  */

  if (cholesky_invert(fitWeight,"U") != 0) return false; // should not happen 
  // fitWeight is now a covariance matrix
  // extract the PSF param part:
  for (unsigned i=0; i<npar; ++i) 
    for (unsigned j=i; j<npar; ++j)
      {
	psfParamsWeight(i,j) = fitWeight(i,j);
	psfParamsWeight(j,i) = fitWeight(i,j);
      }
  // turn it back into a weight matrix
  psfParamsWeight.CholeskyInvert("L"); // should never fail
  psfParamsWeight.Symmetrize("L");


#if (DEBUG>=1)
  cout << " results " << psfFlux << ' ' << psfChi2 << endl;
#endif
  return (niter < maxiter);
}


/*! One more Gauss Newton solver for the star params (x,y,flux). No
provision for fitting several stars at a time (as DAOPhot does). Will
do that in an other life. Or this one if we --really-- need it. The
routine FitPSFParams has some code to find the minimum along the
direction found when solving the linear system. Since we enter here
with already accurate positions and fluxes, it is in practice not
necessary.  Would it become mandatory, it is not difficult to
implement it. We could merge the FitStarParams and FitPSFParams
routines by providing a routine in ImagePSF that does both PSFValue
and PSFValueWithParams.*/
bool PSFStar::FitStarParams(const Image &I, const Image &W, 
			    const ImagePSF &PSF)
{
  int starti,endi,startj,endj;
  PSF.StampLimits(psfX, psfY, starti, endi, startj, endj);
  double gain = PSF.Gain();

  double oldChi2;
  int niter = 0;
#if (DEBUG>=1)
  cout << " fitting star params at " << x << ' ' << y << ' ' << " f = " << flux << endl;
#endif
  psfChi2 = 1e30;
  int maxiter = 20;

  do { // fit iterations
    Vect grad(3);
    Vect b(3);
    Mat m(3,3);
    oldChi2 = psfChi2;
    psfChi2 = 0;
    
    for (int j=startj; j <endj; ++j)
      for (int i=starti ; i < endi; ++i)
	{
	  double w = W(i,j);
	  if (w<=0) continue;
	  double psfVal = PSF.PSFValue(psfX, psfY, i, j, &grad, 0);
	  grad *= psfFlux;
	  grad(2) = psfVal;
	  double res = (I(i,j) - psfFlux*psfVal);
	  /* w only accounts for the sky noise. Add Poisson variance
	  of the star photons. There is a trap here: in the variance
	  associated to the pixel, there should not be any noise
	  correlated to the data: this leads to biased estimates (in
	  particular of the flux). This is why we use the model value
	  rather than the data to compute this Poisson contribution.*/
	  w = 1./(1./w+psfVal*psfFlux/gain);
	  psfChi2 += res*res*w;
	  for (unsigned ki=0; ki<3; ++ki)
	    {
	      b(ki) += grad(ki)*res*w;
	      for (unsigned kj=ki; kj<3; ++kj)
		{
		  m(ki,kj) += grad(ki)*grad(kj)*w;
		}
	    }
	}
#if (DEBUG>=2)
    cout << " niter: " << niter << " chi2: " << psfChi2 
	 << " dchi2 " << oldChi2 - psfChi2 
	 << endl;
#endif
    if (cholesky_solve(m,b,"U") != 0) return false;
    xyfCov = m;
#if (DEBUG>=2)
    if (psfChi2 > oldChi2 + 1e-2)
      {
	cout << "PSFStar::FitStarParams : chi2 rising " 
	     << psfChi2-oldChi2<< endl;
      }
#endif
    // update parameters
    psfX += b(0);
    psfY += b(1);
    psfFlux += b(2);
    niter ++;
  }
  while (niter < maxiter && (fabs(oldChi2-psfChi2)>0.01 || oldChi2<=psfChi2+1e-2) && abs(oldChi2-psfChi2)>1e-3);
  // extract the cov matrix
  cholesky_invert(xyfCov,"U");

#if (DEBUG>=2)
  cout << " FitStarParams: results dx,dy, f, chi2 " 
       << psfX-x << ' ' << psfY-y << ' ' << psfFlux << ' ' << psfChi2 << endl;
#endif
  return (niter < maxiter);
}




/************ PSFStarList ******************/

PSFStarList::PSFStarList(const AperSEStarList &L)
{
  for (AperSEStarCIterator i = L.begin(); i != L.end(); ++i)
    push_back( new PSFStar(**i));
}



#include <fstream>
#include <iomanip>
#include "gtransfo.h"

void PSFStarList::WriteTuple(const string &FileName, const Gtransfo* Wcs, const ImagePSF* PSF) const
{
  if (size() == 0 ) 
    {
      cout << " no tuple written to " << FileName << " because list is empty " << endl;
      return;
    }
  ofstream file(FileName.c_str());
  file << "# ra : " << endl;  
  file << "# dec : " << endl;
  file << "# flux : psf flux" << endl;
  file << "# xpsf : xpsf" << endl;
  file << "# ypsf : ypsf" << endl;
  file << "# ex : expsf" << endl;
  file << "# ey : eypsf" << endl;
  file << "# eflux : e psf flux" << endl;
  file << "# covxy : " << endl;
  file << "# covxf : " << endl;
  file << "# covyf : " << endl;
  file << "# xg : gaussx " << endl;
  file << "# yg : gaussy " << endl;
  file << "# fs : Sex f" << endl;
  file << "# fluxmax : max pixel" << endl;
  file << "# dx : xpsf -x " << endl;  
  file << "# dy : ypsf -y" << endl;  
  file << "# chi2psf : " << endl;
  const PSFStar &first = *front();
  int npar = first.PSFParams().Size();
  for (int k = 0 ; k < npar; k++)
    file << "# p" << k << " : psf param " << k << endl;     
  if (PSF) for (int k = 0 ; k < npar; k++)
    file << "# pg" << k << " : psf param from global fit " << k << endl;     
  file << "#end" << endl;
  file << "# format BaseStar 2"<< endl;
  file << setprecision(8);
  for (PSFStarCIterator it = begin(); it != end(); ++it)
    {
      const PSFStar &s = **it;
      const Mat &cov = s.XYFCov();
      if (Wcs)
	{
	  double ra,dec;
	  Wcs->apply(s.PSFX(), s.PSFY(), ra,dec);
	  file << ra << ' ' << dec << ' ';
	}
      else file << " -1. 100. ";
      file   << s.PSFFlux() << ' '      
	     << s.PSFX() << ' '
	     << s.PSFY() << ' ';
      if (cov.SizeX()>=3)
	file << sqrt(cov(0,0)) << ' '
	     << sqrt(cov(1,1)) << ' '
	     << sqrt(cov(2,2)) << ' '
	     << cov(0,1) << ' '
	     << cov(0,2) << ' '
	     << cov(1,2) << ' ';
      else file << " -1 -1 -1 0 0 0 ";
      file << s.x << ' '
	   << s.y << ' '
	   << s.flux << ' '
	   << s.Fluxmax() << ' '
	   << s.PSFX() - s.x << ' '
	   << s.PSFY() - s.y << ' '
	   << s.PSFChi2() << ' ';
      for (int k = 0 ; k < npar; k++)
	file << s.PSFParams()(k) << ' ';     
      if (PSF)
	{
	  Vect p(PSF->PSFParams(s.PSFX(), s.PSFY()));
	  for (int k = 0 ; k < npar; k++)
	    file << p(k) << ' ';
	}
      file << endl;	
    }
  file.close();
}

#include "dicstar.h"

bool PSFStarList::ReadTuple(const string &FileName)
{
  DicStarList l(FileName);
  for (DicStarCIterator i = l.begin(); i != l.end(); ++i)
    {
      const DicStar &s = **i;
      double x = s.getval("xpsf");
      double y = s.getval("ypsf");
      double f = s.getval("flux");
      PSFStar *psfStar = new PSFStar(x,y,f);
      psfStar->x = s.getval("xg");
      psfStar->y = s.getval("yg");
      psfStar->flux = s.getval("fs");
      psfStar->Fluxmax() = s.getval("fluxmax");
      push_back(psfStar);
    }
  return true;  
}
