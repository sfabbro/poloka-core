#include "psfstar.h"
#include "image.h"
#include "imagepsf.h"
#include "starlistexception.h"

#include <cmath>
#include <assert.h>

#define DEBUG 0


using namespace std;

PSFStar::PSFStar() : fluxmax(0)
{
  psfX = 0;
  psfY = 0;
  psfChi2 = 0;
}

PSFStar::PSFStar(const SEStar &S) : BaseStar(S), fluxmax(S.Fluxmax())
{
  psfX = x;
  psfY = y;
  psfChi2 = 0;
}

PSFStar::PSFStar(const double X, const double Y, const double Flux)
{
  x = X;
  y =Y;
  flux = Flux;
  fluxmax = 0;

  psfX = X;
  psfY = Y;
  psfChi2 = 0;
}

void PSFStar::SetPSFParams(const Vect &Params)
{
  psfParams = Params;
  psfParamsWeight.allocate(Params.Size(), Params.Size());
}


static double sqr(const double &x) { return x*x;}

/* Auxiliary class for Fit :compute chi2, and its first and
   second derivatives. A class here avoids too long and repeated
   argument lists ... */
struct FitStuff
{
  const int starti;
  const int startj;
  const int endi; const int endj;
  const double Gain;
  const Image &I; const Image &W;
  const ImagePSF &PSF;
  const bool fitPos;
  const bool fitParams;

  FitStuff(int Si, int Sj, int Ei, int Ej,
	   const double g,
	   const Image &i, const Image &w,
	   const ImagePSF &psf, const bool Fpos, const bool FitParams) :
    starti(Si), startj(Sj), endi(Ei), endj(Ej),
    Gain(g),
    I(i), W(w),
    PSF(psf), fitPos(Fpos), fitParams(FitParams) {};

  double ComputeAandB(const Vect &Params,  const double FluxForWeight,
		      Mat *A=NULL, Vect *B=NULL) const;
};

double FitStuff::ComputeAandB(const Vect &Params, const double FluxForWeight,
			      Mat *A, Vect *B) const
{
  // The idea here is that the layout of A and B does not depend on
  // what is to be fitted. So, we manipulate "big" matrices (~ 6x6) even if 
  // only fitting pos and flux. 

  unsigned npar = PSF.NPar();
  Vect gradPar(npar); // will remain zero if (!fitParams)
  Vect gradPos(2); // will remain zero if (!fitPos)
  unsigned nparTot = (B)? B->Size(): 1;
  Vect h(nparTot);

  const double flux = Params(npar);
  const double x = Params(npar+1);
  const double y = Params(npar+2);

  double chi2 = 0;
  
  if (A) A->Zero();  if (B) B->Zero();

  for (int j=startj; j <endj; ++j)
    for (int i=starti ; i < endi; ++i)
      {
	double w = W(i,j);
	if (w<=0) continue;
	double psfVal;
	if (fitParams) 
	  psfVal = PSF.PSFValueWithParams(x,y, i, j, Params, 
					  fitPos ? &gradPos : NULL, 
					  fitParams ? &gradPar : NULL);
	else
	  psfVal = PSF.PSFValue(x,y, i, j, 
				fitPos ? &gradPos : NULL, 
				fitParams ? &gradPar : NULL);

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
	    double factDerWeight = +0.5*FluxForWeight*sqr(res*w)/Gain;
	    Vect gradw(nparTot);
	    for (unsigned k = 0; k < npar; ++k)
	      {
		h(k) = gradPar(k) * flux;
		gradw(k) = gradPar(k) * factDerWeight;
	      }
	    h(npar) = psfVal;
	    h(npar+1) = gradPos(0) * flux;
	    h(npar+2) = gradPos(1) * flux;
	    //	    gradw(npar) = 0; useless
	    gradw(npar+1) = gradPos(0) * factDerWeight;
	    gradw(npar+2) = gradPos(1) * factDerWeight;
	    for (unsigned ki=0; ki<nparTot; ++ki)
	      {
		(*B)(ki) += h(ki)*res*w + gradw(ki);
		for (unsigned kj=ki; kj<nparTot; ++kj)
		  {
		    (*A)(ki,kj) += h(ki)*h(kj)*w;
		  }
	      }
	    if (!fitPos) // make the matrix "invertible "
	      { // the corresponding B terms are already zero 
		(*A)(npar+1,npar+1)=1; (*A)(npar+2,npar+2)=1;
	      }
	    if (!fitParams) // make the matrix "invertible "
	      for (unsigned k=0; k< npar; ++k) (*A)(k,k) = 1;
	    if (false && psfVal>0.03) // DEBUG
	      {
		Vect pars(Params);
		double eps = 1e-4;
		pars(0) += eps;
		double newpsfVal = PSF.PSFValueWithParams(x,y, i, j, pars, NULL, NULL);
		double newres = (double(I(i,j)) - flux*newpsfVal);
		double neww = 1./(1./W(i,j)+newpsfVal*FluxForWeight/Gain);
		double numgrad = (sqr(newres)*neww-res*res*w)/eps;
		cout << " dpix grad,numgrad " << h(0)*res*w + gradw(0) << ' ' << -0.5*numgrad << endl;
	      }
	  }
      }
  return chi2;
}

static double sign(const double& a, const double& b) {
  if(b>0) return fabs(a);
  return -fabs(a);
}

//goodies for brent
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define CGOLD 0.3819660
#define ZEPS 1.e-60


/*! Yet another Gauss-Newton minimizer. It fits both the analytical
  PSF parameters, but also the position and flux. Nothing very clever,
  and in particular, there is nothing to fit several (close) stars
  at a time, as DAOPHOT does. Not a (serious) problem for the CFHTLS deep fields. */
bool PSFStar::do_the_fit(const Image &I, const Image &W, 
			 const ImagePSF &PSF, 
			 const bool FitPos, const bool FitPars)
{
  int starti,endi,startj,endj;
  PSF.StampLimits(psfX, psfY, starti, endi, startj, endj);

  unsigned npar = PSF.NPar();
  double gain = PSF.Gain();

  unsigned nparTot = npar+1 /* flux */ + 2 /*position */;
  
#if (DEBUG>=1)
  cout << " fitting star at " << x << ' ' << y << ' ' << " f = " << flux << endl;
#endif


  Vect B(nparTot);
  Mat A(nparTot, nparTot);
  Mat fitWeight;


  // initialize the vector of fitted parameters
  Vect params(nparTot);
  if (FitPars)
    for (unsigned k=0; k < npar; ++k) params(k) = psfParams(k);
  params(npar) = flux;
  params(npar+1) = psfX;
  params(npar+2) = psfY;


  int niter = 0;
  int maxiter = 20;

  /* we carry out 2 minimizations because the flux enters into the error computation (see
  ComputeAandB), and it should not be accounted for in the derivatives. So we do a first "rough" fit
  to get the flux right, and a final fit using the "correct" flux. In between we reinitialize
  in particular the chi2 "monitoring".
  */
  double diff[2] = {0.1,0.01}; // chi2 difference to stop.


  FitStuff fitStuff(starti, startj, endi, endj,  
		    gain, I, W, PSF, FitPos, FitPars);

  for (int loop=0; loop<=1; ++loop)
    {
      double minDiff = diff[loop];
      double oldChi2=1e30;
      psfChi2 = 1e20;
      double fluxForWeight = flux;
      if (niter>=maxiter) break;
      niter = 0;
#if (DEBUG>=2)
      cout << " starting iteration loop with fluxForWeight = " << fluxForWeight << endl;
#endif
      do { // minimization loop 
	oldChi2 = psfChi2;
	
	// use fluxForWeight for computing variances rather than psfFlux to avoid rising chi2.
	double c0 = fitStuff.ComputeAandB(params, fluxForWeight, &A, &B);
	if (0) { //DEBUG
	  Mat A(nparTot, nparTot);
	  Vect B(nparTot);
	  fitStuff.ComputeAandB(params, fluxForWeight, &A, &B); // get b

	  double eps[6]={1e-4,1e-4,1e-4,1,1e-3,1e-3};
	  for (unsigned k=0; k<6; ++k)
	    {
	      Vect pars(params);
	      pars(k)+= eps[k];
	      double cplus = fitStuff.ComputeAandB(pars, fluxForWeight, NULL,NULL);
	      pars(k) -= 2*eps[k];
	      double cminus = fitStuff.ComputeAandB(pars, fluxForWeight, NULL,NULL);
	      cout << " der : k , grad, numgrad " 
		   <<  k << ' ' 
		   << -2*B(k) << ' ' 
		   << 0.5*(cplus-cminus)/eps[k] << ' '
		   << endl;
	    }
	}// END DEBUG	
	// actually solve
	if (cholesky_solve(A,B,"U") != 0) 
	  {
	    psfChi2 = 1e30; return false;
	  }
	fitWeight = A; // to extract covariance once at minimum.



	// trying to bracket a minimum along the found direction : 
	Vect dir = B;
	double x0 = 0;
	double x2 = 2;
	double c2;
	do {
	  c2 = fitStuff.ComputeAandB(params+x2*dir, fluxForWeight);
	  if (isnan(c2)) x2*=0.6;
	}while (isnan(c2));
	double x1 = 0.5*x2;
	double c1 = fitStuff.ComputeAandB(params+x1*dir, fluxForWeight);


	// first stopping test :
	if (niter>0 && x1 == 1 && fabs(c1-c0)<0.1*minDiff) break;
#ifdef STORAGE
	for (int kk=0; kk<10; ++kk)
	  {
	    if (c1<c0) break;
	    x1 -= 0.6*(x1-x0);
	    c1 =  fitStuff.ComputeAandB(params+x1*dir, fluxForWeight);
	  }
	if (c1>c0)
	  {
	    x1 = -1;
	    c1 =  fitStuff.ComputeAandB(params+x1*dir, fluxForWeight);
	    if (c0<c1)
	      {
		swap(c0,c1);
		swap(x0,x1);
	      }
	    else
	      {
		cout << " PSFStar::Fit : it seems that the derivative is going the wrong way" << endl;
		cout << " x0 x1 x2 c0 c1-c0 c2-c0 " << x0 << ' ' << x1 << ' ' << x2 << ' ' << c0 << ' ' << c1-c0 << ' ' << c2-c0 << endl;
		cout << " chi2(x=1)-chi2(x=0) " << fitStuff.ComputeAandB(params+dir, fluxForWeight)-c0 
		     << " niter = " << niter << endl;
		return false;
	      }
	  }
	for (int kk=0; kk<10; ++kk)
	  {
	    if (c2>c1) break;
	    x2 += (x2-x1)*0.2;
	    c2 = fitStuff.ComputeAandB(params+x2*dir, fluxForWeight);
	    if (isnan(c2))
	      {
		cout << " PSFStar::Fit : cannot find the upper bound for line search  without causing FPE " << endl;
		return false;
	      }
	  }
#endif
#if (DEBUG>=1)
	cout << " before brent x0 x1 x2 c0 c1-c0 c2-c0 " << x0 << ' ' << x1 << ' ' << x2 << ' ' << c0 << ' ' << c1-c0 << ' ' << c2-c0 << endl;
#endif
	// we are now set up for brent ( borrowed from GSL, I think)
	double cmin = c1;
	double xmin = x1;
	double bx = x1;
	int iter; int maxiter=20;
	double tol = 0.01; // i.e. 1% of the Gauss step.
	double a,b,d,etemp,fu,fv,fw,fx,tol1,tol2,u,v,w,x,xm;
	double e=0.0;  // This will be the distance moved on the step before last.
	a=(x0 < x2 ? x0 : x2); //a and b must be in ascending order,
	b=(x0 > x2 ? x0 : x2); //but input abscissas need  not be.
	x=w=v=bx; // Initializations...
	fw=fv=fx=c1;
	for (iter=0;iter<maxiter;iter++) { // main Brent loop
	  xm=0.5*(a+b);
	  // I (Pierre) don't understand this convergence criterion. So leave the computations 
	  // of tol1 and tol2, but stop on abcissa differences
	  tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	  //	  if (fabs(x-xm) <= (tol2-0.5*(b-a))) { // Test for done here.
	  if (iter>1 && fabs(x-xm) <=  tol){ // Test for done here.
	    cmin = fx;
	    xmin = x;
	      break;
	    }
	    if (fabs(e) > tol1) { // Construct a trial parabolic fit.
	      double r=(x-w)*(fx-fv);
	      double q=(x-v)*(fx-fw);
	      double p=(x-v)*q-(x-w)*r;
	      q=2.0*(q-r);
	      if (q > 0.0) p = -p;
	      q=fabs(q);
	      etemp=e;
	      e=d;
	      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
		d=CGOLD*(e=(x >= xm ? a-x : b-x));
	      /* The above conditions determine the acceptability of the
		 parabolic fit. Here we take the golden section step into the
		 larger of the two segments.
	      */ 
	      else {
		d=p/q; //Take the parabolic step.
		u=x+d;
		if (u-a < tol2 || b-u < tol2)
		  d=sign(tol1,xm-x);
	      }
	    } else {
	      d=CGOLD*(e=(x >= xm ? a-x : b-x));
	    }
	    u=(fabs(d) >= tol1 ? x+d : x+sign(tol1,d));
	    fu= fitStuff.ComputeAandB(params+u*dir, fluxForWeight);
	    /*This is the one function evaluation per iteration.*/
	    if (fu <= fx) { // Now decide what to do with our function evaluatiON
	      if (u >= x) a=x; else b=x; 
      
	      SHFT(v,w,x,u); // Housekeeping follows:
	      SHFT(fv,fw,fx,fu);
	    } else {
	      if (u < x) a=u; else b=u;
	      if (fu <= fw || w == x) {
		v=w;
		w=u;
		fv=fw;
		fw=fu;
	      } else if (fu <= fv || v == x || v == w) {
		v=u;
		fv=fu;
	      }
	    } //Done with housekeeping. Back for
	  }//another iteration. (end of brent loop)      
	if (xmin<x0 || xmin>x2)
	  {
	    cout << " ERROR : brent output outside input interval " << endl;
	    abort();
	  }

	if (iter == maxiter)
	  {
	    cout << " PSFStar::Fit : Brent over maxiter" << endl;
	    return false;
	  }
#if (DEBUG>=1)
	cout << " xmin cmin loop " <<  xmin << ' ' << cmin << ' ' << loop << endl;
#endif

	// update parameters
	psfChi2 = cmin;
	for (unsigned k=0; k < nparTot; ++k) params(k) += xmin*dir(k);
	
	// some sanity checks flux == nan is really bad
	for (unsigned k=0; k < nparTot; ++k) 
	  if (isnan(params(k)))
	    {
	      cout << " PSFStar::Fit : one parameter read nan " << endl;
	      psfChi2 = 1e30;
	      return false;
	    }
	niter ++;
      } // end of minimization loop
      while (niter < maxiter && (fabs(oldChi2-psfChi2)>minDiff));
      fluxForWeight = params(npar);
    } // end of for (int loop=0 ; ...)

#if (DEBUG>=1)
  if (fitPos)
  cout << " FitPSFParams: results dx,dy, f, chi2 "
       << psfX-params(npar+1) << ' ' << psfY-params(npar+2) << ' ' << psfFlux << ' ' << psfChi2 << endl;
  else
  cout << " FitPSFParams: results f, chi2 "
       << psfFlux << ' ' << psfChi2 << endl;
#endif


  /* We have to extract the weight matrix of the psf parameters (the first npar of params vector).
     This involves "marginalization" over the position and flux of the star. Compute 
     full covariance matrix (from chi2 second derivatives), extract the psf params sub block,
     and invert it back to get a weight matrix
  */

  if (cholesky_invert(fitWeight,"U") != 0) return false; // should not happen 
  fitWeight.Symmetrize("U");
  // fitWeight is now a covariance matrix

  // copy fitted parameters in the right places:
  if (FitPars)
    for (unsigned k=0; k < npar; ++k) psfParams(k) = params(k);

  flux = params(npar); 
  eflux = sqrt(fitWeight(npar,npar));

  if (FitPos)
    {
      x = psfX = params(npar+1);
      y = psfY = params(npar+2);
      vx = fitWeight(npar+1, npar+1);
      vy = fitWeight(npar+2, npar+2);
      vxy = fitWeight(npar+1, npar+2); // hope it is what 'U' means in lapack stuff.
      if (fxyCov.SizeX() != 3) fxyCov.allocate(3,3);
      for (unsigned k=0; k<3; ++k) for (unsigned l=0; l<3 ;++l)
	fxyCov(k,l) = fitWeight(npar+k, npar+l);
    }
  else
    { // it should not have changed
      assert(psfX == params(npar+1)); 
    }



  // extract the PSF param part of the covariance matrix (called fitWeight !)
  if (FitPars)
    {
      for (unsigned i=0; i<npar; ++i) 
	for (unsigned j=i; j<npar; ++j)
	  {
	    psfParamsWeight(i,j) = fitWeight(i,j);
	    psfParamsWeight(j,i) = fitWeight(i,j);
	  }
      // turn it back into a weight matrix
      psfParamsWeight.CholeskyInvert("L"); // should never fail
      psfParamsWeight.Symmetrize("L");
    }
  return (niter < maxiter);
}


bool PSFStar::FitAllParams(const Image &I, const Image &W, 
			   const ImagePSF &PSF)
{
  return do_the_fit(I,W,PSF, /*Fitpos */ true, /* FitPsfParams */ true);
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
  return do_the_fit(I,W,PSF, /* Fitpos */ true, false);
}
  
#ifdef STORAGE
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
#warning : in this routine, the order of parameters is fxy, while it is assumed everywhere else that it fxy.
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
#endif




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
      const Mat &cov = s.FXYCov();
      if (Wcs)
	{
	  double ra,dec;
	  Wcs->apply(s.PSFX(), s.PSFY(), ra,dec);
	  file << ra << ' ' << dec << ' ';
	}
      else file << " -1. 100. ";
      file   << s.flux << ' '      
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

string PSFStar::WriteHeader_(ostream &s, const char* i) const
{
  if (i== NULL) i= "";
  string baseStarFormat =  BaseStar::WriteHeader_(s, i);
  s << "# eflux"<< i <<" : Flux uncertainty " << endl;  
  s << "# fluxmax"<< i <<" : Peak pixel value above background " << endl;
  s << "# psfx" << i << " : " << endl;
  s << "# psfy" << i << " : " << endl;
  s << "# psfchi2" << i << " : " << endl;
  s << "# npar : number of psf params" << endl;
  unsigned npar = psfParams.size();  
  for (unsigned k=0; k < npar; ++k) s << "# p"<< k << i <<" : " << endl;
  for (unsigned k=0; k< npar; ++k) 
    for (unsigned l=0; l<=k; ++l) s << "# wp" << k << l << i << " : " << endl;
  char names[4] = "fxy";
  for (unsigned k=0; k<3; ++k)
    for (unsigned l=0; l<=k; ++l) s << "#cov" << names[k]<<names[l]<<i << " : "<< endl; 
  return baseStarFormat+" PSFStar 1";
}


void PSFStar::writen(ostream &s) const
{
  BaseStar::writen(s);
  s << eflux << ' ' << fluxmax << ' '
    << psfX << ' ' << psfY << ' ' 
    << psfChi2 << ' ';
  unsigned npar = psfParams.size();
  s << npar << ' ';
  for (unsigned k=0; k <npar ; ++k)  s << psfParams(k) << ' ';
  for (unsigned i=0; i< npar; ++i)
    for (unsigned j=0; j<=i; ++j) s << psfParamsWeight(i,j) << ' ';
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<=i; ++j) s << fxyCov(i,j) << ' ';
}
    
#include "fastifstream.h"

PSFStar* PSFStar::read(fastifstream& Rd, const char *Format)
{
  PSFStar *s= new PSFStar();
  s->read_it(Rd,Format);
  return s;
}


void PSFStar::read_it(fastifstream& Rd, const char *Format)
{
  int formatValue = 0;
  if (Format) 
    formatValue = DecodeFormat(Format,"PSFStar");
  if (formatValue == 1)
    {
      BaseStar::read_it(Rd , Format);
      unsigned npar;
      Rd >> eflux >> fluxmax >> psfX >> psfY >> psfChi2 >> npar;
      psfParams.allocate(npar);
      for (unsigned k=0; k <npar ; ++k)  Rd >> psfParams(k);
      psfParamsWeight.allocate(npar,npar);
      for (unsigned i=0; i< npar; ++i)
	for (unsigned j=0; j<=i; ++j) 
	  {Rd >> psfParamsWeight(i,j);psfParamsWeight(j,i) = psfParamsWeight(i,j);}
      
      fxyCov.allocate(3,3);
      for (unsigned i=0; i<3; ++i)
	for (unsigned j=0; j<=i; ++j) 
	  { Rd >> fxyCov(i,j); fxyCov(j,i) = fxyCov(i,j);}
    }
 else throw(StarListException(" Unknown format value for PsfStar "));
}  

//instantiate I/O stuff
#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

template class StarList<PSFStar>; 
