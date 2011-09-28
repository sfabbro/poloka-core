
#include <cmath>
#include "matvect.h"
#include "histo2d.h"


#include "histopeakfinder.h"

static double sq(const double &x) { return x*x;}



/* one more star finder:
1) locate a peak (the maximum) in fwhm%log(flux/fluxmax) histogram

2) try to figure out the span of the peak by looking at the bin contents
   around the peak.

3) "fit" a 2D gaussian on the data itself, using an iterative approach:

   using a guess of the peak extension (vxx, vyy, vxy),
   you may compute weighted second moments (sum(x**2*w)/sum(w)),
   with a gaussian w. You don't get the right answer, but if w
   has the same width as the data, the outcome is a variance
   which is half of the true one. There are several inversions
   of 2D symetric matrices done on the flight to switch from covariance 
   to weight matrix. 


   Once we have located the peak and its extension
   in the plane, we can finally select objects whithin an elliptic
   aperture around the star peak, and call them stars.
   If you don't need the stars, you may however consider the Fwhm
   that comes out from this routine: it is somehow better that the
   clipped average done in seeing_box.cc

*/

/*! models a 2D gaussian from the list contents. The histogram
   is assumed to be filled, and the peak location is expected 
   around XGuess and YGuess. The distance to the average,
   in number of sigmas is written into the input List.
   The Ellipse parameters can be retreived in the returned Ellipse
*/

//#define DEBUG

bool HistoPeakFinder(StarScoresList &List, const Histo2d &H,
		     const double &XGuess, const double &YGuess,
		     Ellipse &Ell, int verbose)
{
  bool ok = true;
  double xBin, yBin;
  H.BinWidth(xBin,yBin);

  double xc = XGuess;
  double yc = YGuess;

  //default values  
  double wxx = 3/xBin;
  double wyy = 3/yBin;
  double wxy = 0;
#ifdef USE_FIT_TO_GUESS
  /*  This code was actually tested (using NBINS=3)
      and debugged but turns out to be useless.
      The fit with the default values seems to behave properly.
      If you think you should try to guess the span of the peak,
      you can try to use it. However, for the CFHTLS deep field
      images, there is not a signle one that has the 9 requested bins
      positive. I never found any that has 25 bins positive.
  */
#define NBINS 5
  double bins[NBINS][NBINS];
  double xx[NBINS][NBINS];
  double yy[NBINS][NBINS];
  bool allpos = true;
  for (int i=0; i <NBINS; ++i)
    for (int j=0; j<NBINS; ++j)
      {
	xx[i][j] = xc +(i-1)*xBin;
	yy[i][j] = yc +(j-1)*yBin;
	bins[i][j] = H.BinContent(xx[i][j], yy[i][j]);
	//	cout << "bins["<<i<<','<<j<<"] = " << bins[i][j] << endl;
	if (bins[i][j] <= 0) allpos = false;
      }
  
  if (allpos)
    {
      // I was not able to find an image for which we enter here
      // so this code is not tested.
      Mat A(6,6);
      Vect B(6);
      Vect h(6);
      for (int i=0; i <NBINS; ++i)
	for (int j=0; j<NBINS; ++j)
	  {
	    double x = xx[i][j];
	    double y = yy[i][j];
	    h(0) = sq(x);
	    h(1) = sq(y);
	    h(2) = x*y;
	    h(3) = x;
	    h(4) = y;
	    h(5) = 1;
	    double data = 2*log(bins[i][j]);
	    for (int ii=0; ii <6; ii++)
	      {
		for (int jj=0; jj<6; ++jj)
		  A(ii,jj) += h(ii)*h(jj);
		B(ii) += h(ii) * data;
	      }
	  }
      /*
> mf(x,y);
              2      2
           a x  + b y  + c x y + d x + e y + f

> solve({diff(mf(x,y),x),diff(mf(x,y),y)},{x,y});
         -c e + 2 b d        2 a e - d c
  {x = - ------------, y = - -----------}
          2                   2
        -c  + 4 b a         -c  + 4 b a

      */




      if (cholesky_solve(A,B)==0)
	{
	  cout << " we finally got into this never tested code...." << endl;
	  double deno = sq(B(2))-4.*B(1)*B(0);
	  xc = -(B(2)*B(4)-2*B(1)*B(3))/deno;
	  yc = -(B(3)*B(2) - 2*B(0)*B(4))/deno;
	  wxx = -B(0);
	  wyy = -B(1);
	  wxy = -0.5*B(2);
	}
    } // end if (allpos)
#endif /* USE_FIT_TO_GUESS */


  if (verbose > 0 )
    {
  cout << "HistoPeakFinder :  starting iterations with xc = " << xc 
       << " yc = " << yc 
       << " w = " << wxx << ' ' << wyy << ' ' << wxy << endl;
    }

  int count = 0;
  while (count < 20)
    {
      count ++;
#ifdef DEBUG
      cout << " iter  wxx,wyy,wxy " << count << ' ' << wxx << ' ' << wyy << ' ' << wxy << endl;
#endif
      double sumx = 0;
      double sumy = 0;
      double sumxx = 0;
      double sumyy = 0;
      double sumxy = 0;
      double sumw  = 0;

      for (StarScoresCIterator i = List.begin(); i != List.end(); ++i)
	{
	  const StarScores &star = **i;
          double x = star.sx - xc;
	  double y = star.sy - yc;
	  double arg = (wxx*x*x + wyy*y*y + 2.*wxy*x*y);
	  if (arg > 25) continue; // avoid overflows
	  double w = exp(-0.5*arg);
	  sumx += w*x;
	  sumy += w*y;
	  sumxx += w*x*x;
	  sumyy += w*y*y;
	  sumxy += w*x*y;
	  sumw += w*star.eventWeight;
	}
      if (sumw==0)
	{
	  cout << "  HistoPeakFinder no event in shape histogram : cannot figure out a peak " << endl;
	  ok = false;
	  Ell = Ellipse(XGuess, YGuess, 12.*xBin, 12.*yBin, 0.);
	  break;
	}      // normalize
      sumx /= sumw;
      sumy /= sumw;
      sumxx /= sumw;
      sumyy /= sumw;
      sumxy /= sumw;
      
      // averages and variances
      xc += sumx;
      yc += sumy;
      sumxx -= sq(sumx);
      sumyy -= sq(sumy);
      sumxy -= sumx*sumy;

      // invert to get weights
      double det = sumxx*sumyy - sq(sumxy);
      if (det <=0 || sumxx < 0 || sumyy < 0)
	{
#ifdef DEBUG
	  cout << " det, sumxx, sumyy " << det << ' ' << sumxx << ' ' << sumyy << endl;
#endif
	  cout << " HistoPeakFinder cannot figure out a peak " << endl;
	  ok = false;
	  Ell = Ellipse(XGuess, YGuess, 12.*xBin, 12.*yBin, 0.);
	  break;
	}
      // new weights
      wxx = 0.5*sumyy/det;
      wyy = 0.5*sumxx/det;
      wxy = -0.5*sumxy/det;
      // handle here the case where all values are almost equal
      double spot_size=sqrt(sqrt(det));
#ifdef DEBUG
      cout << " det, sumxx, sumyy " << det << ' ' << sumxx << ' ' << sumyy << endl;
      cout << " spot_size " << spot_size << endl;
#endif
      if (spot_size<0.0005*(xc+yc)) break;


    }

  if (ok)
    {
      Ell = Ellipse(xc,yc, wxx, wyy, wxy);
      for (StarScoresIterator i = List.begin(); i != List.end(); ++i)
	{
	  StarScores &star = **i;
	  star.nSig = Ell.Distance(star.sx,star.sy);
	}
    }
  return ok ;
}


double Ellipse::Distance(const double &X, const double &Y) const
{ 
  return sqrt(wxx*sq(X-xc)+(Y-yc)*(wyy*(Y-yc)+2.*wxy*(X-xc)));
}

double Ellipse::SigmaX() const
{
  return sqrt(wyy/( wxx*wyy-wxy*wxy));
}

double Ellipse::SigmaY() const
{
  return sqrt(wxx/( wxx*wyy-wxy*wxy));
}

double Ellipse::Corr() const
{
  return (-wxy/sqrt(wxx*wyy));
}

void Ellipse::dump(ostream &s) const
{ 
    s << " center :(" << xc << ',' << yc << ") sigma(x,y,cor) = (" 
      << SigmaX() << ','<< SigmaY()  << ',' 
      << Corr() << ')' << endl;
}
	



ostream & operator <<(ostream &s, const Ellipse &Ell) 
{
  Ell.dump(s);
  return s;
}
