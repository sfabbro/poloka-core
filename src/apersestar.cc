#include <iostream>
#include <cmath> // asin, sqrt
#include <assert.h>

#include "apersestar.h"


#define CHECK_BOUNDS /* for Image bounds */
#include "image.h"


/* concerning the choice of oversampling value, see comment at the end
   of this routine */
#define OVERSAMPLING 9


static double sq(const double &x) { return x*x;}

//! computes the aperture flux and associated scores and add one record to the apers vector.
/*! if you want several aperures, call repeteadly the same routines, with increasing radiuses (far sake of efficiency) */
void Aperture::computeflux(double x, double y, const Image& I, 
			   const Image& W, const Image *pC, const Image *pS, 
			   const double Gain, const double Radius, 
			   int segmentation_number)
{
  int imin = int(floor(x - Radius ));
  int imax = int(ceil(x + Radius ));
  int jmin = int(floor(y - Radius ));
  int jmax = int(ceil(y + Radius ));
  double rad2 = sq(Radius);
  double rmax2 = sq(Radius+0.71);
  double rint2 = sq(Radius-0.70);
  double startOffset = -0.5+0.5/double(OVERSAMPLING);
  double overstep = 1./double(OVERSAMPLING);
  double subpix = sq(overstep);

  double nnbad = 0; // double because there will be fractions of pixels...
  double frac = 0; // fraction of the current pixel within the aperture.
  double var =0;
  double fflux = 0;
  
  double nncos=0;
  double ffcos=0;
  
  int nx = I.Nx();
  int ny = I.Ny();

  double ffother = 0;

  int thisStarNumber = segmentation_number; // the S (Segmentation) image contains the "num" value of pixels that Sextractor attributed to objects

  for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax ; ++i)
      {
	double r2 = sq(i-x)+sq(j-y);
	if (r2 > rmax2) continue;
	if (r2 < rint2) frac = 1;
	else // oversample
	  {
	    double xx = i + startOffset;
	    double rx2;
	    frac = 0;
	    for (int ik = OVERSAMPLING; ik; --ik, xx+=overstep)
	      {
		rx2 = sq(x-xx);
		double yy = j + startOffset;
		for (int jk = OVERSAMPLING; jk; --jk, yy+=overstep)
		  {
		    if (sq(yy-y) + rx2 < rad2) frac += 1;
		  }
	      }
	    frac *= subpix;
	    if (fabs(frac-1)<1e-4) frac = 1;
	  }
	if (i<0 || i>= nx || j<0 || j>=ny)
	  {
	    nnbad += frac;
	    continue;
	  }
	double w = W(i,j);
	double pixVal = I(i,j);
	if (w > 0)
	  {
	    var += frac/w;
	    fflux += frac * pixVal;
	  }
	else {
	  nnbad += frac; // was += frac ???
	  if( (pC != NULL) && ((*pC)(i,j)>0) ) {
	    nncos += frac;
	    ffcos += frac * pixVal;
	  }
	}
	int seg =  0 ;
	if (pS != NULL) seg = int( (*pS)(i,j) ) ;
	if (seg !=0 && seg != thisStarNumber)
	  {
	    ffother+= frac*pixVal;
	  }
      }


  // store the result
  radius = Radius;
  flux = fflux;
  double totalVar = var+fflux/Gain;
  if (totalVar >0)
    eflux = sqrt(totalVar);
  else  eflux = 0;
  nbad = int(ceil(nnbad));
  ncos = int(ceil(nncos));
  fcos = ffcos;
  fother = ffother ;
}



void AperSEStar::ComputeFlux(const Image& I, const Image& W, const Image& C,  const Image& S, 
			     const double Gain, const double Radius, bool sort_radii)
{
  computeflux(I, W, &C, &S, Gain, Radius, sort_radii);
}

void AperSEStar::ComputeFlux(const Image& I, const Image& W, const double Gain, const double Radius, bool sort_radii)
{
  computeflux(I, W, NULL, NULL, Gain, Radius, sort_radii);
}




void AperSEStar::computeflux(const Image& I, const Image& W, const Image *pC, 
			     const Image *pS, 
			     const double Gain, const double Radius, 
			     bool sort_radii)
{
  apers.push_back(Aperture());
  Aperture &aper = apers.back();
  aper.computeflux(x,y,I,W,pC,pS,Gain,Radius,this->num);
  // impose increasing order
  unsigned naper = apers.size();
  if (sort_radii && (naper>1 || Radius < apers[naper-1].radius))
    sort(apers.begin(),apers.end());
}

/* When computing aperture fluxes, one may consider, for pixels on
the aperture border, the fraction of their area that lies within the aperture.

  If one choses oversampling (cut the pixel into smaller squares and count
the fraction that have their center inside the aperture), 
here are a few figures that justify the choice of oversampling by a factor 
of 9 (computed on a sample of 1000):

  aperture(pix)   <area/pi*r^2-1>    <r.m.s(area/pi*r^2-1)>
                over=5   over=9         over=5     over=9
        3       -0.7e-5   0.7e-4       3.7e-3     1.5e-3
        5       0.15e-4   0.18e-4      1.8e-3     0.6e-3
        7      -0.2e-5    0.77e-5      1.0e-3     0.4e-3
        9       0.3e-4    0.7e-5       0.6e-3     0.25e-3
       11       0.13e-5   0.16e-4      0.5e-3     0.22e-3

The r.m.s's of the computed area are of course upper bounds to the contribution
to a flux measurement error introduced by the computation, because
this error is introduced at the border, and most of the object fluxes
comes from the center, where pixels don't have to be cut into sub-pieces.

going to an oversampling of 9 increases the CPU by 50%


*/

//! interpolates existing measurements (avoids going back to pixels)
/*! assumes that apertures are in increasing order, which is
  enforced by ComputeFlux if called with default parameters sort_radii set to true*/
Aperture AperSEStar::InterpolateFlux(const double &Radius) const
{
  unsigned naper = apers.size();
  unsigned k=0;
  for (   ; k < naper; ++k)
    if (apers[k].radius > Radius) break;
  if (k == 0) return apers[0];
  if (k == naper) return apers[k-1];
  k--;
  double x = (Radius-apers[k-1].radius)/(apers[k].radius-apers[k-1].radius);
  Aperture aper;
  aper.radius = Radius;
  aper.flux = (1-x)*apers[k-1].flux+x*apers[k].flux;
  // interpolate errors linearly because for low fluxes, it scales with radius
  aper.eflux = (1-x)*apers[k-1].eflux+x*apers[k-1].eflux;
  aper.nbad = apers[k].nbad;
  return aper;
}



// recompute the star position. 
// generally used on strongly defocussed stars.
// ==> we just compute the first order moment
// 
// It may be possible to also compute the 2d order moments
// and weight them with a defocussed PSF model... 

// should be the same as ComputePos, except that 
// we do not weight w/ a gaussian PSF...
void AperSEStar::ComputePos(const Image& I, const Image& W)
{
  int i,j;
  double radius = 20;
  
  int nbad = 0;
  int iter = 0;
  
  int nx = I.Nx(), ny = I.Ny();
  double wX = x, wY = y;
  double sumx=0., sumy=0., sumw=0.;
  
  int imin = int(floor(x-radius));
  int imax = int(ceil(x+radius));
  int jmin = int(floor(y-radius));
  int jmax = int(ceil(y+radius));
  
  cout << "[ComputePos] Initial position: x=" << x << " y=" << y << endl;
  while(iter < 10) 
    {
      iter++;
      
      for(j=jmin;j<=jmax;j++)
	{
	  if(j<0 || j>=ny) continue;
	  for(i=imin;i<imax;i++) 
	    {
	      if(i<0 || i>=nx) continue;
	      double dx = i-wX;
	      double dy = j-wY;
	      double w = W(i,j);
	      if(w<0) { nbad+=1; continue; }
	      double flx = I(i,j);
	      sumx += flx*dx;
	      sumy += flx*dy;
	      sumw += flx;
	    }
	}
      sumx /= sumw;
      sumy /= sumw;
      cout << "      dx=" << sumx
	   << " dy=" << sumy << endl;
      wX += sumx;
      wY += sumy;
    }
  
  x = wX;
  y = wY;
  cout << "[ComputePos] Final position: x=" << x << " y=" << y << endl;
}

//#define DEBUG

void AperSEStar::ComputeShapeAndPos(const Image&I, const Image &W, 
				    const double &Gain)
{
  gflag &= ~(BAD_GAUSS_MOMENTS);
  double det = Mxx()*Myy() - sq(Mxy());
  double wxx, wyy, wxy, radius;
  if (Mxx() >0 && Myy()>0 && det >0)
    {
      wxx = Myy()/det;
      wyy = Mxx()/det;
      wxy = - Mxy()/det;
      /* 0.5*(mxx+myy) + sqrt(sq(0.5*(mxx-myy))+sq(mxy)) is the largest
	 eigenvalue of the second moment matrix. This formula just
	 says that we integrate up to 4 sigma, and at least up to 3
	 pixels.  this formula is copied later in the same routine.
      */
      double largest_eigenval = 0.5*(Mxx()+Myy())+sqrt(sq(0.5*(Mxx()-Myy()))+sq(Mxy()));
      radius = max(4*sqrt(largest_eigenval),3.);
      radius = min(radius,50.);    
    }
  else 
    {
      wxx = wyy = 0.5; wxy = 0;
      radius = 10.;
    }

  int iter=0;
  int nbad;
  int nx = I.Nx();
  int ny = I.Ny();
  double weightedX = x;
  double weightedY = y;
  double varxx = -1;
  double varyy = -1;
  double varxy = 0;

  while (iter < 50)
    {
      iter++;
      int imin = int(floor(x - radius ));
      int imax = int(ceil(x + radius ));
      int jmin = int(floor(y - radius ));
      int jmax = int(ceil(y + radius ));

      // position
      double sumx = 0;
      double sumy = 0;
      // weighted secon moments
      double sumxx = 0;
      double sumyy = 0;
      double sumxy = 0;
      double sumw = 0;
      // for position variance
      double sumxxw2 = 0;
      double sumyyw2 = 0;
      double sumxyw2 = 0;



      nbad = 0;
      for (int j = jmin; j <= jmax; ++j)
	for (int i = imin; i <= imax ; ++i)
	  {
	    double dx = i-weightedX;
	    double dy = j-weightedY;
	    double wg = wxx*dx*dx + wyy*dy*dy + 2.*wxy*dx*dy;
	    if (wg > 16) continue; // 4 sigmas, and avoids overflows
	    // test now if this useful pixel lies within image bounds
	    if (i<0 || i>= nx || j<0 || j>=ny)
	      {
		nbad += 1;
		continue;
	      }
            // Gaussian weighting "window function":
	    wg = exp(-0.5*wg);
	    double pix = I(i,j);
	    double wgi = wg*pix;
	    sumx += wgi*dx;
	    sumy += wgi*dy;
	    sumxx += wgi*dx*dx;
	    sumyy += wgi*dy*dy;
	    sumxy += wgi*dx*dy;
	    sumw += wgi;
	    double w = W(i,j);
	    if (w <= 0) nbad += 1;
	    else 
	      {
	    // ingredients for the pixel variance
	    /* contrarily to the above computations, this is only useful
	       at the last iteration */
		double pixvar = 1/w + pix/Gain;
		if (pixvar<=0) continue;
		double wv = wg*wg*pixvar;
		sumxxw2 += dx*dx*wv;
		sumyyw2 += dy*dy*wv;
		sumxyw2 += dx*dy*wv;		
	      }
	    
	  } // end loop on pixels
      if (sumw <= 0)
	{
	  gflag |= BAD_GAUSS_MOMENTS;
	  wxx = wyy = 12;
	  wxy = 0;
#ifdef DEBUG
	  cout << " sumw <= 0 " << endl;
#endif
	  break; // end of iterations
	}
      sumx /= sumw;
      sumy /= sumw;
      sumxx /= sumw;
      sumyy /= sumw;
      sumxy /= sumw;

      weightedX += sumx;
      weightedY += sumy;
      sumxx -= sq(sumx);
      sumyy -= sq(sumy);
      sumxy -= sumx*sumy;

      

      double det = sumxx*sumyy - sq(sumxy);
      // handle degenerate or almost degenerate matrices
      // 1/0.0833 = 12.
#define MIN_VAR 0.0833333
#define MIN_DET_VAR (MIN_VAR*MIN_VAR)
      if (det < MIN_DET_VAR)
	{
	  sumxx += MIN_VAR;
	  sumyy += MIN_VAR;
	  det = sumxx*sumyy - sq(sumxy);
	}
      /* the determinant can still be negative here since 
	 nothing forces I*wg to be positive. On low S/N detections
	 it happens that det < 0, even after the single pixel trick
	 just above.
	 An other reason to give up is that the position drifts
	 w.r.t the unweighed scheme. We then keep the original 
	 SExtractor position.
      */
      double drift2 = sq(weightedX-x)+sq(weightedY-y);
      if ( drift2>4 || det <= 0 || sumxx <= 0 || sumyy <= 0 )  
	{
	  gflag |= BAD_GAUSS_MOMENTS;
	  wxx = wyy = 12.;
	  wxy = 0.;
#ifdef DEBUG
	  cout << " drift2>4 || det < 0 || sumxx < 0 || sumyy < 0 " << endl;
#endif
	  break; 
	}
      /* 0.5*(mxx+myy) + sqrt(sq(0.5*(mxx-myy))+sq(mxy)) is the largest
	 eigenvalue of the second moment matrix. This formula just
	 says that we integrate up to 4 sigma, and at least up to 3
	 pixels.  This formula is copied from above */ 
      //      double oldradius = radius;
      double largest_eigenval = 0.5*(sumxx+sumyy)+sqrt(sq(0.5*(sumxx-sumyy))+sq(sumxy));
      radius = max(4*sqrt(largest_eigenval),3.);
      radius = min(radius,50.);

      double new_wxx = 0.5*sumyy/det;
      double new_wyy = 0.5*sumxx/det;
      double new_wxy = -0.5*sumxy/det;
      double diff_weights = fabs(wxx-new_wxx)+fabs(wyy-new_wyy)+fabs(wxy-new_wxy);
      double eps = 1e-10;
      bool stop_iter = (diff_weights<eps);

#ifdef DEBUG
      cout << " wxx wyy wxy " << wxx << ' ' << wyy << ' ' << wxy << endl;
      cout << " delta " << new_wxx-wxx << ' ' << new_wyy-wyy << ' ' << new_wxy-wxy << endl;
#endif
      // convergence of the iterations : speed up if the thing is almost stabilized.
      // we use the radius to decide upon that
      double alpha = 1;
      if (diff_weights<0.1) alpha = 1.8;
      wxx = alpha*new_wxx + (1.-alpha)*wxx;
      wyy = alpha*new_wyy + (1.-alpha)*wyy;
      wxy = alpha*new_wxy + (1.-alpha)*wxy;
      if (wxx<0 || wyy <0 || wxx*wyy<sq(wxy))
	{
	  wxx = new_wxx;
	  wyy = new_wyy;
	  wxy = new_wxy;
	}
      if (stop_iter) 
	{
	  // position variance is not needed to loop
	  varxx = sumxxw2/sq(sumw);
	  varyy = sumyyw2/sq(sumw);
	  varxy = sumxyw2/sq(sumw);
	  if (!(varxx>0 && varyy>0 && sq(varxy)<varxx*varyy))
	    gflag |= BAD_GAUSS_POS_VARIANCE;
	  break;
	}
    }// end of iterations.
#ifdef DEBUG
  cout << " niter = " << iter << endl;
#endif
  det = wxx*wyy-sq(wxy);
  // fill AperSestar scores :
  gmxx = wyy/det;
  gmyy = wxx/det;
  gmxy = -wxy/det;
  if (gmxx<0 || gmyy<0) // should never happen
    {
      cout << " something went wrong in ComputeShapeAndPos : negative second moment" << endl;
      abort();
      // If you feel like removing this test, please consider finding
      // what causes negative second moments.

  }
  if (nbad != 0) gflag |= (BAD_GAUSS_MOMENTS+BAD_GAUSS_POS_VARIANCE);
  if (varxx == -1) gflag |= BAD_GAUSS_POS_VARIANCE;
  if ((gflag & (BAD_GAUSS_MOMENTS+BAD_GAUSS_POS_VARIANCE)) == 0)
    { // can then use the "fitted" position, and its errors
      x = weightedX;
      y = weightedY;
      vx = varxx;
      vy = varyy;
      vxy = varxy;
      assert(vx>0 && vy>0 && sq(vxy)<vx*vy);
    }
}
  


ostream& operator << (ostream &stream, const Aperture &A)
{
  stream << " rad= " << A.radius << " f=" << A.flux 
	 << " ef=" << A.eflux << " nbad=" << A.nbad;
  return stream;
} 


void AperSEStar::SetNeighborScores(const BaseStar &Neighbor)
{
  neighborDist = Neighbor.Distance(*this);
  neighborFlux = Neighbor.flux;
}


std::string AperSEStar::WriteHeader_(ostream & pr, 
				     const char* i) const
{
  // write "global" values:
  std::string sestarFormat = SEStar::WriteHeader_(pr, i);
  if (!i) i="";
  pr << "#neid"<<i<< " : distance to nearest neighbor" << endl;
  pr << "#neifl"<<i<< " : flux of nearest neighbor" << endl;
  pr << "#gmxx"<<i<< " : gaussian filtered 2nd moment (xx)" << endl;
  pr << "#gmyy"<<i<< " : gaussian filtered 2nd moment (yy)" << endl;
  pr << "#gmxy"<<i<< " : gaussian filtered 2nd moment (xy)" << endl;
  pr << "#gflag"<<i<< " : flag for computation of gaussian filtered 2nd moments" << endl;
  pr << "#naper"<<i<< " : number of apertures" << endl;
  for (unsigned k=0; k < apers.size(); ++k)
    {
      char kk[8];
      sprintf(kk,"%-d",k);
      pr << "#rad"  << kk << i << " : aperture radius"<< endl;
      pr << "#apfl" << kk << i << " : aperture flux"<< endl;
      pr << "#eapfl"<< kk << i << " : error on aperture flux"<< endl;
      pr << "#apnb" << kk << i << " : nb of bad pixels in aper (including cosmics)"<< endl;
      pr << "#apnc" << kk << i << " : nb of pixels flagged as cosmics" << endl;
      pr << "#apfc" << kk << i << " : flux of the pixels above" << endl;
      pr << "#apfo" << kk << i << " : flux of the pixels attributed to other stars" << endl;
    }
  return sestarFormat+ "AperSEStar 3"; 
}



void AperSEStar::writen(ostream& s) const
{
  SEStar::writen(s);
  s << neighborDist << ' '
    << neighborFlux << ' '
    << gmxx << ' '
    << gmyy << ' '
    << gmxy << ' '
    << gflag << ' '
    << apers.size() << ' ';
  for (unsigned k=0; k < apers.size(); ++k)
    {
      const Aperture &aper = apers[k];
      s <<  aper.radius << ' ';
      s << aper.flux << ' ';
      s << aper.eflux << ' ';
      s << aper.nbad << ' ';
      s << aper.ncos << ' ';
      s << aper.fcos << ' ';
      s << aper.fother << ' ';
    }
}


#include "fastifstream.h"

void AperSEStar::read_it(fastifstream& r, const char* Format)
{
  SEStar::read_it(r, Format);
  int format = DecodeFormat(Format, "AperSEStar");
  if(format != 2 && format != 3) {
    cerr << "AperSEstar::read_it() format: " << format << ". The current version is 2 or 3!" << endl;
    abort();
  }
  r >> neighborDist 
    >> neighborFlux;

  r >> gmxx >> gmyy >> gmxy >> gflag ; 

  unsigned naper;
  r >> naper;
  apers.reserve( naper );
  for (unsigned k=0; k < naper; ++k)
    {
      apers.push_back(Aperture());
      Aperture &aper = apers.back();
      r >> aper.radius;
      r >> aper.flux;
      r >> aper.eflux;
      r >> aper.nbad;
      r >> aper.ncos;
      r >> aper.fcos;
      if (format >=3) r >> aper.fother;
    }
}

AperSEStar* AperSEStar::read(fastifstream& r, const char* Format)
{
  AperSEStar *pstar = new AperSEStar();  
  pstar->read_it(r, Format);
  return(pstar);
}
  


#include "starlist.cc"

//instantiate all template routines from StarList
template class StarList<AperSEStar>;

// just copy and initialize
AperSEStarList::AperSEStarList(const SEStarList &L)
{
  for (SEStarCIterator i = L.begin(); i != L.end(); ++i)
    push_back(new AperSEStar(**i));
}


#include "fastfinder.h"

void AperSEStarList::SetNeighborScores(const BaseStarList &List, 
				       const double Maxdist)
{
  FastFinder finder(List);
  for (AperSEStarIterator i= begin(); i != end(); ++i)
    {
      AperSEStar &current = **i;
      const BaseStar *closest = NULL;
      // as specified, 'List' contains current objects, and we hence
      // look for second closest.

      const BaseStar *neighbor = finder.SecondClosest(current, Maxdist, closest);
      if (!neighbor) continue;
      /* if everything goes well we have a closest because we have a
	 second closest (neighbor). Hence no test on closest!=NULL */
      if (current.Distance(*closest) > 2)
	{
	  cout << "SetNeighborScores: missing object in neighbor list? object" << endl
	       << BaseStar(current) 
	       << " first closest is : " << endl
	       << *closest 
	       << " dist = " << current.Distance(*closest) << endl;
	}
      current.SetNeighborScores(*neighbor);
    }
}

int AperSEStarList::write(const std::string &FileName) const
{
  if (!empty())
    {
      unsigned naper = front()->apers.size();
      for (AperSEStarCIterator i = begin(); i !=end(); ++i)
	{
	  if ((*i)->apers.size() != naper)
	    {
	      cerr << " writing a AperSEStarList  with different number of apertures" << endl;
	      cerr << " hope you know what you are doing " << endl;
	    }
	}
    }
  
  return StarList<AperSEStar>::write(FileName);
}

/*************************************/

#include "histopeakfinder.h"
#include "histo2d.h"

bool FindStarShapes(const AperSEStarList &List, const double MinSN,
		    double &XSize, double &YSize, double &Corr, 
		    StarScoresList &scores)
{
  Histo2d histo(30,0,10.,30,0.,10.);
  XSize = 0;
  YSize = 0;
  Corr = 0;
  for (AperSEStarCIterator i = List.begin(); i != List.end(); ++i)
    {
      const AperSEStar &s = **i;
      /* basic quality cut : gflag = 0 means 2nd moments 
	 (gmxx & co) computed sucessfully. Flag is the SExtractor 
	 extraction flags, and Flagbad refers to bad pixels.
      */
      if (s.gflag || s.Flag() || s.FlagBad()) continue;
      /* if the s/n is poor, the measurement of 2nd moments
	 is also poor: cut on S/N before entering into the
	 "star clump finder".
      */
      if (s.flux<0 || s.EFlux() < 0 || s.flux< MinSN * s.EFlux()) continue;
      double xx = sqrt(s.gmxx);
      double yy = sqrt(s.gmyy);
      histo.Fill(xx,yy, s.Fluxmax());
      scores.push_back(new StarScores (xx,yy,&s));
    }

  unsigned nobj = scores.size();
  if (nobj<5)
    {
      cout << " FindStarShapes : only " << nobj << " objects passing the cuts " << endl;
      return false;
    }

  Ellipse ellipse;
  bool ok;
  for (unsigned k=0;k<3; ++k)
    {
      if (k>=1) cout << " INFO : FindStarShape trying histo max # " << k << endl;
      double xMax, yMax;
      histo.MaxBin(xMax, yMax);
      ok = HistoPeakFinder(scores, histo, xMax, yMax, ellipse);
      if (!ok)
	{
	  cout << " FindStarShape could not parametrize a peak using HistoPeakFinder with histo max" << endl;
      /* we could do that, if it is useful, but we then have to remove the 
	 XSize = 0 upstream... */
      //      if (XSize != 0) ok =  HistoPeakFinder(scores, histo, XSize, YSize, ellipse);
	}
      else break;
      histo.ZeroBin(xMax,yMax);
    }
  cout << " shape ellipse " << ellipse << endl;
  ellipse.GetCenter(XSize, YSize);
  if (!ok)
    {
      cout << " could not parametrize a peak using HistoPeakFinder" << endl
	   << " (very) crude shape Ellipse " << endl;
      Corr = 0;
      return false;
    }
  else 
    {
      //compute average correlation coefficient
      int count = 0;
      double sum = 0; 
      for (StarScoresIterator i = scores.begin(); i != scores.end(); ++i)
	{
	  StarScores &star = **i;
	  if (star.nSig > 4) continue;
	  const AperSEStar & a = dynamic_cast<const AperSEStar&>(*star.star);
	  sum += a.gmxy;
	  count ++;
	}
      Corr = sum/(count*XSize*YSize);
    }
  return true;
}




/************ converters *************/
BaseStarList* AperSE2Base(AperSEStarList * This)
{ return (BaseStarList *) This;}

const BaseStarList* AperSE2Base(const AperSEStarList * This)
{ return (const BaseStarList *) This;}

BaseStarList& AperSE2Base(AperSEStarList &This)
{ return (BaseStarList &) This;}

const BaseStarList& AperSE2Base(const AperSEStarList &This)
{ return (const BaseStarList &) This;}

/* ================================== */

SEStarList* AperSE2SE(AperSEStarList * This)
{ return (SEStarList *) This;}

const SEStarList* AperSE2SE(const AperSEStarList * This)
{ return (const SEStarList *) This;}

SEStarList& AperSE2SE(AperSEStarList &This)
{ return (SEStarList &) This;}

const SEStarList& AperSE2SE(const AperSEStarList &This)
{ return (const SEStarList &) This;}


