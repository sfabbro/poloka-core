#include <iostream>
#include <cmath> // asin, sqrt


#include "apersestar.h"


#define CHECK_BOUNDS
#include "image.h"




// premier probleme, quelle est la fraction de pixel couverte par un disque de rayon R ?
// x_pix, y_pix coord du centre du pixel (de taille 1x1 par convention)

static double sq(const double x) {return x*x;}

static double fraction(double x_centre, double y_centre, double R, double x_pix, double y_pix)
{
  double x_min, x_max;
	// Je m'arrange pour que le pixel soit dans le premier quadrant par rapport au centre du cercle
	
	if (x_pix < x_centre) {x_pix = 2.*x_centre - x_pix;};
	if (y_pix < y_centre) {y_pix = 2.*y_centre - y_pix;};
       
	// Je m'arrange même pour qu'il soit dans la zone y>x
	
	if ((x_pix-x_centre)>(y_pix-y_centre))
	{
		double prov = x_pix; x_pix = y_pix; y_pix = prov;
		prov = x_centre; x_centre = y_centre; y_centre = prov;
	}
	
#ifdef DEBUG
	cout << x_centre << " " << y_centre << " " << R << " " << x_pix << " " << y_pix << endl;
#endif

	// je definis les coordonnees par rapport au coin inf gauche du pixel

	x_centre -= x_pix - 0.5;
	y_centre -= y_pix - 0.5;

#ifdef DEBUG
	cout << "coord centre utilisees : " << x_centre << " " << y_centre << endl;
#endif
	// si la distance entre centre et bord inf gauch du pix est > R, le pixel est dehors (faux pour R petit)

	if (sq(x_centre)+sq(y_centre) > pow(R,2.)) return 0.;

	// si la distance entre centre et bord sup droit du pix est < R, le pixel est dedans (faux pour R petit)

	if (pow(x_centre-1.,2.)+pow(y_centre-1.,2.) < R*R) return 1.;

	// L'equation du cercle est y = y_centre + sqrt(R^2-(x-x_centre)^2)
	
	// 
	// On calcule les intersections de ce cercle avec les bords du pixel.
	//
	
		
	// en quels points du bord inferieur du pixel passe le cercle.

	x_min = x_centre - sqrt(R*R - sq(y_centre));
	x_max = x_centre + sqrt(R*R - sq(y_centre));

	if (x_min > x_max) { double prov = x_min; x_min=x_max; x_max=prov;}
	
#ifdef DEBUG
	cout << "intersections bord inf : " << x_min << " " << x_max << endl;
#endif

	if (x_min < 0.) x_min = 0.;
	if (x_max > 1.) x_max = 1.;


	// Si le point UL du pixel n'est pas dans le cercle, il suffit de calculer l'integrale

	double complement;
	if (sq(x_centre)+sq(y_centre+1.) < R*R)
	  {
#ifdef DEBUG
	    cout << "cas 1 : complement = 0" << endl;
#endif
	    complement = 0.;
	  }
	else
	// Sinon il faut recalculer x_min, le point ou le cercle coupe le bord sup du pixel
	  {
	    double x_min1 = x_centre - sqrt(R*R - sq(y_centre-1.));
	    double x_min2 = x_centre + sqrt(R*R - sq(y_centre-1.));

#ifdef DEBUG
	    cout << "intersections bord sup : " << x_min1 << " " << x_min2 << endl;
#endif

	    x_min = x_min1;
	    if (x_min<x_min2) x_min=x_min2;
	    complement = x_min;
#ifdef DEBUG
	    cout << "cas 2 :complement = " << x_min << endl;
#endif
	  }

	//int(sqrt(R^2-x^2),x) = 1/2 x sqrt(R^2-x^2) + 1/2 R^2 arcsin (x/R)

	double integrale = (x_max-x_centre) * sqrt(R*R - sq(x_max-x_centre))/2. 
	  + R*R* asin ((x_max-x_centre)/R)/2. + y_centre*(x_max-x_centre);
	integrale -= (x_min-x_centre) * sqrt(R*R - sq(x_min-x_centre))/2. 
	  + R*R* asin ((x_min-x_centre)/R)/2. + y_centre*(x_min-x_centre);
	integrale += complement;
	  
#ifdef DEBUG
	cout << integrale << endl;
#endif
	
	return x_min + integrale;	
}




/* concerning the choice of oversampling value, see comment at the end
   of this routine */
#define OVERSAMPLING 9


//! computes the aperture flux and associated scores and add one record to the apers vector.
/*! if you want several aperures, call repeteadly the same routines, with increasing radiuses (far sake of efficiency) */
void AperSEStar::ComputeFlux(const Image& I, const Image& W, const Image& C, 
			     const double Gain, const double Radius)
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

  double nbad = 0; // double because there will be fractions of pixels...
  double frac = 0; // fraction of the current pixel within the aperture.
  double var =0;
  double flux = 0;
  
  double ncos=0;
  double fcos=0;
  
  int nx = I.Nx();
  int ny = I.Ny();

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
	    nbad += frac;
	    continue;
	  }
	double w = W(i,j);
	if (w > 0)
	  {
	    var += frac/w;
	    flux += frac * I(i,j);
	  }
	else {
	  nbad += frac; // was += frac ???
	  if( C(i,j)>0 ) {
	    ncos += frac;
	    fcos += frac * I(i,j);
	  }
	}
      }

  // store the results
  apers.push_back(Aperture());
  Aperture &aper = apers.back();
  // copy
  aper.radius = Radius;
  aper.flux = flux;
  double totalFlux = var+flux/Gain;
  if (totalFlux >0)
    aper.eflux = sqrt(var+flux/Gain);
  else  aper.eflux = 0;
  aper.nbad = int(ceil(nbad));
  aper.ncos = int(ceil(ncos));
  aper.fcos = fcos;

  // impose increasing order
  unsigned naper = apers.size();
  if (naper>1 || Radius < apers[naper-1].radius)
    sort(apers.begin(),apers.end());
}

/* When computing aperture fluxes, one may consider, for pixels on
the aperture border, the fraction of their area that lies within the aperture.
t
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
  enforced by ComputeFlux*/
AperSEStar::Aperture AperSEStar::InterpolateFlux(const double &Radius) const
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




void AperSEStar::ComputeShapeAndPos(const Image&I, const Image &W)
{
  gflag &= ~(BAD_GAUSS_MOMENTS);
  double det = Mxx()*Myy() - sq(Mxy());
  double wxx, wyy, wxy, radius;
  if (det!=0)
    {
      wxx = Myy()/det;
      wyy = Mxx()/det;
      wxy = - Mxy()/det;
      double half_trace = 0.5*(Mxx()+Myy());
      /* half_trace + sqrt(sq(half_trace)-det) is the largest
	 eigenvalue of the second moment matrix. This formula just
	 says that we integrate up to 4 sigma, and at least up to 3
	 pixels.  this formula is copied later in the same routine */
      radius = max(4*sqrt(half_trace + sqrt(sq(half_trace) - det)),3.);
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

  while (iter < 50)
    {
      iter++;
      int imin = int(floor(x - radius ));
      int imax = int(ceil(x + radius ));
      int jmin = int(floor(y - radius ));
      int jmax = int(ceil(y + radius ));

      double sumx = 0;
      double sumy = 0;
      double sumxx = 0;
      double sumyy = 0;
      double sumxy = 0;
      double sumw = 0;

      nbad = 0;
      for (int j = jmin; j <= jmax; ++j)
	for (int i = imin; i <= imax ; ++i)
	  {
	    double dx = i-weightedX;
	    double dy = j-weightedY;
	    double wg = wxx*dx*dx + wyy*dy*dy + 2.*wxy*dx*dy;
	    if (wg > 16) continue; // 4 sigmas, and avoids overflows
	    wg = exp(-0.5*wg);
	    if (i<0 || i>= nx || j<0 || j>=ny)
	      {
		nbad += 1;
		continue;
	      }
	    double w = W(i,j);
	    if (w <= 0)
	      {
		nbad += 1;
	      }
	    wg *= I(i,j);
	    sumx += wg*dx;
	    sumy += wg*dy;
	    sumxx += wg*dx*dx;
	    sumyy += wg*dy*dy;
	    sumxy += wg*dx*dy;
	    sumw += wg;
	  }
      if (sumw <= 0)
	{
	  gflag |= BAD_GAUSS_MOMENTS;
	  wxx = wyy = 12;
	  wxy = 0;
	  break; 
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
      if ( drift2>4 || det < 0 || sumxx < 0 || sumyy < 0 )  
	{
	  gflag |= BAD_GAUSS_MOMENTS;
	  wxx = wyy = 12;
	  wxy = 0;
	  break; 
	}
      double half_trace = 0.5*(sumxx + sumyy);
      /* half_trace + sqrt(sq(half_trace)-det) is the largest
	 eigenvalue of the second moment matrix. This formula just
	 says that we integrate up to 4 sigma, and at least up to 3
	 pixels. This formula is copied from above */
      radius = max(4*sqrt(half_trace + sqrt(sq(half_trace) - det)),3.);
      radius = min(radius,50.);

      wxx = 0.5*sumyy/det;
      wyy = 0.5*sumxx/det;
      wxy = -0.5*sumxy/det;
    }// end of iteration.

  
  det = wxx*wyy-sq(wxy);
  // fill AperSestar scores :
  gmxx = wyy/det;
  gmyy = wxx/det;
  gmxy = -wxy/det;
  if (gmxx <0) abort(); // hopefully never happens
  if (nbad != 0) gflag |= BAD_GAUSS_MOMENTS;
  if ((gflag & BAD_GAUSS_MOMENTS) == 0)
    { // can then use the "fitted" position
      x = weightedX;
      y = weightedY;
    }
}
  


ostream& operator << (ostream &stream, const AperSEStar::Aperture &A)
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
    }
  return sestarFormat+ "AperSEStar 2"; 
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
    }
}


void AperSEStar::read_it(istream& r, const char* Format)
{
  SEStar::read_it(r, Format);
  int format = DecodeFormat(Format, "AperSEStar");
  if(format != 2) {
    cerr << "AperSEstar::read_it() version mismatch. The current version is 2!" << endl;
    abort();
  }
  r >> neighborDist 
    >> neighborFlux;

  r >> gmxx >> gmyy >> gmxy >> gflag ; 

  unsigned naper;
  r >> naper;
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
    }
}

AperSEStar* AperSEStar::read(istream& r, const char* Format)
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

bool FindStarShapes(const AperSEStarList &List, double &XSize, 
		    double &YSize, double &Corr)
{
  Histo2d histo(30,0,10.,30,0.,10.);
  StarScoresList scores;
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
      if (s.flux<0 || s.EFlux() < 0 || s.flux< 10 * s.EFlux()) continue;
      double xx = sqrt(s.gmxx);
      double yy = sqrt(s.gmyy);
      histo.Fill(xx,yy, s.Fluxmax());
      scores.push_back(new StarScores (xx,yy,&s));
    }

  double xMax, yMax;
  histo.MaxBin(xMax, yMax);
  Ellipse ellipse;
  bool ok = HistoPeakFinder(scores, histo, xMax, yMax, ellipse);
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


