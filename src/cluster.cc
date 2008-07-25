#include "cluster.h"
#include "reducedimage.h"
#include "fitsimage.h"

#ifndef MAXCLUSTERSIZE
/* was chosen  at 80000, I changed it to this large value to 
   accomodate saturated bright stars */
#define MAXCLUSTERSIZE  200000
#endif

#ifndef MY_PI
#define MY_PI  3.14159265358979
#endif

double sqr(const double & x){return x*x;};


Cluster::Cluster(long c) 
{
  color = c ;
  size = 0 ;
  xsum = ysum = x2sum = xysum = y2sum = 0.0 ;
}



long Cluster::browse_and_color(const Image& src, long x0, long y0, 
		      double threshold, Image& map) 
{
  long xsize = src.Nx() ; long ysize = src.Ny() ;
  
  //      cerr << "debug: killsatellites: browse_and_color: begin "
  //  	 << "(color=" << color << ")" 
  //  	 << "((x0 y0)=(" << x0 << " " << y0 << ")" << ")" 
    //  	 << endl ;

    if (size < 0) {
      return -1 ;
    }

    if ((x0 < 0) || (x0 >= xsize) ||
	(y0 < 0) || (y0 >= ysize)) {
      return 0 ;
    }
    
    if (map(x0, y0) != 0) {
      return 0 ;
    }
    
    if (src(x0, y0) < threshold) {
      return 0 ;
    }
    
    map(x0, y0) = color ;
    
    size += 1 ;
    xsum += x0 ;
    ysum += y0 ;
    x2sum += x0*x0 ;
    xysum += x0*y0 ;
    y2sum += y0*y0 ;
    
    long r ;

    if (size > MAXCLUSTERSIZE) {
      // warning: killsatellites: too many recursive calls.
      r = -1 ;
    } else {
      long north, south, east, west ;
      north = browse_and_color(src, x0, y0+1, threshold, map) ;
      south = browse_and_color(src, x0, y0-1, threshold, map) ;
      east = browse_and_color(src, x0+1, y0, threshold, map) ;
      west = browse_and_color(src, x0-1, y0, threshold, map) ;
      if ((north < 0) || (south < 0) || (east < 0) || (west < 0)) {
	r = -1 ;
      } else {
	r = 1 + north + south + east + west ;
      }
    }
    
    //      cerr << "debug: killsatellites: browse_and_color: end "
    //  	 << "(color=" << color << ")" 
    //  	 << "(r=" << r << ")" 
    //  	 << endl ;
    
    return r ;
}


ClusterList::ClusterList(ReducedImage & inrim, int debuglev)
{
  debuglevel = 1;//debuglev;
  cout << "#### " << inrim.FitsName() << endl;
  FitsImage infits(inrim.FitsName());
  {
  FitsImage weight(inrim.FitsWeightName());
  FitsImage satur(inrim.FitsSaturName());
  for(long y=0 ; y<infits.Ny() ; y++) 
    for(long x=0 ; x<infits.Nx() ; x++) 
      if (weight(x,y)== 0  || satur(x,y) ==1) 
	{
	  infits(x,y) = 0.0 ;
	}
  }
  double nsigma = 1 ;
  double threshold =0;
  
  const int RELATIVE = 2 ;
  int giventhreshold = RELATIVE ;


  long xsize = infits.Nx() ;
  long ysize = infits.Ny() ;
  //
  colors = new Image(xsize, ysize) ;


  //double ccddiagonal = sqrt(double(xsize * xsize + ysize * ysize)) ;
  int bigcluster = 1;
  
  // we should give sigmaback and backlevel from reducedimage
  double fondciel = inrim.BackLevel();
  double sigfond = inrim.SigmaBack();
  threshold = nsigma * sigfond;

  cout << fondciel << " " << sigfond << endl;
  while (bigcluster) {
    
    bigcluster = 0 ;
    
    if (debuglevel >= 1) 
      { 
	cerr << "debug: killsatellites: threshold type = "
	     << (giventhreshold == RELATIVE ? "relative" : "absolute")
	     << endl
	     << "debug: killsatellites: threshold = " << threshold
	     << endl ;
      }
    
    
    // -- Looking for clusters
    
    (*colors) = 0.0 ;
    
    // Main Loop on pixels
    
    long  x0, y0 ;
    long nextcolor = 10 ;
    long npixels ;
    
    for(y0=0 ; (y0<ysize) && (bigcluster == 0) ; y0++) {
      for(x0=0 ; (x0<xsize) && (bigcluster == 0) ; x0++) {
	
	float val = infits(x0,y0) ;
	
	// 0. is this point above threshold ?
	
	if (val < threshold) 
	  continue ;
	
	// 1. is this point already part of a known cluster ?
	//    (i.e. is it already colored ?)
	
	if ((*colors)(x0,y0) != 0) 
	  continue ;
	
	Cluster current(nextcolor) ;
	
	npixels = current.browse_and_color(infits, x0, y0, 
					    threshold, (*colors)) ;
	
	if (npixels < 0) { // Error
	  if (giventhreshold == RELATIVE) 
	    {
	      cerr << "warning: killsatellites: too many recursive calls." 
		   << endl
		   << "warning: killsatellites: there is an enormous cluster..." 
		   << endl
		   << "warning: killsatellites: perhaps threshold is too low." 
		   << endl
		   << "warning: killsatellites: we will increase nsigma by 0.1"
		   << endl ;
	      nsigma += 0.1 ;
	      threshold = fondciel + nsigma * sigfond ;
	      bigcluster = 1 ;
	      // since we restart from the beginning (fix by P.A on july-08)
	      clear();
	    } 
	  else 
	    {
	      cerr << "warning: killsatellites: too many recursive calls." 
		   << endl
		   << "warning: killsatellites: there is an enormous cluster..." 
		   << endl
		   << "warning: killsatellites: perhaps threshold is too low." 
		   << endl ;
	      exit (3) ;
	    }
	} 
	else 
	  {
	    if (npixels) 
	    {
	      push_back(current) ;
	      nextcolor++ ;
	    }
	  }
	
      }
    }
    
    
  }
}


bool ClusterList::Cut()
{
  cout << "Number of cluster before cut " << size() << endl;
  double xmean, ymean, x2mean, xymean, y2mean ;
  double a, b, theta, elongation ;
  double mm, mn, mo ;
  double mp;
  double ccddiagonal = sqrt(double(sqr(colors->Nx()) + sqr(colors->Ny()))) ;
  trackfound = false;
  
  for(ClusterIterator i = begin() ; i != end();) 
    {
      
      xmean = i->xsum / i->size ;
      ymean = i->ysum / i->size ;
      x2mean = i->x2sum / i->size - xmean * xmean ;
      xymean = i->xysum / i->size - xmean * ymean ;
      y2mean = i->y2sum / i->size - ymean * ymean ;
      
      mm = ((x2mean + y2mean) / 2) ;
      mn = ((x2mean - y2mean) / 2) ;
      mo = sqrt(mn*mn + xymean * xymean) ;
      
      a = 2.0 * sqrt(mm + mo) ; b = 2.0 * sqrt(mm - mo) ;
      
      elongation = a / b ;
      
      mp = x2mean - y2mean ;
      theta = 0.5 * atan2(2*xymean, mp) ; 
      
      if (debuglevel >= 2) {
	cout << "Cluster #" << i->color << "  size = " << i->size 
	     << "  xmean = " << xmean << "  ymean = " << ymean
	     << "  x2mean = " << x2mean 
	     << "  xymean = " << xymean 
	     << "  y2mean = " << y2mean 
	     << "  (A,B,theta) = ( " << a << " " << b << " " << theta << " )"
	     << "  elongation = " << elongation << endl ;
      }
      bool keep = false;

      if ((*i).size >= (ccddiagonal / 10)) 
	{
	  if ((elongation >=  2) && (a >= (ccddiagonal / 50))) 
	    //if (fabs(fabs(theta) - (MY_PI/2)) >= 0.01) 
	      {
		// To avoid dead columns
		cout << "Cluster #" << i->color << "  size = " << i->size 
		     << "  xmean = " << xmean << "  ymean = " << ymean
		     << "  x2mean = " << x2mean 
		     << "  xymean = " << xymean 
		     << "  y2mean = " << y2mean 
		     << "  (A,B,theta) = ( " << a 
		     << " " << b << " " << theta << " )"
		     << "  elongation = " << elongation << endl ;
		trackfound = true ;
		keep = true;
		//ColorTagged.push_back(i->color);
	      }
	}

      if(!keep)
	{
	  //cout << "erasing " << i->Color() << endl;
	  i = erase(i);
	}
      else ++i;
    }
  cout << "Number of cluster after cut " << size() << endl;

  return(trackfound);
}


Image ClusterList::Mask()
{
  cout<< "Constructing the mask " << endl;
  Image Mask(colors->Nx(), colors->Ny()) ;
  //FitsImage("colors.fits",*colors);

  if(!trackfound) return Mask;

  for (int j =0; j!=colors->Ny()-1; ++j)
    for (int i =0; i!=colors->Ny()-1; ++i)
      {
	if((*colors)(i,j)==0) continue;
	for ( ClusterCIterator it = begin(); it != end(); ++it)
	  {
	    //cout << it->Color() << endl;
	    if (it->Color() == (*colors)(i,j))
	      Mask(i,j) = 1;
	  }
      }
  cout<< "Mask Done" << endl;
  return(Mask);
  
}



