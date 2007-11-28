#include <iostream>
#include "string"
#include "reducedimage.h"
#include "polokaexception.h"

using namespace std;

#include "fitsslice.h"
#include "apersestar.h"
#include "histopeakfinder.h"
#include "fastfinder.h"
#include "wcsutils.h"
#include "dicstar.h"
#include "gtransfo.h"

/* What is this code intended to do?
   produce a clean list of isolated stars in the deep fields.
   1) first step : select stars from their shape on the "large references"
     1.1) measure all second moments of objects detected by SExtractor
          I use mkcat2 for that
     1.2 ) split the catalog into tiles and find the star locus
        in each of them (with a rather high S/N cut)
     1.3) smooth the PSF second moment variation as a function of position
          with a polynomial
     1.4) find objects close to the smoothed PSF second moments
          in all the image (these are hopefully the stars) 
     1.5) restart from 1.2 with the output from 1.4 (actually, the first time
          I use half the S/N cut that I use the second time).

   2) merge with tertiaries
      bright stars are saturated on stacks and hence do not provide 
      reliable 2nd moments. There is however no very good reason to drop
      them at this level, this is why objects brighter than 18 present
      in the tertiaries catalogs and isolated are added to the list from 
      the previous steps

    3) outputs : 
          - a "short" star file with only a few scores (~15) per star
          - the same entries as apersestars.
          - the 2nd moments found on tiles and the smoothed version (shapes.list)
                 
*/





#ifdef STORAGE
static bool ComputeSecondMoments(const ReducedImage &Ri, AperSEStarList &Cat)
{
  //this code now "slices" the input images in order to accomodate monsters

  int ySliceSize = 500;
  int sliceOverlap = 110;


  // first pass on the pixels : compute shapes (and positions)
  {

    FitsParallelSlices slices(ySliceSize,sliceOverlap);

    slices.AddFile(Ri.FitsName());

    if (!bool(slices[0]->KeyVal("BACK_SUB")))
      {
	cout << " ReducedImage::MakeAperCat : this code cannot accomodate images with background left" << endl;
	return false;
      }
    
    slices.AddFile(Ri.FitsWeightName());
    
    
    double yStarMin = 0;
    do  // loop on image slices
      {
	double yStarMax = (slices.LastSlice())? slices.ImageJ(ySliceSize): slices.ImageJ(ySliceSize)-55;
	
	double offset = slices.ImageJ(0);
	for (AperSEStarIterator i = Cat.begin(); i != Cat.end(); ++i)
	  {
	    AperSEStar &s = **i;
	    if (s.y < yStarMax && s.y >= yStarMin)
	      {
		s.y -= offset;
		s.ComputeShapeAndPos(*slices[0],*slices[1]);
		s.y += offset;
	      }
	  }
	yStarMin = yStarMax;
      }
    while(slices.LoadNextSlice());
  } // end of first pass
  return true;
}
#endif

class StarShape
{
 public :
  double x,y;
  double mxx,myy,mxy;

  StarShape(double X, double Y, double Mxx, double Myy, double Mxy) 
    : x(X), y(Y), mxx(Mxx), myy(Myy), mxy(Mxy) {};

};


#include "poly2.h"

class StarShapeList : public list<StarShape>
{
  Poly2 mxx_model, myy_model, mxy_model;
  

 public :
  bool Fit(const Frame &F, const int Deg = 4);
  void Write(const string& FileName) const;
  double MxxModel(const double &x, const double &y) const
  {  return mxx_model.Value(x,y);}
  double MyyModel(const double &x, const double &y) const
  {  return myy_model.Value(x,y);}
  double MxyModel(const double &x, const double &y) const
  {  return mxy_model.Value(x,y);}

  bool FindStars(const AperSEStarList &Cat,
		 AperSEStarList &Stars, double MinSN) const;

  void dump(ostream& s) const;


};

void StarShapeList::dump(ostream &s) const
{
  s << "# x: " << endl;
  s << "# y: " << endl;
  s << "# gmxx: measured " << endl;
  s << "# gmyy: " << endl;
  s << "# gmxy: " << endl;
  s << "# gmxxs: fit " << endl;
  s << "# gmyys: " << endl;
  s << "# gmxys: " << endl;
  s << "#end" << endl;
  for (const_iterator i = begin(); i!= end(); ++i)
    {
      const StarShape &ss = *i;
      s << ss.x << ' ' << ss.y << ' ' 
	<< ss.mxx << ' ' 
	<< ss.myy << ' ' 
	<< ss.mxy << ' '
	<< MxxModel(ss.x, ss.y) << ' ' 
	<< MyyModel(ss.x, ss.y) << ' ' 
	<< MxyModel(ss.x, ss.y)	<< endl;
    }
}

bool StarShapeList::Fit(const Frame &F, const int Deg)
{

  Poly2 model(F.xMin, F.yMin, F.xMax, F.yMax, Deg);
  unsigned s = model.NTerms();
  Mat A(s,s);
  Mat B(s,3);
  Vect monomials(s);
  for (const_iterator i = begin(); i!= end(); ++i)
    {
      const StarShape &ss = *i;
      model.Monomials(ss.x,ss.y,monomials);
      for (unsigned k1=0; k1<s; ++k1)
	{
	  for (unsigned k2=0; k2<s; ++k2)
	    A(k1,k2)+=monomials(k1)*monomials(k2);
	  B(k1,0)+= monomials(k1)*ss.mxx;
	  B(k1,1)+= monomials(k1)*ss.myy;
	  B(k1,2)+= monomials(k1)*ss.mxy;
	}
    }
  for (unsigned k=0; k<3;++k)
    {
      Mat AA(A);
      Vect BB(s);
      for (unsigned k1=0; k1<s; ++k1) BB(k1) = B(k1,k);
      if (cholesky_solve(AA,BB,"U") != 0)
	{
	  return false;
	}
      Poly2 *p;
      if (k==0) p = &mxx_model;
      if (k==1) p = &myy_model;
      if (k==2) p = &mxy_model;
      *p = model;
      p->SetCoeffs(BB);
    }
  return true;      
}

#include <fstream>
void StarShapeList::Write(const string& FileName) const
{
  ofstream s(FileName.c_str());
  mxx_model.Write(s);
  myy_model.Write(s);
  mxy_model.Write(s);
  s.close();
}



static double sq(const double &x) {return x*x;}


#include <algorithm>
static StarShapeList FindStarShapes(const AperSEStarList &Cat, 
				    unsigned NxSlice=9, 
				    unsigned NySlice=9)
{
  StarShapeList l;
  double xmin=1e30,xmax = -1e30, ymin = 1e30, ymax=-1e30;
  for (AperSEStarCIterator i = Cat.begin(); i != Cat.end(); ++i)
    {
      const AperSEStar &s = **i;
      xmin = min(xmin,s.x);
      xmax = max(xmax,s.x);
      ymin = min(ymin,s.y);
      ymax = max(ymax,s.y);
    }
  double xstep = (xmax-xmin)/NxSlice;
  double ystep = (ymax-ymin)/NySlice;
  for (unsigned kx=0;kx<NxSlice; kx++)
    {
      double xmins = xmin+kx*xstep;
      double xmaxs = xmins+xstep;
      Frame ySliceFrame(xmins,ymin,xmaxs,ymax);
      AperSEStarList ySliceCat;
      Cat.ExtractInFrame(ySliceCat, ySliceFrame);
      for (unsigned ky=0; ky<NySlice; ky++)
	{
	  double ymins = ymin+ky*ystep;
	  double ymaxs = ymins+ystep;
	  Frame tile(xmins,ymins,xmaxs,ymaxs);
	  AperSEStarList tileCat;
	  ySliceCat.ExtractInFrame(tileCat, tile);
	  double xSize,ySize, corr;
	  StarScoresList scores;
	  if (FindStarShapes(tileCat, 20, xSize, ySize, corr, scores))
	    {
	      double mxx = sq(xSize);
	      double myy = sq(ySize);
	      double mxy = corr*xSize*ySize;
	      l.push_back(StarShape(0.5*(xmaxs+xmins),
				    0.5*(ymaxs+ymins),
				    mxx,myy,mxy));
	    }
	}
    }
  Frame wholeFrame(xmin,ymin,xmax,ymax);
  l.Fit(wholeFrame);
  return l;
}




#include "histo2d.h"
bool StarShapeList::FindStars(const AperSEStarList &Cat,
			      AperSEStarList &Stars, double MinSN) const
{
  // inspired from FindStarShapes in apersestar.cc

  Histo2d histo(30,-5., 5., 30, -5., 5.);
  double NShapeSig = 5;


  StarScoresList scores;
  for (AperSEStarCIterator i = Cat.begin(); i != Cat.end(); ++i)
    {
      const AperSEStar &s = **i;
      /* basic quality cut : gflag = 0 means 2nd moments 
	 (gmxx & co) computed sucessfully. Flag is the SExtractor 
	 extraction flags, and Flagbad refers to bad pixels.
      */
      if (s.gflag || s.FlagBad()) continue;
      /* if the s/n is poor, the measurement of 2nd moments
	 is also poor: cut on S/N before entering into the
	 "star clump finder".
      */
      if (s.flux<0 || s.EFlux() < 0 || s.flux< MinSN * s.EFlux()) continue;

      double xx = sqrt(s.gmxx)  - sqrt(MxxModel(s.x,s.y));
      double yy = sqrt(s.gmyy) -  sqrt(MyyModel(s.x,s.y));
      histo.Fill(xx,yy, s.Fluxmax());
      scores.push_back(new StarScores (xx,yy,&s));
    }

  unsigned nobj = scores.size();
  if (nobj<5)
    {
      cout << " FindStars : only " << nobj << " objects passing the cuts " << endl;
      return false;
    }

  double xMax, yMax;
  histo.MaxBin(xMax, yMax);
  Ellipse ellipse;
  bool ok = HistoPeakFinder(scores, histo, xMax, yMax, ellipse);
  cout << "average ellipse " << ellipse;

  for (StarScoresCIterator i = scores.begin(); i != scores.end(); ++i)
    {
      const StarScores &scores = **i;
      if (scores.nSig<NShapeSig) // it lies within the star cluster
	{
	  const AperSEStar &a = 
	    dynamic_cast< const AperSEStar &>(*scores.star);
	  Stars.push_back(&a);
	}
    }
  cout << " found " << Stars.size() << " stars " << endl;
  GlobalVal &glob = Stars.GlobVal();
  glob = Cat.GlobVal(); // copy input keys
  glob.AddKey("STARSIGSHAPECUT", NShapeSig);
  glob.AddKey("STARSIGSNCUT", MinSN);
  glob.AddKey("PRODUCER", "extract_stars_from_large_ref");
  return ok;
}


static void SelectIsolatedStars(AperSEStarList &Stars, const unsigned AperRank,
				const double MaxFrac, const double Zp)
{
  int polluted = 0;
  for (AperSEStarIterator si = Stars.begin(); si != Stars.end(); )
    {
      AperSEStar &a = **si;
      const Aperture &aper = a.apers[AperRank];
      if (aper.fother > a.flux*MaxFrac)
	{
	  polluted++;
#ifdef STORAGE
	  if (polluted == 1) 
	    cout << " polluted stars (x,y,mag,frac):"  << endl;
	  cout << a.x << ' ' << a.y << ' ' 
	       << -2.5*log10(a.flux)+Zp << ' '
	       <<aper.fother/a.flux << endl;
#endif
	  si = Stars.erase(si);
	  continue;
	}
      ++si;
    }
  cout << " dropped " << polluted << " polluted stars with a cut at " << MaxFrac << endl;
}

#include <map>

struct MapSE : public map<const BaseStar *,int>
{
  int operator() (const BaseStar * s) const
  {
    MapSE::const_iterator i = find(s);
    if (i == end()) return -1;
    return i->second;
  }
};



#include <iomanip>
static void WriteSummaryCatalog(const string &FileName, 
				const AperSEStarList &L,
				const MapSE &Origin,
				const StarShapeList &Ssl,
				const Gtransfo &Wcs,
				const double Zp)
{
  ofstream s(FileName.c_str());

  vector<string> globs(L.GlobVal().OutputLines());
  for (unsigned k=0; k<globs.size() ; ++k) s << '@' << globs[k] << endl;

  s << "# ra :" << endl;
  s << "# dec :" << endl;
  s << "# x: coord in large image" << endl;
  s << "# y: coord in large image" << endl;
  s << "# mag : not precise at all" << endl;
  s << "# flag : Sextractor extraction flag on large ref" << endl;
  s << "# pol : pollution by other objects in aperture (flux fraction)" << endl;
  s << "# origin : 1=stars in largeRef, 2=first iteration tertiaries, 3=both" << endl;
  s << "# gmxx : second moment of this star " << endl;
  s << "# gmyy :" << endl;
  s << "# gmxy :" << endl;
  s << "# gmxxs : spatially smoothed gmxx" << endl;
  s << "# gmyys :" << endl;
  s << "# gmxys :" << endl;
  s << "# end" << endl;
  for (AperSEStarCIterator i = L.begin(); i!= L.end(); ++i)
    {
      const AperSEStar &star = **i;
      Point raDec = Wcs.apply(star);
      double mag = (star.flux) ? -2.5*log10(star.flux)+Zp : 99; 
      s << setprecision(10) << raDec.x << ' ' << raDec.y << ' ';
      s << setprecision(9) << star.x << ' ' << star.y << ' ';
      s << setprecision(6) << mag << ' ' << star.Flag() << ' ';
      s << star.apers[6].fother/star.flux << ' ';
      s << Origin(&star) << ' ';
      s << star.gmxx << ' ' << star.gmyy << ' ' << star.gmxy << ' ';
      s << Ssl.MxxModel(star.x, star.y) << ' ';
      s << Ssl.MyyModel(star.x, star.y) << ' ';
      s << Ssl.MxyModel(star.x, star.y) << ' ';
      s << endl;
    }
  s.close();
}
      
    



static void usage(const char *prog)
{
  cerr << "usage : " << endl
       << prog << "  -r <large ref> [-t <tertiaries>] -o <outFileName>" << endl;
  cerr << " -p (output in pix coordinates) " << endl;
  cerr << " -m <minSig2Noise> " << endl;
  cerr << " -f <maxFrac> (max fraction of aperture flux due to other stuff)" << endl;
  exit(EXIT_FAILURE);
}


int main(int nargs,char **args)
{
  string riName;
  string tertiariesName;
  string outFileName;
  double minSN = 100;
  double maxFrac = 0.01; // max fraction of flux in aperture attributed to other stuff

  for (int i=1; i<nargs;++i)
    {
      char *arg = args[i];
      if (arg[0] == '-')
	switch (arg[1])
	  {
	  case 'r' : i++;riName = args[i]; break;
	  case 't' : i++;tertiariesName = args[i]; break;
	  case 'o' : i++;outFileName = args[i]; break;
	  case 'm' : i++; minSN = atof(args[i]); break;
	  case 'f' : i++; maxFrac = atof(args[i]); break;
	  default : usage(args[0]);
	  }
      else usage(args[0]);
    }
  if ((riName == "") || (outFileName == "")) 
    usage(args[0]);
  bool ok = true;
  try
    {
      ReducedImage ri(riName);
      // load wcs to make sure it is there 
      FitsHeader head(ri.FitsName());
      Gtransfo *wcs;
       if (!WCSFromHeader(head, wcs))
	 throw(PolokaException(" No valid WCS in " + ri.FitsName()));

       // select stars in aper catalog of large ref
      string aperseName(ri.AperCatalogName());
      cout << " reading " << aperseName << "..." << endl;
      AperSEStarList apercat(aperseName);
      // run star finder on tiles
      StarShapeList ssl1 = FindStarShapes(apercat);
      AperSEStarList stars1;
      // select stars based on local PSF 2nd moments
      ok &= ssl1.FindStars(apercat,stars1,minSN/2.);
      // since it is a large ref, there are lot of galaxies: rerun
      // the "PSF" modeling and star finder
      StarShapeList ssl = FindStarShapes(stars1);
      AperSEStarList stars;
      // find stars based on local PSF 2nd moments
      ok &= ssl.FindStars(stars1,stars,minSN);

      /* The output may mix entries from the large ref catalog and
       entries from a previous tertiary catalog (if provided in
       input). To keep track of where the stars come from we use a bit
       pattern (stored into a MapSE) */
      MapSE origin;
      for (AperSEStarIterator i = stars.begin(); i != stars.end(); ++i)
	origin[*i] = 1;

      // build a name like "mr" for mag in r band
      string magName("m"+StringToLower(string(head.KeyVal("TOADBAND"))));      	
      if (tertiariesName != "")
	{
	  DicStarList tertiaries(tertiariesName);
	  Gtransfo *inverseWCS = wcs->InverseTransfo(0.01,Frame(head));
	  tertiaries.ApplyTransfo(*inverseWCS);
	  FastFinder starFinder(AperSE2Base(stars));
	  FastFinder aperFinder(AperSE2Base(apercat));
	  double maxDist = 3; // pixels
	  double magSplit = 18; // tertiaries fainter than this are ignored
	  // now all 3 lists are in pixel coordinates
	  /* for each bright tertiary:
	     o  check if it is already in "stars"
	     o  if not, find it in the large apercat
	     o  if found, copy it into "stars"
	  */
	  unsigned missed = 0;
	  for (DicStarCIterator ti=tertiaries.begin(); 
	       ti != tertiaries.end(); ++ti)
	    {
	      const DicStar &t = **ti;
	      const BaseStar* bs = starFinder.FindClosest(t,maxDist);
	      if (bs) // is already in "stars"
		{
		  origin[bs] = 3;
		  continue; 
		}
	      /* in the tertiaries, there are objects with mag = 0
		 ignore them here */
	      double terMag = t.getval(magName);
	      if (terMag > magSplit || terMag < 1) continue;
	      bs = aperFinder.FindClosest(t,maxDist);
	      if (!bs)
		{
		  missed++;
		  continue;
		}
	      const AperSEStar *ap = dynamic_cast<const AperSEStar*>(bs);
	      stars.push_back(ap);
	      origin[ap] = 2;
	    }
	  cout << " number of bright tertiairies not found in " << aperseName 
	       << " : " << missed << endl;
	}// end of (if tertiaries...)
      unsigned aperRank = 6;
      double zp = head.KeyVal("ZP");
      SelectIsolatedStars(stars, aperRank, maxFrac, zp);
      cout << "writing " << stars.size() << " stars to " 
	   << outFileName << endl;
      stars.write(outFileName+".full");
      WriteSummaryCatalog(outFileName, stars, origin, ssl, *wcs, zp); 
      ofstream shapes("shapes.list");
      ssl.dump(shapes);
      shapes.close();
    }
  catch(PolokaException p)
    {
      p.PrintMessage(cout);
      ok = false;
    }
  return ( ok? EXIT_SUCCESS : EXIT_FAILURE);
}

  
