#include <sstream>
#include <fstream>

#include <poloka/reducedutils.h>
#include <poloka/photoratio.h>
#include <poloka/fitsimage.h>
#include <poloka/sestar.h>
#include <poloka/vutils.h>
#include <poloka/listmatch.h>
#include <poloka/polokaexception.h>

static double sq(const double &x) { return x*x; }

static const double EMAG_SCALE = 0.92103403719761845;

static bool SM_DecFluxMax(const StarMatch &S1, const StarMatch &S2)
{
  const SEStar *p1 = (const SEStar *) (const BaseStar *) S1.s1;
  const SEStar *p2 = (const SEStar *) (const BaseStar *) S2.s1;
  return  DecFluxMax(p1,p2);
}

static bool goodForRatio(const SEStar* star) {
  return (star->IsOK() && !star->IsSaturated());
}

static StarMatchList* PhotoRatioList(const ReducedImage& Im, const ReducedImage& Ref, const Gtransfo* Im2Ref) {

  GtransfoRef im2ref = GtransfoRef();
  if (!Im2Ref) {
    im2ref = FindTransfo(Im, Ref);
    if (!im2ref)
      throw("bad match between " + Im.Name() + " and " + Ref.Name());
  } else 
    im2ref = Im2Ref;
  
  if (!Im.HasCatalog() || !Ref.HasCatalog())
    throw PolokaException("no similar catalogues available for " + Im.Name() + " or " + Ref.Name());

  SEStarList imList(Im.CatalogName()); BaseStarList *bimList  = SE2Base(&imList);
  SEStarList refList(Ref.CatalogName()); BaseStarList *brefList = SE2Base(&refList); 
  
  const size_t nmax = 500;
  const double maxDist = 1.;

  StarMatchList *matchList = ListMatchCollect(*bimList, *brefList, im2ref, maxDist);
  matchList->sort(&SM_DecFluxMax);
  for (StarMatchIterator it = matchList->begin(); it != matchList->end() && matchList->size() < nmax; )  {
    if (goodForRatio(it->s1) && goodForRatio(it->s2))
      ++it;
    else
      it = matchList->erase(it);
  }
  
  if (matchList->empty())
    throw PolokaException("empty match list between " + Im.Name() + " and " + Ref.Name());
  
  return matchList;
}

static void LoadPhotoRefCatalog(const string& RefCatalogFile,
				const ReducedImage& Im,
				BaseStarList& BList,
				const Gtransfo* Im2Cat=0) {
  
  GtransfoRef raDec2Pix = Im2Cat;
  if (!raDec2Pix) {
    raDec2Pix = Im.RaDecToPixels();
    if (!raDec2Pix)
      throw PolokaException("no good WCS for " + Im.Name());
  }

  GtransfoRef pix2RaDec = Im2Cat;
  if (!pix2RaDec) {
    pix2RaDec = Im.PixelsToRaDec();
    if (!pix2RaDec)
      throw PolokaException("no good invert WCS for " + Im.Name());
  }
  
  Frame pixFrame = Im.UsablePart();
  Frame raDecFrame = ApplyTransfo(pixFrame, *pix2RaDec);

  ifstream in(RefCatalogFile.c_str());
  char c;
  string line;

  while (in >> c) {
    in.unget();
    if (c == '@' || c == '#') continue;
    if (!getline(in, line)) break;
    double ra,dec,mag,x,y,emag;
    istringstream iline(line);
    if (!(iline >> ra >> dec >> mag)) continue;
    if (!raDecFrame.InFrame(ra,dec)) continue;
    raDec2Pix->apply(ra, dec, x, y);
    if (!pixFrame.InFrame(x,y)) continue;
    double flux = pow(10., 0.4*(ZP_REF-mag));
    double eflux = 0;
    if (iline >> emag) eflux = EMAG_SCALE * emag * flux;
    BaseStar *star = new BaseStar(x, y, flux, eflux);
    BList.push_back(star);
  }
}

/*! Notes about photometric ratios:

   Measuring a flux ratio is not that straightforward. The "obvious"
   least-squares approach is for example facing serious problems.
   Imagine you want to fit y=ax on a sample of x's and y's, i.e. estimate the
   ratio a. If you minimize sum_i(y_i-a * x_i)^2/sy_i^2, the estimate
   of a will be unbiased only if there are no uncertainties on x.
   In the opposite case, there is a bias : <a_min> = <a_true> *
   sum(x_i^2)/ sum(x_i^2+sx_i^2). There is no single way to get around
   this problem. Getting back to photometric ratios, there is usually
   no problem about optimality, because the statistical accuracy is
   usually far too good, compared to systematics (e.g. non-uniformity
   of response due to PSF and/or flatfielding).

   AveragePhotoRatio computes a weighted average of the ratios. This is
   less biased than straight LS, but still biased because <1/x> !=
   1/<x> for symetric pdf's. One improvement is to compute both ratios
   (r1 and r2); if the estimator is unbiased, then <r1*r2> = 1; one
   way to compensate for a bias of r1 and r2 is to return the
   geometric mean sqrt(r1/r2); This is not coded (yet).

   One estimator with a very small bias, and almost optimal is to use
   the functional relationship framework. Namely, you use
   least-squares:

      chi2 = sum (y_i-a xp_i)^2/sy_i^2 + (x_i-xp_i)^2/sx_i^2

    where a and xp_i are parameters. When you solve the normal
    equations, you find that it is equivalent to minimize:

       chi2 = sum (y_i - a x_i)^2/(sy_i^2+(a sx_i)^2).

    The minimization is tedious (non linearities are important), but
    doable. This is done in TLSPhotoRatio
*/

double AveragePhotoRatio(const StarMatchList& MatchList, double& Error) {

  double sumwt=0, sum=0;
  for (StarMatchCIterator it = MatchList.begin(); it != MatchList.end(); ++it) {      
    const BaseStar *cur = it->s1;
    const BaseStar *ref = it->s2;
    const double ratio = cur->flux / ref->flux;
    const double weight = 1./(sq(ratio*(cur->eflux/cur->flux + ref->eflux/ref->flux)));
    sumwt += weight;
    sum += weight * ratio;
  }

  if (sumwt == 0 || sum == 0)
    throw PolokaException(" only found zero fluxes and errors");
  
  Error = 1/sqrt(sumwt);
  return sum / sumwt;
}

double AveragePhotoRatio(const ReducedImage &Im,
			 const ReducedImage &Ref, 
			 double &Error,
			 const Gtransfo* Im2Ref) {


  StarMatchList* matchList = PhotoRatioList(Im, Ref, Im2Ref);
  double ratio = AveragePhotoRatio(*matchList, Error);
  delete matchList;
  return ratio;
}

double MedianPhotoRatio(const StarMatchList& MatchList) {

  // This is unweighted. If we weight it, we might consider 
  // averaging between the direct and reversed photom ratio

  double *ratios = new double[MatchList.size()];
  double *pr = &ratios[0];
  for (StarMatchCIterator it = MatchList.begin(); it != MatchList.end(); ++it)
    *pr++ = it->s1->flux / it->s2->flux;
  double med = DArrayMedian(ratios, MatchList.size());
  delete [] ratios;
  return med;
}

double MedianPhotoRatio(const StarMatchList& MatchList, double& Error) {

  // This is unweighted. If we weight it, we might consider 
  // averaging between the direct and reversed photom ratio

  vector<double> ratios(MatchList.size());
  vector<double>::iterator ir = ratios.begin();
  for (StarMatchCIterator it = MatchList.begin(); it != MatchList.end(); ++it, ++ir)
    *ir = it->s1->flux / it->s2->flux;
  return median_mad(ratios, Error);
}


double MedianPhotoRatio(const ReducedImage &Im,
			const ReducedImage &Ref, 
			double &Error,
			const Gtransfo* Im2Ref) {


  StarMatchList* matchList = PhotoRatioList(Im, Ref, Im2Ref);
  double ratio = MedianPhotoRatio(*matchList, Error);
  delete matchList;
  return ratio;
}

//! see notes above if you wonder what is done here.
double TLSPhotoRatio(const StarMatchList& MatchList,
		    double& Error,
		    const double& NSigChi2Cut) {
  int niter = 0;
  double chi2;
  double chi2Old = 1e30;
  size_t oldCount = MatchList.size();
  double *chi2Vals = new double[oldCount];
  double chi2Cut = 1e30;
  const int nitermax = 50;

  double ratio = MedianPhotoRatio(MatchList);
  double var;

  do {
    double num = 0;
    double deno = 0;
    chi2 = 0;
    int count = 0;
    for (StarMatchCIterator it = MatchList.begin(); it != MatchList.end(); ++it) {
      double y = it->s1->flux;
      double x = it->s2->flux;
      double sx2 = sq(it->s2->eflux);
      double sy2 = sq(it->s1->eflux);
      double d = (sy2 + sq(ratio)*sx2);
      double xp = (ratio*y*sx2+x*sy2)/d;
      double thisChi2 = sq(y-ratio*x)/d;
      if (thisChi2 > chi2Cut) continue;
      chi2Vals[count++] = thisChi2;
      chi2 += thisChi2;
      num += xp*(y-ratio*xp)/sy2;
      deno += (sq(xp)-sx2*sq(y-2.*ratio*xp)/d)/sy2;
    }

    // no protection for dividing by zero
    if (fabs(deno) < 1.e-20)
      throw PolokaException("dividing by zero");

    ratio += num/deno;
    var = 1/deno;
    double chi2Mean, chi2Med, chi2Sig;
    Dmean_median_sigma(chi2Vals, count, chi2Mean, chi2Med, chi2Sig);
    chi2Cut = chi2Med+NSigChi2Cut*chi2Sig;
    // in case we got twice the same list as input
    if (chi2Cut <1e-2) chi2Cut = 1e-2;
    bool outliers = false;
    for (size_t k=0; k < count; ++k)
      if (chi2Vals[k] > chi2Cut) { outliers = true; break; }
    if (!outliers && oldCount == count && chi2Old - chi2 < 1e-3) break;
    chi2Old = chi2;
    oldCount = count;
    if (outliers) niter++;
  } while (niter < nitermax);

  delete [] chi2Vals;
  
  if (isnan(ratio))
    throw PolokaException("photometric ratio is NaN value");

  if (niter>=nitermax)
    throw PolokaException("reached max number of iterations without converging");

  Error = sqrt(var);
  return ratio;
}

double TLSPhotoRatio(const ReducedImage &Im,
		     const ReducedImage &Ref, 
		     double &Error,
		     const Gtransfo* Im2Ref) {


  StarMatchList* matchList = PhotoRatioList(Im, Ref, Im2Ref);
  if (!matchList)
    throw PolokaException("could not find a decent match");
  double ratio = TLSPhotoRatio(*matchList, Error);
  delete matchList;
  return ratio;
}

double TLSPhotoRatio(const ReducedImage &Im,
		     const string& RefCatalogFile,
		     double &Error,
		     const Gtransfo* Im2Cat) {

  if (!Im.HasCatalog())
    throw PolokaException("no catalog available for " + Im.Name());

  if (!FileExists(RefCatalogFile))
    throw PolokaException("reference catalog " + RefCatalogFile + " not valid");

  BaseStarList refList;

  LoadPhotoRefCatalog(RefCatalogFile, Im, refList, Im2Cat);
  
  SEStarList imList(Im.CatalogName());
  BaseStarList *bimList = SE2Base(&imList);
  StarMatchList *matchList = ListMatchCollect(*bimList, refList, 2);
  if (!matchList || matchList->empty())
    throw PolokaException("could not find a decent match");
  double ratio = TLSPhotoRatio(*matchList, Error);
  delete matchList;
  return ratio;
}

double ZpPhotoRatio(const double& Zp, const double& SigZp, 
		    const double& ZpRef, const double& SigZpRef,
		    double& Sig) {

  double ratio = pow(10., 0.4*(ZpRef - Zp));
  // var(r)/ r = (log(10)*0.4) ** 2
  Sig = ratio * EMAG_SCALE * sqrt(sq(SigZp) + sq(SigZpRef));
  return ratio;
}

double ZpPhotoRatio(const ReducedImage& Im, const ReducedImage& Ref, double& Error) {
  
  // use fitsheader until zero points in reducedimage stabilize...
  if (!Im.HasImage() || !Ref.HasImage())
    throw PolokaException("missing fits files for "+Im.Name()+" or "+Ref.Name());

  FitsHeader imhead(Im.FitsName());
  FitsHeader refhead(Ref.FitsName());
  
  double errdef = 0.05;
  Error = errdef;

  if (imhead.HasKey("ZP_PHOT") && refhead.HasKey("ZP_PHOT"))
    return ZpPhotoRatio(imhead.KeyVal("ZP_PHOT"),
			imhead.KeyVal("EZP_PHOT"),
			refhead.KeyVal("ZP_PHOT"),
			refhead.KeyVal("EZP_PHOT"),
			Error);
  else if (imhead.HasKey("ZP") && refhead.HasKey("ZP"))
    return ZpPhotoRatio(imhead.KeyVal("ZP"),
			errdef,
			refhead.KeyVal("ZP"),
			errdef,
			Error);
  else if (imhead.HasKey("ZPTOADS") && refhead.HasKey("ZPTOADS"))
    return ZpPhotoRatio(imhead.KeyVal("ZPTOADS"),
			errdef,
			refhead.KeyVal("ZPTOADS"),
			0,
			Error);  
  else if (imhead.HasKey("ZEROUSNO") && refhead.HasKey("ZEROUSNO"))
    return ZpPhotoRatio(imhead.KeyVal("ZEROUSNO"),
			imhead.KeyVal("DZEROUSN"),
			refhead.KeyVal("ZEROUSNO"),
			refhead.KeyVal("DZEROUSN"),
			Error);
  else
    throw PolokaException("could not find consistent zero point in headers of " + Im.Name() + " and " + Ref.Name());
}

double ZpPhotoRatio(const ReducedImage& Im, const double& ZpRef) {
  
  // use fitsheader until zero points in reducedimage stabilize...
  if (!Im.HasImage())
    throw PolokaException("missing fits file for "+Im.Name());

  FitsHeader imhead(Im.FitsName());
  
  double err, errdef=0.05;

  if (imhead.HasKey("ZP_PHOT"))
    return ZpPhotoRatio(imhead.KeyVal("ZP_PHOT"),
			imhead.KeyVal("EZP_PHOT"),
			ZpRef,
			0.,
			err);
  else if (imhead.HasKey("ZPTOADS"))
    return ZpPhotoRatio(imhead.KeyVal("ZPTOADS"),
			errdef,
			ZpRef,
			0.,
			err);  
  else if (imhead.HasKey("ZP"))
    return ZpPhotoRatio(imhead.KeyVal("ZP"),
			errdef,
			ZpRef,
			0.,
			err);
  else if (imhead.HasKey("ZEROUSNO"))
    return ZpPhotoRatio(imhead.KeyVal("ZEROUSNO"),
			imhead.KeyVal("DZEROUSN"),
			ZpRef,
			0.,
			err);
  else
    throw PolokaException("could not find zero point in header of " + Im.Name());
}


double PhotoRatio(const ReducedImage& Im, const ReducedImage& Ref, double& Error,
		  const Gtransfo* Im2Ref, const PhotoScalingMethod Method) {

  if (Im == Ref) {
    Error = 0.;
    return 1.; 
  }
  
  switch (Method) {
  case ZeroPointDiff:
    return ZpPhotoRatio(Im, Ref, Error);
  case TotalLeastSquares:
    return TLSPhotoRatio(Im, Ref, Error, Im2Ref);
  case AverageRatio:
    return AveragePhotoRatio(Im, Ref, Error, Im2Ref);
  case NoScaling:
    Error = 0.; return 1.;
  default:
    return MedianPhotoRatio(Im, Ref, Error, Im2Ref);      
  }
  return MedianPhotoRatio(Im, Ref, Error, Im2Ref);
}

