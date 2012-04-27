#include "photoratio.h"
#include "fitsimage.h"
#include "sestar.h"
#include "vutils.h"
#include "toadscards.h"
#include "listmatch.h"
#include "polokaexception.h"

static double sq(const double &x) { return x*x; }

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

  string imListName = Im.CatalogName();
  string refListName = Ref.CatalogName();

  if (!Im.HasCatalog() || !Ref.HasCatalog())
    throw PolokaException("no similar catalogues available for " + Im.Name() + " or " + Ref.Name());

  SEStarList imList(imListName); BaseStarList *bimList  = SE2Base(&imList);
  SEStarList refList(refListName); BaseStarList *brefList = SE2Base(&refList); 

  if (!Im2Ref)
    Im2Ref = ListMatch(*bimList, *brefList, MatchConditions(DefaultDatacards()));
  if (!Im2Ref)
    throw("bad match between " + Im.Name() + " and " + Ref.Name());

  const size_t nmax = 500;
  const double maxDist = 1.;

  StarMatchList *matchList = ListMatchCollect(*bimList, *brefList, Im2Ref, maxDist);
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

       chi2 = sum (y_i - a x_i)^2/(sy_i^2+(a x_i)^2).

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


static double median_mad(vector<double>& x, double& disp) {
  size_t n = x.size();
  sort(x.begin(), x.end());
  double med = (n & 1) ? x[n/2] : (x[n/2-1] + x[n/2])*0.5;  
  for (vector<double>::iterator it = x.begin(); it != x.end(); ++it) {
    *it = fabs(*it - med);
  }
  sort(x.begin(), x.end());
  double mad = (n & 1) ? x[n/2] : (x[n/2-1] + x[n/2])*0.5;  
  disp = 1.4826 * mad; // robust estimator of standard deviation
  return med;
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
  double ratio = TLSPhotoRatio(*matchList, Error);
  delete matchList;
  return ratio;
}


double ZpPhotoRatio(const double& Zp, const double& SigZp, 
		    const double& ZpRef, const double& SigZpRef,
		    double& Sig) {

  double ratio = 1;
  ratio = pow(10.,0.4*(ZpRef - Zp));
  // var(r)/ r = (log(10)*0.4) ** 2
  Sig = ratio * 0.92103403719761845 * sqrt(sq(SigZp)+sq(SigZpRef));
  return ratio;
}

double ZpPhotoRatio(const ReducedImage& Im, const ReducedImage& Ref, double& Error) {
  
  // use fitsheader until zero points in reducedimage stabilize...
  if (!Im.HasImage() || !Ref.HasImage())
    throw PolokaException("missing fits files for "+Im.Name()+" or "+Ref.Name());

  FitsHeader imhead(Im.FitsName());
  FitsHeader refhead(Ref.FitsName());
  
  Error = 0.05; // default error

  if (imhead.HasKey("ZP_PHOT") && imhead.HasKey("ZP_PHOT"))
    return ZpPhotoRatio(imhead.KeyVal("ZP_PHOT"),
			imhead.KeyVal("EZP_PHOT"),
			refhead.KeyVal("ZP_PHOT"),
			refhead.KeyVal("EZP_PHOT"),
			Error);
  else if (imhead.HasKey("ZPTOADS") && imhead.HasKey("ZPTOADS"))
    return ZpPhotoRatio(imhead.KeyVal("ZPTOADS"),
			0,
			refhead.KeyVal("ZPTOADS"),
			0,
			Error);  
  else if (imhead.HasKey("ZP") && imhead.HasKey("ZP"))
    return ZpPhotoRatio(imhead.KeyVal("ZP"),
			0,
			refhead.KeyVal("ZP"),
			0,
			Error);
  else if (imhead.HasKey("ZEROUSNO") && imhead.HasKey("ZEROUSNO"))
    return ZpPhotoRatio(imhead.KeyVal("ZEROUSNO"),
			imhead.KeyVal("DZEROUSN"),
			refhead.KeyVal("ZEROUSNO"),
			refhead.KeyVal("DZEROUSN"),
			Error);
  else
    throw PolokaException("could not find decent zero point for "+Im.Name()+" or "+Ref.Name());
}

double PhotoRatio(const ReducedImage& Im, const ReducedImage& Ref, double& Error,
		  const Gtransfo* Im2Ref, const PhotoRatioMethod Method) {

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
  }

  return MedianPhotoRatio(Im, Ref, Error, Im2Ref);
}
