#include <iostream>
#include <string>

#include "apersestar.h"
#include "polokaexception.h"



struct WeightedVal
{
  double val;
  double w;
  WeightedVal( const double &V, const double &W): val(V), w(W) {};
  static bool CompareVal(const WeightedVal &L, const WeightedVal &R)
  { return L.val < R.val;}

};



#include <vector>


typedef std::vector<WeightedVal> WeightedVals;

static double sq(const double &x) { return x*x;}


double WeightedMedian(WeightedVals &Vals, double &Sig, const double NSig)
{
  if (Vals.size() == 0) {Sig= -1; return -1;}
  sort(Vals.begin(), Vals.end(), WeightedVal::CompareVal);
  double med = 0;
  Sig = 1e30;
  double oldsig = Sig;

  for (int iter=0; iter<10; ++iter)
    {

      unsigned ilow = 0;
      unsigned ihigh = Vals.size() -1;
      if (iter>0)
	{
	  double cut = med-NSig*Sig;
	  while (Vals[ilow].val < cut)  ilow++;
	  cut = med+NSig*Sig;
	  while (Vals[ihigh].val > cut)  ihigh--;
	  oldsig = Sig;
	  //DEBUG	  
	}
      double lowsum = 0;
      double highsum = 0;
      double sum2 = 0;
      double sum = 0;

      int lastk = ihigh;
      for (int k=ilow; k<= lastk; ++k)
	{
	  if (lowsum <= highsum)
	    {
	      lowsum += Vals[ilow].w;
	      ilow++;
	    }
	  else
	    {
	      highsum += Vals[ihigh].w;
	      ihigh--;
	    }
	  sum += Vals[k].w*Vals[k].val;
	  sum2 += Vals[k].w*sq(Vals[k].val);
	}
      // here ilow == ihigh+1
      ilow--;
      ihigh++;
      double vlow = Vals[ilow].val;
      double vhigh = Vals[ihigh].val;
      med = (0.5*(vlow+vhigh) + 
	     0.5*(highsum-lowsum)*(vhigh - vlow)
	     /((highsum<lowsum)? Vals[ilow].w : Vals[ihigh].w));
      sum /= (highsum+lowsum);
      Sig = sum2/(highsum+lowsum)-sq(sum);
      if (Sig >= 0) Sig = sqrt(Sig);
      else {Sig = -1;return 1e30;}
      if (iter>0 && fabs(oldsig-Sig)<0.001*Sig) break;
    }
  return med;



}

static double weighted_mean(const WeightedVals &Vals)
{
  if (Vals.size() == 0) return -1;
  double sumw = 0;
  double sum = 0;
  //  double sum2=0;
  WeightedVals::const_iterator p = Vals.begin();
  WeightedVals::const_iterator pstop = Vals.end();
  for ( ; p!= pstop; ++p)
    {
      sumw += p->w;
      sum += p->val*p->w;
    }
  return sum/sumw;
}
  

using namespace std;

int main(int nargs, char **args)
{
  if (nargs < 4)
    {
      cout << " usage :" << endl
	   << args[0] << " <rank1> <rank2> <aperse.list ...>" << endl
	   << "- You may use standalone_stars.list as input" << endl
	   << "- Output : average and sigma of (apfl<rank2>-apfl<rank1>)/apfl<rank1>"  << endl;
      return EXIT_FAILURE;
    }
  bool ok = true;
  unsigned irad1 = atoi(args[1]);
  unsigned irad2 = atoi(args[2]);  
  for (int k=3; k<nargs;++k)
    {
      string catName=args[k];
      try
	{
	  AperSEStarList cat(catName);
	  WeightedVals aperCorrs;
	  aperCorrs.reserve(cat.size());
	  for (AperSEStarCIterator i = cat.begin(); i != cat.end(); ++i)
	    {
	      const AperSEStar &s = **i;
	      if (s.Flag() != 0 
		  || s.FlagBad() != 0 ) continue;
	      // cuts to apply:
	      /* apfc6/apfl6<0.001.and.apfo6/apfl6<0.001.
		 and.apfl6/eapfl6>10.and.flagbad=0 */
	      const Aperture aper2 = s.apers[irad2];
	      const Aperture aper1 = s.apers[irad1];
	      if (aper1.fother > 0.001*aper1.flux || aper1.fcos >= 0.001*aper1.flux)
		continue;
	      double diff = aper2.flux-aper1.flux;
	      /* the minus sign in the next line is NOT a mistake. It is really
		 the uncertainty on the flux difference, because
		 cov(aper2.flux,aper1.flux) = Var(aper1.flux)
	      */
	      double diff_var = aper2.eflux*aper2.eflux - aper1.eflux*aper1.eflux;
	      if (diff_var <0 ) continue;
	      double val = diff/aper1.flux;
	      double val_var =  (diff_var+sq(val*aper1.eflux))/sq(aper1.flux);
	      aperCorrs.push_back(WeightedVal(val,1/val_var));
	    }
	  double sig;
	  double med = WeightedMedian(aperCorrs, sig, 3);
	  cout << catName << ' ' << med << ' ' << sig << endl;
	  ok = true;
	}
      catch (PolokaException p)
	{
	  p.PrintMessage(cout);
	  ok=false;
	}
    }
  return (ok) ?  EXIT_SUCCESS : EXIT_FAILURE;
}
