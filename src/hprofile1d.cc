#include "hprofile1d.h"

#include <cmath>
#include <iostream>
#include <algorithm>


HProfile1D::HProfile1D(const int NBins, 
		       const double &XMin, const double &XMax, 
		       const double &YMin, const double &YMax)
{
  xMin = XMin;
  xMax = XMax;
  yMin = YMin;
  yMax = YMax;
  nbins = 0;
  if (xMin == xMax)
    {
      cout << " cannot handle a HProfile1D with XMin == XMax !!! " << endl;
      return;
    }
  nbins = NBins;
  if (xMax < xMin) swap(xMax, xMin);
  scalex = nbins/(xMax-xMin);
  bins.resize(nbins);
}

void HProfile1D::Bin::AddVal(const double &Val, const double &W)
{
  binOK = false;
  values.push_back(WeightedVal(Val,W));
  sumw += W;
  count++;
}


static double weighted_median(list<WeightedVal> &values, const double sumw )
{
  values.sort();
  double halfSumWeight = sumw*0.5;
  double runningSum = 0;
  list<WeightedVal>::iterator  i = values.begin(); 
  int count = 0;
  for (;i != values.end(); ++i)
    {
      runningSum += i->w;
      if (runningSum> halfSumWeight) break;
      count++;
    }
  if (count == 0) return 0;
  if (count == 1) return i->val;
  list<WeightedVal>::iterator other = i;
  if (i==values.begin()) other++; else other--;
  // not sure about this ....
  return (other->val - i->val)*runningSum/sumw + i->val;
}  

static double sq(const double &x) { return x*x;}


void HProfile1D::Bin::Stat(double &Med, double &Err, int &Count)
{
  if (!binOK)
    {
      median = weighted_median(values, sumw);
      list<WeightedVal> squares;
      for (list<WeightedVal>::const_iterator i=values.begin(); 
	   i != values.end(); ++i)
	squares.push_back(WeightedVal(sq(i->val-median),i->w));
      double med_squares = weighted_median(squares,sumw);
      /* this 2.19 was determined experimentally to get the right
	 answer for a gaussian distribution */
      error = sqrt(med_squares*2.19);
      binOK = true;
    }
  Med = median;
  Err = error;  
  Count = count;
}

double  HProfile1D::BinCenter(int IBin) const
{
  if (IBin <0 || IBin >= nbins)
    {
      cerr << "HProfile1D::BinCenter : bin index out of range " << endl;
      return -1e30;
    }
  return (xMin + (IBin+0.5)/scalex);
}

void  HProfile1D::BinStat(const int IBin, 
			  double &X, double &Mean, double &Err, int &Count)
{
  if (IBin <0 || IBin >= nbins)
    {
      cerr << "HProfile1D::BinStat : bin index out of range " << endl;
      Mean =0;
      Err = -1;
      Count = 0;
      X = -1e30;
      return;
    }
  X = BinCenter(IBin);
  bins[IBin].Stat(Mean,Err,Count);
}
  
/* for a gaussian, med(x**2)-med(x)**2 = 1.54*var */


void HProfile1D::Fill(const double &X, const double &Y, const double &W)
{
  int bin = int(floor((X - xMin)*scalex)); 
  if (bin<0 || bin>=nbins) return; 
  if (Y<yMin || Y>yMax) return;
  bins[bin].AddVal(Y,W);
}


#include <fstream>
void HProfile1D::write(const string &FileName, 
		       const char* XTag, const char *YTag)
{
  ofstream s(FileName.c_str());
  string xtag((XTag)? XTag : "x");
  string ytag((YTag)? YTag : "y");
  s << "# " << xtag  << " : " << endl;
  s << "# " << ytag  << " : " << endl;
  s << "# e" << ytag << " : spread" << endl;
  s << "# c" << ytag << " : event count" << endl;
  s << "#end " << endl;
  int nbin = NBins();
  for (int k=0; k<nbin; ++k)
    {
      double x,mean,err;
      int count;
      BinStat(k, x,mean,err,count);
      s << x << ' ' << mean << ' ' << err << ' ' << count 
	<< endl;
    }
  s.close();
}
