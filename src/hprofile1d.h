#ifndef HPROFILE1D__H
#define HPROFILE1D__H


/*! a class to create profile histograms using medians for
  the bin averages 
*/

#include <list>
#include <vector>
#include <string>

using namespace std;


struct WeightedVal
{
  double val,w;
  WeightedVal(const double Val, const double W) : val(Val), w(W) {};
  friend bool operator < (const WeightedVal &w1, const WeightedVal &w2)
  { return (w1.val<w2.val);}
};





class HProfile1D
{
  public :
    HProfile1D(const int NBins, const double &XMin, const double &XMax, const
	     double &YMin = -1e30, const double &Ymax = 1e30);


  struct Bin
  {
    bool binOK;
    double median;
    double error;
    int count;
    double sumw;
    list<WeightedVal> values;
    Bin () : binOK(false), median(1e30), error(-1), count(0), sumw(0) {};
    void AddVal(const double &Val, const double &W);
    void Stat(double &Med, double &Err, int &Count);
  };

  vector<Bin> bins;
  int nbins;
  double xMin, xMax, scalex;
  double yMin, yMax;

  int NBins() const { return nbins;}
  
  double BinCenter(int IBin) const;

  void BinStat(const int IBin, 
	       double &X, double &Mean, double &Err, int &Count);

  void Fill(const double &X, const double &Y, const double &W = 1);

  void write(const string &FileName, const char* XTag=NULL, const char *YTag=NULL);

};
  




#endif /* HPROFILE1D__H */
