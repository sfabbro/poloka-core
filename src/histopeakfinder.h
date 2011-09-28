#ifndef HISTOPEAKFINDER__H
#define HISTOPEAKFINDER__H

#include <iostream> // for cout <<

#include <list>
#include "countedref.h"
#include "basestar.h"


struct StarScores : public RefCount
{
  double sx,sy;
  double nSig;
  double eventWeight;
  CountedRef <BaseStar> star;


  StarScores(const double Sc1, const double Sc2, const BaseStar *S):
       sx(Sc1), sy(Sc2), nSig(-1), eventWeight(1), star(S) {};
  StarScores(const double Sc1, const double Sc2, const double W,
	     const BaseStar *S):
       sx(Sc1), sy(Sc2), nSig(-1), eventWeight(W), star(S) {};

};


class StarScoresList : public std::list<CountedRef<StarScores> > 
{
  
};

typedef StarScoresList::const_iterator StarScoresCIterator;
typedef StarScoresList::iterator StarScoresIterator;




class Ellipse
{
 private:
  double xc,yc,wxx,wyy,wxy;
  
 public:

  Ellipse() { wxx=wyy=wxy = 0;}

  Ellipse(double  Xc, double Yc, double Wxx, double Wyy, double Wxy) :
    xc(Xc), yc(Yc), wxx(Wxx), wyy(Wyy), wxy(Wxy) {};

  void GetCenter(double &Xc, double &Yc) { Xc = xc; Yc = yc;}

  double Distance(const double &X, const double &Y) const;

  double SigmaX() const;
  
  double SigmaY() const;

  //! correlation coefficient
  double Corr() const;

  void dump(ostream &s = cout) const;

};

ostream & operator <<(ostream &s, const Ellipse &);





class Histo2d;
bool HistoPeakFinder(StarScoresList &List, const Histo2d &H,
		     const double &XGuess, const double &YGuess,
		     Ellipse &Ell, int verbose=0);

#endif /* HISTOPEAKFINDER__H */
