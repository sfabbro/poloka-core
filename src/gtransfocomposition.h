#ifndef GTRANSFOCOMPOSITION__H
#define GTRANSFOCOMPOSITION__H

#include "gtransfo.h"

class GtransfoComposition : public Gtransfo {
  private :
    Gtransfo* first, *second;
  public :
    GtransfoComposition(const Gtransfo *Second, const Gtransfo *First);
    void apply(const double Xin, const double Yin, double &Xout, double &Yout) const;
    void dump(ostream &stream = cout) const; 
    double fit(const StarMatchList &List, const Gtransfo *PriorTransfo = NULL , const Gtransfo *PosteriorTransfo = NULL);
    Gtransfo *Clone() const;
    ~GtransfoComposition();

    //#ifndef SWIG
    //  ClassDef(GtransfoComposition,1);
    //#endif /*SWIG */
};


#endif
