#ifndef ROLLINGSTAR__H
#define ROLLINGSTAR__H

#include "persistence.h"

#include "basestar.h"
#include "image.h"
#include "nstarmatch.h"
#include "sestar.h"
#include "reducedimage.h"

class RollingStar : public BaseStar {
  CLASS_VERSION(RollingStar,1);
  #define RollingStar__is__persistent
  //class RollingStar {
 private:
  double Faper_ref;
  int Nlists;
 public:
  RollingStar() {}
  RollingStar(int & n);

  
  //Concern the candidate

  // Concern the possible host
  BaseStar host;
  double dass_ref;


  //Concern every point on the lightcurve, should add the date.
  double *flux;
  double *eflux;
  double *stnoise;
  double *prcti;
  double *chi2;
  //~RollingStar(){};
  
  
  //double Faper_ref()  const {return faper_ref;};
  
#ifdef STORAGE

  /* DOCF for dump with NO end-of-line*/
  virtual void dumpn(ostream& s = cout) const;

  /* DOCF for dump  */
  virtual void dump(ostream& s = cout) const ;
  
  /* DOCF for write with NO end-of-line*/
  virtual void writen(ostream& s = cout) const ;
  
  /* DOCF to read once the object is created */
  virtual void    Read(istream& r, const char *Format); 
  
  /* DOCF to read and create the object  */
  
  static RollingStar* read(istream& r, const char* Format); 
#endif
  /* DOCF  to write the StarList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */
  void  WriteHeader(ostream & pr ) ;
  
  void write(ostream &pr) const;

  
  /* DOCF to write the StarList header */
  static const char *TypeName() { return "RollingStar";}
  

  
};




#include "starlist.h"
//class RollingStarList : public list<RollingStar*> 
class RollingStarList : public StarList<RollingStar> 
 
{ 
  CLASS_VERSION(RollingStarList,1);
  #define RollingStarList__is__persistent
 private:
  int Nlists;
 public:

  void  write(const string &FileName) ;
  
  RollingStarList(){Nlists=0;}; // for persistence
  RollingStarList(const NStarMatchList & nstml, const ReducedImage &ImRef);
  
/*   ~RollingStarList() */
/*     { for (RollingStarIterator s=begin(); s!=end(); ++s) */
/*       {delete *s; s=erase(s);} */
/*     } */

};

typedef RollingStarList::const_iterator RollingStarCIterator;
typedef RollingStarList::iterator RollingStarIterator;


//BaseStarList * Rolling2Base(RollingStarList * This);
//const BaseStarList* Rolling2Base(const RollingStarList * This);

#endif /* ROLLINGSTAR__H */
