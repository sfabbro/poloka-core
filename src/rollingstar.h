#ifndef ROLLINGSTAR__H
#define ROLLINGSTAR__H

#include "basestar.h"
#include "image.h"
#include "nstarmatch.h"
#include "sestar.h"
#include "reducedimage.h"

class RollingStar : public BaseStar {
  //class RollingStar {
 private:
  double Faper_ref;
  int Nlists;
 public:
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


typedef list<RollingStar*>::const_iterator RollingStarCIterator;
typedef list<RollingStar*>::iterator RollingStarIterator;

//#include "starlist.h"
class RollingStarList : public list<RollingStar*> 
{
 private:
  int Nlists;
 public:

  void  write(const string &FileName) ;
  

  RollingStarList(const NStarMatchList & nstml, const ReducedImage &ImRef);
  
  ~RollingStarList()
    { for (RollingStarIterator s=begin(); s!=end(); ++s)
      {delete *s; s=erase(s);}
    }

};




//BaseStarList * Rolling2Base(RollingStarList * This);
//const BaseStarList* Rolling2Base(const RollingStarList * This);

#endif /* ROLLINGSTAR__H */
