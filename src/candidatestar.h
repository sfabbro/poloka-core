// This may look like C code, but it is really -*- C++ -*-


#ifndef  CANDISTAR__H 
#define  CANDISTAR__H


#include "sestar.h"
#include "image.h"
#include "frame.h"


class CandidateStar : public SEStar {

public:
  CandidateStar();
  CandidateStar(SEStar const & sestar);
  
  
  int Numero() const {return numero;};
  double FluxCv()  const {return fluxcv;};
  double Noise()  const {return noise;};
  double SigToNoise()  const {return sigtonoise;};
  double Faper_ref()  const {return faper_ref;};
  double PrctIncrease()  const {return prctincrease;} ;
  double Sigx() const {return sigx;};
  double Sigy() const {return sigy;};
  double Sigxy() const {return sigxy;};
  double Chi2() const {return chi2;};

  int& Numero()  {return numero;};
  double& FluxCv() {return fluxcv;};
  double& Noise()  {return noise;};
  double& SigToNoise()  {return sigtonoise;};
  double& Faper_ref()   {return faper_ref;};
  double& PrctIncrease()   {return prctincrease;} ; 
  double& Sigx() {return sigx;};
  double& Sigy() {return sigy;};
  double& Sigxy() {return sigxy;};
  double& Chi2() {return chi2;};
  
  bool  IsGood(double cut_sigtonoise);
  bool  IsGood(double cut_sigtonoise, double frac_faper_ref);
  void ComputeNoise(Image & imgref, Image & imgnew, double rayon_ref, double rayon_new,
		    double sigma_ref, double sigma_new, double fd_ref, double fd_new);
  void  ComputeFaperRef(Image & imgref, double rayon,double  fd_ref);
  
  
  void write_scan(ostream & pr) const;
  static void WriteHeader_Scan(ostream &pr);
  
  /* DOCF for dump with NO end-of-line*/
  virtual void dumpn(ostream& s = cout) const;
  /* DOCF for dump  */
  virtual void dump(ostream& s = cout) const ;
  
  /* DOCF for write with NO end-of-line*/
  virtual void writen(ostream& s = cout) const ;
  
  /* DOCF to read once the object is created */
  virtual void    read_it(istream& r, const char *Format); 
  
  /* DOCF to read and create the object  */
  
  static CandidateStar* read(istream& r, const char* Format); 
  
  /* DOCF  to write the StarList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */
  string WriteHeader_(ostream & pr = cout, const char*i = NULL) const;
  
 
  
  /* DOCF to write the StarList header */
  static const char *TypeName() { return "CandidateStar";}
  
private:
  // uniquement appelee par le constructeur, donc private  
  void Set_to_Zero();
  
protected:
  int numero;
  double fluxcv ;
  double noise ;
  double sigtonoise;
  double faper_ref;
  double prctincrease; 
  double sigx;  
  double sigy;
  double sigxy;
  double chi2;
};

/* what concerns the CandidateStarList's : */
#include "starlist.h"

typedef StarList<CandidateStar> CandidateStarList;
typedef StarList<CandidateStar>::const_iterator CandidateStarCIterator;
typedef StarList<CandidateStar>::iterator CandidateStarIterator;

BaseStarList* Candidate2Base(CandidateStarList * This);
const BaseStarList* Candidate2Base(const CandidateStarList * This);

#endif
