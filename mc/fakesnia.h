#ifndef FAKESNIA__H
#define FAKESNIA__H


#include <iostream>
//#include "basestar.h"
//#include "sestar.h"
#include "point.h"
#include "countedref.h"
#include "cosmology.h"
#include "saltmodel.h"
//#include "fakeobject.h"
#include "instrument.h"
#include <cstdio>

using namespace std;

//using saltmodel

int DecodeFormat(const char *FormatLine, const char *StarName);

class FakeSNIa : public Point, public RefCount

{

protected:
// sn data
  double redshift;
  double d0; // MJDday of the max mag
  double stretch;
  double alpha; // correction de strectch  
  double color;
  double beta; // correction de couleur
  double dispersion;
  double MW_rv; // host extibction rv
  double MW_extinction; //  
// host data
  double h_extinction;
  double h_rv;
  double h_ra;
  double h_dec;
  double h_imag;



 private:
 void Set_to_Zero();

 public:
 
 double Ra() const {return x;}
 double& Ra()  {return x;}
 
 double Dec() const {return y;}
 double& Dec()  {return y;}
 
 double Redshift() const {return redshift;}
 double& Redshift(){return redshift;}
 
 double D0() const {return d0;}
 double& D0(){return d0;}
 
 double Stretch() const {return stretch;}
 double& Stretch(){return stretch;}
 
 double MW_Extinction() const {return MW_extinction;}
 double& MW_Extinction(){return MW_extinction;}
 
 double H_Extinction() const {return h_extinction;}
 double& H_Extinction(){return h_extinction;}
 

 double Color()const {return color;}
 double& Color(){return color;}

 double Dispersion()const {return dispersion;}
 double& Dispersion(){return dispersion;}

 double Alpha()const {return alpha;}
 double& Alpha(){return alpha;}
 
 double Beta()const {return beta;}
 double& Beta(){return beta;}
 
 double H_Rv()const {return h_rv;}
 double& H_Rv(){return h_rv;}
 
 double MW_Rv()const {return MW_rv;}
 double& MW_Rv(){return MW_rv;}
 
   
 double H_Ra()const {return h_ra;}
 double& H_Ra(){return h_ra;}
 
 double H_Dec()const {return h_dec;}
 double& H_Dec(){return h_dec;} 
 
 double H_Imag()const {return h_imag;}
 double& H_Imag(){return h_imag;}

 
 
 double Magnitude(const GeneralCosmo &cosmo_model, Instrument &instrum, const
 string &filter,const double day, const string magsys, const double H0);
// no random
FakeSNIa(double xx, double yy, double RS,double day,double istretch,double ih_extinction, double icolor, 
	 double idispersion, double ialpha, double ibeta, double ih_rv, double iMW_rv, double iMW_extinction);

 FakeSNIa(double xx, double yy);
 FakeSNIa();
  /* DOC \noindent {\bf For read & write}: */
  
  //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;

  //! for dump
  virtual void    dump(ostream& s = cout) const ;
  
  virtual void write(ostream& s = cout) const ;
  
  //! for write with NO end-of-line
  virtual void    writen(ostream& s = cout) const ;

  //! to read once the object is created 
  void   read_it(istream& r, const char *Format); 

   //! to read and create the object  

  static FakeSNIa* read(istream& r, const char *Format); 

  /* DOCF  to write the FakeList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */

  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  
  virtual void WriteHeader (ostream & stream =cout) const; 

  static const char *TypeName() { return "FakeSNIa";}

};


#include "fakelist.h"
#include <list>

#ifdef USE_ROOT

typedef FakeListWithRoot<FakeSNIa> FakeSNIaList;

#else

typedef FakeList<FakeSNIa> FakeSNIaList;

#endif /* USE_ROOT */



typedef FakeSNIaList::const_iterator FakeSNIaCIterator;
typedef FakeSNIaList::iterator FakeSNIaIterator;
typedef CountedRef<FakeSNIa> FakeSNIaRef;

#endif
