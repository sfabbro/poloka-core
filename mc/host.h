#ifndef HOST__H
#define HOST__H

#include <iostream>
#include "point.h"
#include "countedref.h"
#include <cstdio>

using namespace std;


class Host : public Point , public RefCount 

#ifdef USE_ROOT
	, public TObject
#endif
{
protected:
double z;
double err_z;
double gal_type;



double a;
double b;
double theta;
double i_mag;

 public:
 
 double Ra() const {return x;}
 double& Ra()  {return x;}
 double Dec() const {return y;}
 double& Dec()  {return y;}
 double Z() const {return z;}
 double& Z()  {return z;}
 double Err_Z() const {return err_z;}
 double& Err_Z()  {return err_z;} 
 double Gal_Type() const {return gal_type;}
 double& Gal_Type()  {return gal_type;}
 
 double A() const {return a;}
 double& A()  {return a;}
 double B() const {return b;}
 double& B()  {return b;} 
 double Theta() const {return theta;}
 double& Theta()  {return theta;}
 double I_Mag() const {return i_mag;}
 double& I_Mag()  {return i_mag;}
 
Host(double Ra, double Dec);
Host();

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

  static Host* read(istream& r, const char *Format); 

  /* DOCF  to write the FakeList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */

  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  
  virtual void WriteHeader (ostream & stream =cout) const; 

  static const char *TypeName() { return "Host";}

};


#include "fakelist.h"
#include <list>

#ifdef USE_ROOT

typedef FakeListWithRoot<Host> HostList;

#else

typedef FakeList<Host> HostList;

#endif /* USE_ROOT */



typedef HostList::const_iterator HostCIterator;
typedef HostList::iterator HostIterator;
typedef CountedRef<Host> HostRef;

#endif // HOST__H
