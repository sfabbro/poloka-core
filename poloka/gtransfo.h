// -*- C++ -*-
//
// \file gtransfo.h
// \brief Geometrical transformations (of 2D points)
// 
// 
#ifndef GTRANSFO_H
#define GTRANSFO_H

#include <iostream>
#include <string>
#include <vector>

#include <poloka/fatpoint.h>
#include <poloka/countedref.h>
#include <poloka/frame.h> // for whichTransformed


class StarMatchList;
class GtransfoLin;

//! a virtual (interface) class for geometric transformations. 
/*! We implement here One Gtransfo interface class, and actual derived
    classes. Composition in the usual (mathematical) sense is provided
    using GtransfoCompose(), and some classes (e.g. GtransfoLin)
    handle a * operator.  Generic inversion by iteration exists, but
    it is at least 10 times slower than the corresponding "direct
    transformation". If a transfo has an analytical inverse, then
    providing InverseTransfo is obviously a very good idea. Before
    resorting to InverseTransfo, consider using
    StarMatchList::InverseTransfo().  GtransfoLin::invert() and
    TanPix2RaDec::invert() exist.
    The classes also provide derivation and linear approximation.

*/



class Gtransfo: public RefCount{
public:
  
  //!
  virtual void  apply(const double Xin, const double Yin, 
		      double &Xout, double &Yout) const = 0 ;
  
  //! applies the tranfo to Pin and writes into Pout. Is indeed virtual.
  void apply(const Point &Pin, Point &Pout) const 
      {apply(Pin.x, Pin.y, Pout.x, Pout.y);}

  //!
  Point apply(const Point &Pin) const 
      {double xout, yout; apply(Pin.x, Pin.y, xout,yout); 
      return Point(xout,yout);}

  //! dumps the transfo coefficients to stream. 
  virtual void dump(ostream &stream = cout) const = 0;


  //! fits a transfo to a list of Point pairs (p1,p2, the Point fields in StarMatch).
  /*! After the fit this(p1) yields approximately p2. 
    The returned value is the sum of squared residuals.
    If you want to fit a partial transfo (e.g. such that 
    this(T1(p1)) = T2(p2), use StarMatchList::ApplyTransfo beforehand. */
  virtual double fit(const StarMatchList &List) = 0;

  //! allows to write MyTransfo(MyStar)
  void TransformStar(FatPoint &In) const { TransformPosAndErrors(In,In); }

  //! returns the local jacobian.
  virtual double Jacobian(const Point &P) const {return Jacobian(P.x, P.y);}

  //! returns a copy (allocated by new) of the transformation.
  virtual Gtransfo *Clone() const = 0;

  //! to be overloaded by derived classes if they can really "reduce" the composition (e.g. composition of Polynomial can be reduced)
  virtual Gtransfo *ReduceCompo(const Gtransfo *Right) const;


  //! returns the local jacobian.
  virtual double Jacobian(const double x, const double y) const;

  //! Computes the local Derivative of a transfo. Step is used for numerical derivation.
  virtual void Derivative(const Point &Where, GtransfoLin &Der, 
			  const double Step = 0.01) const;

  //! linear (local) approximation.
  virtual GtransfoLin LinearApproximation(const Point &Where, 
					  const double step = 0.01) const;

  virtual void TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const;


  //! transform errors (represented as double[3] in order V(xx),V(yy),Cov(xy))
  virtual void TransformErrors(const Point &Where, const double *VIn, 
			       double *VOut) const;

  //! returns an inverse transfo. Numerical if not overloaded.
  /*! Precision and Region refer to the "input" side of this, 
    and hence to the output side of the returned Gtransfo. */
  virtual Gtransfo* InverseTransfo(const double Precision,
				   const Frame& Region) const;


  //! Params should be at least Npar() long
  void GetParams(double *Params) const;

  //! 
  void SetParams(const double *Params);

  //!
  virtual double ParamRef(const int i) const;

  //!
  virtual double& ParamRef(const int i);

  //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar for y.
  virtual void ParamDerivatives(const Point &Where, double *Derivatives) const;

  //! Rough inverse. 
  /*! Stored by the numerical inverter to guess starting point 
     for the trials. Just here to enable overloading. */
  virtual Gtransfo* RoughInverse(const Frame &Region) const;

  //! returns the number of parameters (to compute chi2's)
  virtual int Npar() const {return 0;}

  void Write(const std::string &FileName) const;

  virtual void Write(ostream &stream) const;
  
  virtual ~Gtransfo() {};


};



typedef CountedRef<Gtransfo> GtransfoRef;

//! allows 'stream << Transfo;' (by calling T.dump(stream)). 
ostream & operator << (ostream &stream, const Gtransfo & T);


//! Returns a pointer to a composition. if Left->ReduceCompo(Right) return NULL, builds a GtransfoComposition and returns it. deletion of returned value to be done by caller

Gtransfo *GtransfoCompose(const Gtransfo *Left, const Gtransfo *Right);


/*=============================================================*/
//! A do-nothing transformation. It anyway has dummy routines to mimick a Gtransfo


class GtransfoIdentity : public Gtransfo {

public:
    //! constructor.
    GtransfoIdentity() {}

  //! Xout = Xin; Yout = Yin !
    void apply(const double Xin, const double Yin, 
	       double &Xout, double &Yout) const 
      {Xout = Xin; Yout = Yin;}; // to speed up

    double fit(const StarMatchList &List)
      {std:: cerr << "GtransfoIdentity cannot be fitted for list : "  
		  << &List << std::endl;
      return -1;
      }

    Gtransfo* ReduceCompo(const Gtransfo *Right) const { return Right->Clone();}
    void dump(ostream &stream = cout) const 
         { stream << "x' = x\ny' = y"<< endl;}

    int Npar() const {return 0;}
    Gtransfo *Clone() const { return new GtransfoIdentity;}

    void Derivative(const Point &Where, GtransfoLin &Derivative, 
		    const double Step = 0.01) const;

    //! linear approximation.
    virtual GtransfoLin LinearApproximation(const Point &Where, 
					    const double Step = 0.01) const;

  void Write(ostream &s) const;

  void Read(istream &s);


    //    ClassDef(GtransfoIdentity,1)
};

//! Shorthand test to tell if a transfo belongs to the GtransfoIdentity class. 
bool IsIdentity(const Gtransfo *a_transfo);

//! Shorthand test to tell if a transfo is a simple integer shift
bool IsIntegerShift(const Gtransfo *a_transfo);

/*====================   GtransfoPoly  =======================*/

//! Polynomial transformation class.
class GtransfoPoly : public Gtransfo
{						
private:
  unsigned deg; // the degree
  unsigned nterms; // number of parameters per coordinate
  vector<double> coeffs; // the actual coefficients
  vector<double> monomials; // value of monomials for the last call (should not go to IOs if any)

  /* use vector rather than double * to avoid 
     writing copy constructor and "operator =".
     Vect would work as well but introduces a dependence 
     that can be avoided */

public :
  //! Default transfo : identity for all degrees (>=1 ). The degree refers to the highest total power (x+y) of monomials. 
  GtransfoPoly(const unsigned Deg=1) ;

  //sets the polynomial degree. 
  void SetDegree(const unsigned Deg);

  void apply(const double Xin, const double Yin, 
	     double &Xout, double &Yout) const;

  Point apply(const Point &P) const
  { double xout, yout; apply(P.x, P.y, xout, yout); return Point(xout,yout);}

  //! specialised analytic routine
  void Derivative(const Point &Where, GtransfoLin &Der, 
		  const double Step = 0.01) const;

  //! a mix of apply and Derivative
  virtual void TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const;

  //! returns degree
  unsigned Degree() const { return deg;}


  //! total number of parameters
  int Npar() const { return 2*nterms;}

  //! print out of coefficients in a readable form.
  void dump(ostream &stream = cout) const;

  //! guess what
  double fit(const StarMatchList &List);

  //! Composition (internal stuff in quadruple precision)
  GtransfoPoly operator*(const GtransfoPoly &Right) const;

  Gtransfo *ReduceCompo(const Gtransfo *Right) const;

  Gtransfo *Clone() const {return new GtransfoPoly(*this);}


  //! access to coefficients (read only)
  double Coeff(const unsigned Powx, const unsigned Powy,
	       const unsigned WhichCoord) const;

  //! write access
  double& Coeff(const unsigned Powx, const unsigned Powy,
	       const unsigned WhichCoord);

  //! access to coefficients via variable names in "old" classes. Here to ease the transition
  double Coeff(const char *CoeffName) const;

  double& Coeff(const char *CoeffName);

  //! tells whether this name actually exists (given the degree)
  bool HasCoeff(const char* Name) const;

  double Determinant() const;

  //!
  double ParamRef(const int i) const;

  //!
  double& ParamRef(const int i);


  //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar for y.
  void ParamDerivatives(const Point &Where, double *Derivatives) const;


  void Write(ostream &s) const;
  void Read(istream &s);


private :
  double do_the_fit(const StarMatchList &List, const Gtransfo &InTransfo,
		    const bool UseErrors);

};

// approximates the inverse by a polynomial, up to required precision.
GtransfoPoly *InversePolyTransfo(const Gtransfo &Direct, const Frame &F, const double Prec);



/*=============================================================*/
//! implements the linear transformations (6 real coefficients). 
class GtransfoLin : public GtransfoPoly {
  
 public:

  //! the default constructor constructs the do-nothing transformation.
  GtransfoLin() : GtransfoPoly(1) {};


  //! This triggers an exception if P.Degree() != 1
  explicit GtransfoLin(const GtransfoPoly &P);

  //!  enables to combine linear tranformations: T1=T2*T3 is legal. 
  GtransfoLin  operator*(const  GtransfoLin &T2) const;
  
  //! returns the inverse: T1 = T2.invert(); 
  GtransfoLin  invert() const;
  
  
  // useful?    double Jacobian(const double x, const double y) const { return Determinant();}
  
  //!
  void Derivative(const Point &Where, GtransfoLin &Derivative, 
		  const double Step = 0.01) const;
  //!
  GtransfoLin LinearApproximation(const Point &Where, 
				  const double step = 0.01) const;
  
  
  //  void dump(ostream &stream = cout) const;
  
  // double fit(const StarMatchList &List);
    
  //! the constructor that enables to set all parameters independently. Not very useful. 
  GtransfoLin(const double ox, const double oy , const double aa11, 
	      const double aa12, const double aa21, const double aa22);
  
  //! Handy converter:
  GtransfoLin(const GtransfoIdentity &T) : GtransfoPoly(1) 
  { if (&T) {} /* avoid a warning */};
  
  Gtransfo* Clone() const { return new GtransfoLin(*this);}
  
  Gtransfo* InverseTransfo(const double Precision,
			   const Frame& Region) const;

  double A11() const { return Coeff(1,0,0);}
  double A12() const { return Coeff(0,1,0);}
  double A21() const { return Coeff(1,0,1);}
  double A22() const { return Coeff(0,1,1);}
  double Dx()  const { return Coeff(0,0,0);}
  double Dy()  const { return Coeff(0,0,1);}

protected :

  
  double& a11() { return Coeff(1,0,0);}
  double& a12() { return Coeff(0,1,0);}
  double& a21() { return Coeff(1,0,1);}
  double& a22() { return Coeff(0,1,1);}
  double& dx()  { return Coeff(0,0,0);}
  double& dy()  { return Coeff(0,0,1);}


  friend class Gtransfo;
  friend class GtransfoIdentity; // for Gtransfo::Derivative
  friend class GtransfoPoly; // // for Gtransfo::Derivative
  
protected:
  
};


/*=============================================================*/

//! just here to provide a specialized constructor, and fit.
class GtransfoLinShift : public GtransfoLin
{
public:
    //! Add ox and oy.
    GtransfoLinShift(double ox =0., double oy =0.) : GtransfoLin(ox,oy,1.,0.,0.,1.) {}
    GtransfoLinShift( const Point &P) : GtransfoLin(P.x, P.y, 1., 0. ,0. ,1.) {};
    double fit(const StarMatchList &List);

    int Npar() const {return 2;}
};

/*=============================================================*/
//! just here to provide a specialized constructor, and fit.
class GtransfoLinRot : public GtransfoLin {
  
 public: 
    GtransfoLinRot() : GtransfoLin() {};
    GtransfoLinRot(const double AngleRad, const Point *Center=NULL, 
		   const double ScaleFactor=1.0);
    double fit(const StarMatchList &List);

    int Npar() const {return 4;}
};


/*=============================================================*/

//! just here to provide specialized constructors. GtransfoLin fit routine.
class GtransfoLinScale :  public GtransfoLin {
  
 public: 
    //!
    GtransfoLinScale(const double Scale=1) : GtransfoLin(0.0, 0.0, Scale, 0.,0.,Scale) {};
    //!
    GtransfoLinScale(const double ScaleX, const double ScaleY) : 
	GtransfoLin(0.0, 0.0, ScaleX, 0.,0.,ScaleY) {};    

    int Npar() const {return 2;}
};


/*====================================================================*/

class TanRaDec2Pix; // the inverse of TanPix2RaDec.

//! the transformation that handles pix to sideral transfos (Gnomonic, possibly with polynomial distortions).
class TanPix2RaDec : public Gtransfo {
  

  GtransfoLin linPix2Tan; // pixels to tangent plane (internally in radians)
  GtransfoPoly *corr;
  double ra0, dec0; // in radians
  double cos0, sin0; // cos(dec0), sin(dec0) 
  
public:
    //! Pix2Tan describes the transfo from pix to tangent plane (in degrees). TangentPoint in degrees. Corrections are applied between Lin and deprojection parts (as in Swarp).
  TanPix2RaDec(const GtransfoLin &Pix2Tan, const Point &TangentPoint, 
	       const GtransfoPoly* Corrections = NULL);

  TanPix2RaDec(const TanPix2RaDec &Original);

#ifndef SWIG    
  void operator = (const TanPix2RaDec &);
#endif

    
  TanPix2RaDec();

  void apply(const double Xin, const double Yin, 
	     double &Xout, double &Yout) const;
  
  Point apply(const Point &Pin) const 
  {double xout, yout; apply(Pin.x, Pin.y, xout,yout); return Point(xout,yout);}

    //! composition with GtransfoLin
  TanPix2RaDec operator *(const GtransfoLin &Right) const;

  Gtransfo *ReduceCompo(const Gtransfo *Right) const;


    //! approximate inverse : it ignores corrections;
    TanRaDec2Pix invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    Gtransfo* RoughInverse(const Frame &Region) const;

    //! Inverse transfo: returns a TanRaDec2Pix if there are no corrections, or the iterative solver if there are.
    Gtransfo* InverseTransfo(const double Precision,
			     const Frame& Region) const;


    //! Sets the corrections (it can be a cubic ocrrection)
    void SetCorrections(const GtransfoPoly *Corrections);
    
    Gtransfo *Clone() const;

    void dump(ostream &stream) const;

    //! The tangent point (in degrees)
    Point TangentPoint() const;

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin LinPart() const;

    //! the "correction" 
    const GtransfoPoly* Corr() const; 

    //! transfo from pix to tangent plane ((*corr)*LinPart)
    GtransfoPoly Pix2TangentPlane() const;

    //! the CRPIX values (this is WCS jargon)
    Point CrPix() const;

    //!
    double fit(const StarMatchList &List);


    ~TanPix2RaDec();

};

//! This one is the Tangent Plane (called gnomonic) projection (from celestial sphere to tangent plane)
/*! this transfo does not implement corrections, since 
   they are defined the other way around (from pixels to sky), 
   and not invertible analytically. The inversion of tangent
   point WCS (TanPix2RaDec) is obtained via InverseTransfo().
*/

class TanRaDec2Pix : public Gtransfo 
{

    double ra0, dec0; //tangent point (internally in radians)
    double cos0,sin0;
    GtransfoLin linTan2Pix; // tangent plane to pixels (internally in radians)

 public:
    //! assume degrees everywhere.
    TanRaDec2Pix(const GtransfoLin &Tan2Pix, const Point &TangentPoint);
    
    //!
    TanRaDec2Pix();

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin LinPart() const;

    //! tangent point coordinates (in degrees)
    Point TangentPoint() const;

    //!
    void apply(const double Xin, const double Yin, double &Xout, double &Yout) const;

    //! transform with analytical derivatives
    void TransformPosAndErrors(const FatPoint &In, 
			       FatPoint &Out) const;



    //! exact typed inverse:
    TanPix2RaDec invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    Gtransfo* RoughInverse(const Frame &Region) const;

    //! Inverse transfo: returns a TanPix2RaDec.
    Gtransfo* InverseTransfo(const double Precision,
			   const Frame& Region) const;

    void dump(ostream &stream) const;

    Gtransfo * Clone() const; 

    double fit(const StarMatchList &List);


};


//! signature of the user-provided routine that actually does the coordinate transfo for UserTransfo.
typedef void (GtransfoFun)(const double, const double, 
			   double &, double &, const void*);


//! a run-time transfo that allows users to define a Gtransfo with minimal coding (just the transfo routine).
class UserTransfo : public Gtransfo
{
  private :
  
  GtransfoFun *userFun;
  const void *userData;

 public:
  //! the transfo routine and extra data that it may need.
  UserTransfo(GtransfoFun &Fun, const void *UserData);

  void apply(const double Xin, const double Yin, 
			  double &Xout, double &Yout) const;

  void dump(ostream &stream = cout) const;

  double fit(const StarMatchList &List);

  Gtransfo *Clone() const;

};


//! The virtual constructor from a file
Gtransfo* GtransfoRead(const std::string &FileName);
//! The virtual constructor from a file
Gtransfo* GtransfoRead(istream &s);


//! assumes that Transfo is a shift or involves a 'simple rotation'
Frame ApplyTransfo(const Frame& inputframe, const Gtransfo &T,
		   const WhichTransformed W = SmallFrame);


#endif /* GTRANSFO__H */
