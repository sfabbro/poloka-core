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

#include "point.h"
#include "countedref.h"

class StarMatchList;
class Frame;
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
  Point operator()(const Point &In) const { return apply(In);};

  //! allow composition of transformations regardless of their actual types.see GtransfoCompose() for a user callable entry.
  //! returns the local jacobian.
  virtual double Jacobian(const Point &P) const {return Jacobian(P.x, P.y);}

  virtual Gtransfo *ReduceCompo(const Gtransfo *Right) const;

  //! returns a copy (allocated by new) of the transformation.
  virtual Gtransfo *Clone() const = 0;

  //! returns the local jacobian.
  virtual double Jacobian(const double x, const double y) const;

  //! Computes the local Derivative of a transfo. Step is used for numerical derivation.
  virtual void Derivative(const Point &Where, GtransfoLin &Der, 
			  const double Step = 0.01) const;

  //! linear (local) approximation.
  virtual GtransfoLin LinearApproximation(const Point &Where, 
					  const double step = 0.01) const;

  //! transform errors (represented as double[3] in order V(xx),V(yy),Cov(xy))
  virtual void TransformErrors(const Point &Where, const double *VIn, 
			       double *VOut) const;

  //! returns an inverse transfo. 
  /*! Precision and Region refer to the "input" side of this, 
    and hence to the output side of the returned Gtransfo. */

  //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar for y.
  virtual void ParamDerivatives(const Point &Where, double *Derivatives) const;

  //! Params should be at least Npar() long
  void GetParams(double *Params) const;

  //! 
  void SetParams(const double *Params);

  //!
  virtual double ParamRef(const int i) const;

  //!
  virtual double& ParamRef(const int i);


  virtual Gtransfo* InverseTransfo(const double Precision,
				   const Frame& Region) const;

  //! Rough inverse. 
  /*! Stored by the numerical inverter to guess starting point 
     for the trials. Just here to enable overloading. */
#ifndef SWIG
  virtual Gtransfo* RoughInverse(const Frame &Region) const;

  //! returns the number of parameters (to compute chi2's)
  virtual int Npar() const {return 0;}
#endif /* SWIG */
  virtual ~Gtransfo() {};

};

typedef CountedRef<Gtransfo> GtransfoRef;

//! allows 'stream << Transfo;' (by calling T.dump(stream)). 
ostream & operator << (ostream &stream, const Gtransfo & T);


//! Returns a pointer to a composition. if Left->ReduceCompo(Right) return NULL, builds a GtransfoComposition and returns it. deletion of returned value to be done by caller

Gtransfo *GtransfoCompose(const Gtransfo *Left, const Gtransfo *Right);


/*=============================================================*/
//! A do-nothing transformation. It anyway has dummy routines to mimick a GTransfo


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
  
    //    ClassDef(GtransfoIdentity,1)
};

//! Shorthand test to tell if a transfo belongs to the GtransfoIdentity class. 
bool IsIdentity(const Gtransfo *a_transfo);

#ifndef SWIG
class GtransfoQuad;
class GtransfoCub;
#endif /* SWIG */

/*=============================================================*/
//! implements the linear transformations (6 real coefficients). 
class GtransfoLin : public Gtransfo {
  
private:

  
 public:

  //! the default constructor constructs the do-nothing transformation.
  GtransfoLin() {identity();}
  
  //!  enables to combine linear tranformations: T1=T2*T3 is legal. 
  GtransfoLin  operator*(const  GtransfoLin &T2) const;
  
  //! returns the inverse: T1 = T2.invert(); 
  GtransfoLin  invert() const;
  
  //!
  void  apply(const double Xin, const double Yin, double &Xout, double &Yout) const
  {
    Xout =  dx + a11*Xin + a12*Yin;
    Yout =  dy + a21*Xin + a22*Yin;
  }
  
  //! 
  double Determinant() const {return (a11*a22-a12*a21);}

  
  // useful?    double Jacobian(const double x, const double y) const { return Determinant();}
  
  //!
  void Derivative(const Point &Where, GtransfoLin &Derivative, 
		  const double Step = 0.01) const;
  //!
  GtransfoLin LinearApproximation(const Point &Where, 
				  const double step = 0.01) const;
  
  
  Point apply(const Point &Pin) const 
  { return Point(dx + a11*Pin.x + a12*Pin.y, dy + a21*Pin.x + a22*Pin.y);}
  
  void dump(ostream &stream = cout) const;
  
  double fit(const StarMatchList &List);
  
  
  //! the constructor that enables to set all parameters independently. Not very useful. 
  GtransfoLin(double ox, double oy , double aa11, double aa12, double aa21, double aa22) :
    dx(ox), dy(oy), a11(aa11), a12(aa12), a21(aa21), a22(aa22) {}
  
  //! Handy converter:
  GtransfoLin(const GtransfoIdentity &T)
  { if (&T) {} /* avoid a warning */ identity();}
  
  
  friend GtransfoCub operator*(const GtransfoLin &L, const GtransfoCub &R);
  
  friend GtransfoCub operator*(const GtransfoCub &L, const GtransfoLin &R);
  
  Gtransfo* Clone() const { return new GtransfoLin(*this);}
  
  Gtransfo* ReduceCompo(const Gtransfo *Right) const;

  Gtransfo* InverseTransfo(const double Precision,
			   const Frame& Region) const;

  double  A11() const { return a11;};
  double A12() const { return a12;};
  double A21() const { return a21;};
  double A22() const { return a22;};
  double dX()  const { return dx ;};
  double dY()  const { return dy ;};
  virtual int Npar() const {return 6;}
  
  double  ParamRef(const int i) const;
  double& ParamRef(const int i);
  void ParamDerivatives(const Point &Where, double *Derivatives) const;

  virtual int Degree() const {return 1;}
  
  friend class Gtransfo;
  friend class GtransfoIdentity; // for Gtransfo::Derivative
  
protected:
  double  dx,dy;
  double a11,a12,a21,a22;
  
  void identity() {dx=dy=a12=a21=0; a11=a22=1.;}
  
private:
  int do_the_fit(double &Chi2, const StarMatchList &List);

  // parameter serializer tools :
  struct LinParams { double GtransfoLin::*Item;};
  static LinParams Params[];
  
  //    ClassDef(GtransfoLin,1)
};


#ifndef SWIG
/*=============================================================*/
#endif /*SWIG */
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


/*=============================================================*/

//!implements the quadratic transformations (12 real coefficients).

/* Maybe it seems convenient to make it inherit from 
   GtransfoLin, but for quite a few methods (InverseTransfo,
   Derivative, LinearApproximation), GtransfoQuad
   should not resort to the GtransfoLin ones ! 
*/
class GtransfoQuad : public GtransfoLin {
  
public :
  //! the default constructor constructs the do-nothing transformation. 
    GtransfoQuad() {identity();}

  //! upgrade a linear transfo to a Quad one 
    GtransfoQuad(const GtransfoLin & Lin);

    //!
    GtransfoQuad operator*(const  GtransfoLin &R) const;

    void  apply(const double Xin, const double Yin, double &Xout, double &Yout) const
	{
	    Xout  =  x_tr(Xin,Yin);
	    Yout =  y_tr(Xin,Yin);
	}

    //void dump(ostream &stream = cout) const ;

    void dump(ostream &stream = cout) const;

    double fit(const StarMatchList &List);


    Gtransfo* Clone() const { return new GtransfoQuad(*this);}
    Gtransfo* ReduceCompo(const Gtransfo *Right) const;
    
    Gtransfo* InverseTransfo(const double Precision,
			     const Frame& Region) const;
    void Derivative(const Point &Where, GtransfoLin &Derivative, 
		    const double Step = 0.01) const;
    
    GtransfoLin LinearApproximation(const Point &Where, 
				    const double step = 0.01) const;
    
    double A1X2() const { return a1x2;};
    double A1XY() const { return a1xy;};
    double A1Y2() const { return a1y2;};
    double A2X2() const { return a2x2;};
    double A2XY() const { return a2xy;};
    double A2Y2() const { return a2y2;};
    int Npar() const {return 12;}
    // parameter serializer:
    double  ParamRef(const int i) const;
    double& ParamRef(const int i);
    void ParamDerivatives(const Point &Where, double *Derivatives) const;




    virtual int Degree() const {return 2;}

 protected:
    double a1x2,a1xy,a1y2,a2x2,a2xy,a2y2;

    void identity() {dx=dy=a12=a21=a1x2=a1xy=a1y2=a2x2=a2xy=a2y2=0; a11=a22=1.;}

 public: 

  /* the constructor that enables to set all parameters independently. Not very useful. */
    GtransfoQuad(double ox, double oy , double aa11, double aa12, double aa21, 
		 double aa22,double aa1x2, double aa1xy, double aa1y2,double aa2x2,
		 double aa2xy,double aa2y2) :
	GtransfoLin(ox,oy,aa11,aa12,aa21,aa22),
	a1x2(aa1x2), a1xy(aa1xy), a1y2(aa1y2), 
	a2x2(aa2x2), a2xy(aa2xy), a2y2(aa2y2) {}



    friend GtransfoQuad operator*(const GtransfoLin &L, const GtransfoQuad &R);

private:
    double x_tr(const double Xin, const double Yin) const 
	{ return  dx + a11*Xin + a12*Yin + a1x2*Xin*Xin + a1xy*Xin*Yin + a1y2*Yin*Yin;}

    struct QuadParams { double GtransfoQuad::*Item;};
    static QuadParams Params[];



  double y_tr(const double Xin, const double Yin) const 
  { return  dy + a21*Xin + a22*Yin + a2x2*Xin*Xin + a2xy*Xin*Yin + a2y2*Yin*Yin;}
    GtransfoQuad truncated_product(const GtransfoQuad &R) const;  
    //    ClassDef(GtransfoQuad,1);
};

GtransfoQuad operator*(const GtransfoLin &L, const GtransfoQuad &R);

/*=============================================================*/
//! very simple stuff to associate names and values. Used to I/O transfos to fits headers.
struct NamedValue {
  string name;
  double value;
  NamedValue(const string &a_name, const double a_value) : name(a_name), value(a_value) {};
};


#include <vector>



/* GtransfoCub*/

//! implements the cubic transformations (20 real coefficients).
class GtransfoCub : public GtransfoQuad {
  
 public:
    //! the default constructor constructs the do-nothing transformation.
    GtransfoCub() {identity();}
    
    //! upgrade a linear transfo to a cubic one 
    GtransfoCub(const GtransfoLin & Lin);
    
    //! upgrade a quadratic transfo to a cubic one 
    GtransfoCub(const GtransfoQuad & Quad);

    //! Cub*Lin
    friend GtransfoCub operator*(const GtransfoCub &L, const GtransfoLin &R);

    //! Lin*Cub
    friend GtransfoCub operator*(const GtransfoLin &L, const GtransfoCub &R);

    void  apply(const double Xin, const double Yin, double &Xout, double &Yout) const
	{
	    Xout =  x_tr3(Xin,Yin);
	    Yout =  y_tr3(Xin,Yin);
	}

    Point apply(const Point &Pin)
	{ return Point(x_tr3(Pin.x,Pin.y), y_tr3(Pin.x,Pin.y));}

    void dump(ostream &stream = cout) const;

    /*! fits a transfo to a list of star pairs (p1,p2). After the fit
       this(PriorTransfo(p1)) yields approximately p2. The returned value is the chi2.*/
    double fit(const StarMatchList &List);

   //!
    Gtransfo* Clone() const { return new GtransfoCub(*this);}
    Gtransfo* ReduceCompo(const Gtransfo *Right) const;


    
    double dX() const { return dx ;};
    double dY() const { return dy ;}; 

    double A11() const { return a11;};
    double A12() const { return a12;};
    double A21() const { return a21;};
    double A22() const { return a22;};

    double A1X2() const { return a1x2;};
    double A1XY() const { return a1xy;};
    double A1Y2() const { return a1y2;};
    double A2X2() const { return a2x2;};
    double A2XY() const { return a2xy;};
    double A2Y2() const { return a2y2;};

    double A1X3() const { return a1x3;};
    double A1X2Y() const { return a1x2y;};
    double A1XY2() const { return a1xy2;};
    double A1Y3() const { return a1y3;};
    double A2X3() const { return a2x3;};
    double A2X2Y() const { return a2x2y;};
    double A2XY2() const { return a2xy2;};
    double A2Y3() const { return a2y3;};

    typedef enum {Old, New} OldOrNew;
    void GetValues(vector<NamedValue> &Values, const string &KeyHeader = "DV") const;
    void SetValues(vector<NamedValue> &Values, const string &KeyHeader);
    int Npar() const {return 20;}

    virtual int Degree() const {return 3;}
  
     // parameter serializer
    double  ParamRef(const int i) const;
    double& ParamRef(const int i);
    void ParamDerivatives(const Point &Where, double *Derivatives) const;
  

protected:
    //    double dx,dy; //2
    //    double a11,a12,a21,a22; //4 
    //    double a1x2,a1xy,a1y2,a2x2,a2xy,a2y2; //6
    double a1x3,a1x2y,a1xy2,a1y3,a2x3,a2x2y,a2xy2,a2y3; //8

    void identity() {dx=dy=a12=a21=a1x2=a1xy=a1y2=a2x2=a2xy=a2y2=
			 a1x3=a1x2y=a1xy2=a1y3=a2x3=a2x2y=a2xy2=a2y3=0; a11=a22=1.;}


private:


  // arrays used to associate names to class data members once for all.
    typedef double GtransfoCub::*GtransfoCubItem;
    typedef struct { const char *name; GtransfoCubItem value; }  CubAssoc;
    static CubAssoc WCS3Names[];
    static CubAssoc DJNames [];
    static CubAssoc DVNames [];



    double x_tr3(const double Xin, const double Yin) const 
	{ return  dx + a11*Xin + a12*Yin + a1x2*Xin*Xin + a1xy*Xin*Yin + a1y2*Yin*Yin +
	      a1x3*Xin*Xin*Xin + a1x2y*Xin*Xin*Yin + a1xy2*Xin*Yin*Yin + 
	      a1y3*Yin*Yin*Yin;}

    struct CubParams { double GtransfoCub::*Item;};
    static CubParams Params[];


    double y_tr3(const double Xin, const double Yin) const 
	{ return  dy + a21*Xin + a22*Yin + a2x2*Xin*Xin + a2xy*Xin*Yin + a2y2*Yin*Yin +
	      a2x3*Xin*Xin*Xin + a2x2y*Xin*Xin*Yin + a2xy2*Xin*Yin*Yin + 
	      a2y3*Yin*Yin*Yin;} 

    //  ClassDef(GtransfoCub,1);

};

/*====================================================================*/

class TanRaDec2Pix; // the inverse of TanPix2RaDec.

//! the transformation that handles pix to sideral transfos (Gnomonic, possibly with polynomial distortions).
class TanPix2RaDec : public Gtransfo {
  

    GtransfoLin linPix2Tan; // pixels to tangent plane (internally in radians)
    GtransfoQuad *corr;
    double ra0, dec0; // in radians
    double cos0, sin0; // cos(dec0), sin(dec0) 
  
 public:
    //! Pix2Tan describes the transfo from pix to tangent plane (in degrees). TangentPoint in degrees. Corrections are applied between Lin and deprojection parts (as in Swarp).
    TanPix2RaDec(const GtransfoLin &Pix2Tan, const Point &TangentPoint, 
		 const GtransfoQuad* Corrections = NULL);

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
    void SetCorrections(const GtransfoQuad *Corrections);
    
    Gtransfo *Clone() const;

    void dump(ostream &stream) const;

    //! The tangent point (in degrees)
    Point TangentPoint() const;

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin LinPart() const;

    //! the correction (can be more than quadratic) 
    const GtransfoQuad* Corr() const; 

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
    void apply(const double Xin, const double Yin, double &Yout, double &Yout) const;

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

//! probably obsolete. use LinearApproximation instead
GtransfoLin *GtransfoToLin(const Gtransfo* transfo);


//#include "gtransfocomposition.h"


#endif /* GTRANSFO__H */
