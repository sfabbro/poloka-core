#include <iostream>
#include <iomanip>
#include <iterator> /* for ostream_iterator */
#include <cmath> // for sin and cos and may be others
#include <fstream>
#include <sstream>
#include <cassert>

#include <poloka/starmatch.h>
#include <poloka/matvect.h>
#include <poloka/gtransfo.h>
#include <poloka/frame.h>
#include <poloka/polokaexception.h>


bool IsIdentity(const Gtransfo *a_transfo)
{ return (dynamic_cast<const GtransfoIdentity*>(a_transfo) != NULL);}

bool IsIntegerShift(const Gtransfo *a_transfo)
{
  const GtransfoPoly* shift = dynamic_cast<const GtransfoPoly*>(a_transfo);
  if (shift == NULL) return false;

  static const double eps = 1e-5;

  double dx = shift->Coeff(0,0,0);
  double dy = shift->Coeff(0,0,1);
  double a11 = floor(shift->Coeff(1,0,0)+0.5);
  double a22 = floor(shift->Coeff(0,1,1)+0.5);

  static Point dumb(4000, 3000);
  if (fabs(dx - int(floor(dx+0.5))) < eps && 
      fabs(dy - int(floor(dy+0.5))) < eps &&
      fabs(a11*dumb.x+dx - shift->apply(dumb).x) < eps &&
      fabs(a22*dumb.y+dy - shift->apply(dumb).y) < eps)    
    return true;

  return false;
}

/********* Gtransfo ***********************/

Gtransfo* Gtransfo::ReduceCompo(const Gtransfo *Right) const 
{// by default no way to compose
  if (Right) {} // avoid a warning
  return NULL;
}

double Gtransfo::Jacobian(const double X, const double Y) const
{
double x2,y2;
double eps=X*0.01;
if (eps == 0) eps = 0.01;
apply(X,Y, x2, y2);
double dxdx, dydx;
apply(X+eps, Y, dxdx, dydx);
dxdx -= x2; dydx -= y2;
double dxdy, dydy;
apply(X, Y+eps, dxdy, dydy);
dxdy -= x2; dydy -= y2;
return ((dxdx * dydy - dxdy * dydx)/(eps*eps));
}

/*! the Derivative is represented by a GtransfoLin, in which
  (hopefully), the offset terms are zero. Derivative should 
  transform a vector of offsets into a vector of offsets. */
void Gtransfo::Derivative(const Point &Where,
			  GtransfoLin &Der, const double Step) const
{
  double x = Where.x;
  double y = Where.y;
  double xp0,yp0;
  apply(x, y, xp0, yp0);

  double xp,yp;
  apply(x+Step, y, xp, yp);
  Der.a11() = (xp-xp0)/Step;
  Der.a21() = (yp-yp0)/Step;
  apply(x, y + Step, xp, yp);
  Der.a12() = (xp-xp0)/Step;
  Der.a22() = (yp-yp0)/Step;
  Der.dx() = 0;
  Der.dy() = 0;
}  

GtransfoLin Gtransfo::LinearApproximation(const Point &Where, 
					  const double Step) const
{
  Point outWhere = apply(Where);
  GtransfoLin der;
  Derivative(Where, der, Step);
  return GtransfoLinShift(outWhere.x, outWhere.y)*der*GtransfoLinShift(-Where.x, -Where.y);
}


void Gtransfo::TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const
{
  FatPoint res; // in case In and Out are the same address...
  res = apply(In);
  GtransfoLin der;
  // could save a call here, since Derivative needs the transform of where that we already have
  // 0.01 may not be a very good idea in all cases. May be we should provide a way of altering that.
  Derivative(In, der, 0.01);
  double a11 = der.A11();
  double a22 = der.A22();
  double a21 = der.A21();
  double a12 = der.A12();
  res.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  res.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  res.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;
  Out = res;
}


void Gtransfo::TransformErrors(const Point &Where,
			       const double *VIn, double *VOut) const
{
  GtransfoLin der;
  Derivative(Where, der, 0.01);
  double a11 = der.A11();
  double a22 = der.A22();
  double a21 = der.A21();
  double a12 = der.A12();

  /*      (a11 a12)          (vxx  vxy)
     M =  (       )  and V = (        )
          (a21 a22)          (xvy  vyy)
     
     Vxx = Vin[0], vyy = Vin[1], Vxy = Vin[2];
     we want to compute M*V*tp(M)
     A lin alg light package would be perfect...
  */
  int xx = 0;
  int yy = 1;
  int xy = 2;
  // M*V :
  
  double b11 = a11 * VIn[xx] + a12 * VIn[xy];
  double b22 = a21 * VIn[xy] + a22 * VIn[yy];
  double b12 = a11 * VIn[xy] + a12 * VIn[yy];
  double b21 = a21 * VIn[xx] + a22 * VIn[xy];

  // (M*V) * tp(M)

  VOut[xx] = b11 * a11 + b12 * a12;
  VOut[xy] = b11 * a21 + b12 * a22;
  VOut[yy] = b21 * a21 + b22 * a22;
}

Gtransfo* Gtransfo::RoughInverse(const Frame &Region) const
{
  // "in" and "out" refer to the inverse direction.
  Point centerOut = Region.Center();
  Point centerIn = apply(centerOut);
  GtransfoLin der;
  Derivative(centerOut,der,sqrt(Region.Area())/5.);
  der = der.invert();
  der = GtransfoLinShift(centerOut.x, centerOut.y)
    *der
    *GtransfoLinShift(-centerIn.x, -centerIn.y);
  return new GtransfoLin(der);
}


/* implement one in Gtransfo, so that all derived 
   classes do not need to provide one... */


/* the routines that follow are used for ea generic parameter 
   transformation serialization, used e.g. for fits. Enables
   to manipulate transformation parameters as vectors.
*/


// not dummy : what it does is virtual because ParamRef is virtual.
void Gtransfo::GetParams(double *Params) const
{
  int npar = Npar();
  for (int i=0; i<npar ; ++i) Params[i] = ParamRef(i);
}


void Gtransfo::SetParams(const double *Params)
{
  int npar = Npar();
  for (int i=0; i<npar ; ++i) ParamRef(i) = Params[i];
}  

double Gtransfo::ParamRef(const int i) const
{
  if (i) {} // warning killer;
  // should  throw
  cout << "Gtransfo::ParamRef should never be called " << endl;
  abort();
  return 0; 
}

double &Gtransfo::ParamRef(const int i)
{
  if (i) {} // warning killer;
  // should  throw
  cout << "Gtransfo::ParamRef should never be called " << endl;
  abort();
  return *(double *)(NULL); 
}

void Gtransfo::ParamDerivatives(const Point &Where, 
				double *Derivatives) const
{
  if (Where.x || Derivatives) {} // compilation warning killer
  // should  throw
  cout << "Gtransfo::ParamDerivatives() should never be called " << endl;
  abort();
}
  

ostream & operator << (ostream &stream, const Gtransfo & T)
           {T.dump(stream); return stream;}


void Gtransfo::Write(const std::string &FileName) const
{
  ofstream s(FileName.c_str());
  Write(s);
  bool ok = !s.fail();
  s.close();
  if (!ok)
    throw(PolokaException("  Gtransfo::Write, something went wrong for file " + FileName )); 
}

void Gtransfo::Write(ostream &stream) const
{
  throw(PolokaException(" Gtransfo::Write should never be called : miss an implementation in some derived class "));
}


/******************* GTransfoInverse ****************/
/* inverse transformation, solved by iterations. Before using 
   it (probably via Gtransfo::InverseTransfo), consider 
   seriously StarMatchList::InverseTransfo */
class GtransfoInverse : public Gtransfo {

private:
  Gtransfo *direct;
  Gtransfo *roughInverse;
  double precision2;
  

public:
  GtransfoInverse(const Gtransfo* Direct, 
		  const double Precision,
		  const Frame& Region);

  //! implements an iterative (Gauss-Newton) solver. It resorts to the Derivative function: 4 calls to the direct transfo per iteration.
  void apply(const double Xin, const double Yin, 
	     double &Xout, double  &Yout) const;

  void dump(ostream &stream) const;

  double fit(const StarMatchList &List);

  virtual Gtransfo *Clone() const;

  GtransfoInverse(const GtransfoInverse&);

  //! Overload the "generic routine" 
  Gtransfo* RoughInverse(const Frame &Region) const
  {
    if (&Region) {} //
    return direct->Clone();
  }

  //! Inverse transfo: returns the direct one!
  Gtransfo* InverseTransfo(const double Precision,
			   const Frame& Region) const
  {
    if (&Region || Precision ) {} //
    return direct->Clone();
  }



  ~GtransfoInverse();

private:
  void operator = (const GtransfoInverse &);

};

Gtransfo* Gtransfo::InverseTransfo(const double Precision,
				   const Frame& Region) const
{
  return new GtransfoInverse(this,Precision,Region);
}


GtransfoInverse::GtransfoInverse(const Gtransfo* Direct, 
				 const double Precision,
				 const Frame& Region)

{
  direct = Direct->Clone();
  roughInverse = Direct->RoughInverse(Region);
  precision2 = Precision*Precision;
}

GtransfoInverse::GtransfoInverse(const GtransfoInverse& Model) : Gtransfo()
{
  direct = Model.direct->Clone();
  roughInverse = Model.roughInverse->Clone();
  precision2 = Model.precision2;
}

GtransfoInverse::~GtransfoInverse()
{
  delete direct;
  delete roughInverse;
}

void GtransfoInverse::operator = (const GtransfoInverse & Model)
{
  if (direct) delete direct; direct = Model.direct->Clone();
  if (roughInverse) delete roughInverse; 
  roughInverse = Model.roughInverse->Clone();
  precision2 = Model.precision2;
}

void GtransfoInverse::apply(const double Xin, const double Yin, 
			    double &Xout, double  &Yout) const
{
  Point in(Xin,Yin);
  Point outGuess = roughInverse->apply(in);
  GtransfoLin directDer, reverseDer;
  int loop = 0;
  int maxloop = 20;
  double move2;
  do
    {
      loop++;
      Point inGuess = direct->apply(outGuess);
      direct->Derivative(outGuess,directDer);
      reverseDer = directDer.invert();
      double xShift, yShift;
      reverseDer.apply(Xin - inGuess.x, Yin - inGuess.y, xShift, yShift);
      outGuess.x += xShift;
      outGuess.y += yShift;
      move2 = xShift*xShift+yShift*yShift;
    } while (( move2 > precision2) && (loop < maxloop));
  if (loop == maxloop)
    {
      cerr << " troubles with Gtransfo inversion at " << in << endl;
    }
  Xout = outGuess.x;
  Yout = outGuess.y;
}

void GtransfoInverse::dump(ostream &stream) const
{
  stream << " GtransfoInverse of  :" << endl
	 << *direct << endl;
}

double GtransfoInverse::fit(const StarMatchList &List)
{
  if (&List) {} // warning killer  
  std::cerr << " Trying to fit a GtransfoInverse... \
try to use StarMatchList::inverseTransfo instead " 
	    << std::endl;
  return -1;
}

Gtransfo *GtransfoInverse::Clone() const
{
  return new GtransfoInverse(*this);
}


/************* GtransfoComposition **************/


// This class was done to allow composition of Gtransfo's, without specifications of their types.
// does not need to be public. Invoked  by GtransfoCompose(Left,Right)
// TODO : use CountedRefs instead of pointers



//! Private class to handle Gtransfo compositions (i.e. piping). Use the routine GtransfoCompose if you need this functionnality.
class GtransfoComposition : public Gtransfo {
  private :
    Gtransfo* first, *second;
  public :
    //! will pipe transfos
    GtransfoComposition(const Gtransfo *Second, const Gtransfo *First);

    //! return Second(First(Xin,Yin))
    void apply(const double Xin, const double Yin, double &Xout, double &Yout) const;
    void dump(ostream &stream = cout) const; 

    //!
    double fit(const StarMatchList &List);

    Gtransfo *Clone() const;
    ~GtransfoComposition();

    //#ifndef SWIG
    //  ClassDef(GtransfoComposition,1);
    //#endif /*SWIG */
};


GtransfoComposition::GtransfoComposition(const Gtransfo *Second, const Gtransfo *First)
{
first =  First->Clone(); 
second = Second->Clone();
}

void GtransfoComposition::apply(const double Xin, const double Yin, double &Xout, double &Yout) const
{
double xout,yout;
first->apply(Xin,Yin, xout,yout);
second->apply(xout,yout,Xout,Yout);
}

void GtransfoComposition::dump(ostream &stream) const
{
first->dump(stream); second->dump(stream);
}

double GtransfoComposition::fit(const StarMatchList &List)
{
  /* fits only one of them. could check that first can actually be fitted... */
  return first->fit(List);
}

Gtransfo *GtransfoComposition::Clone() const
{
return new GtransfoComposition(second,first);
}

GtransfoComposition::~GtransfoComposition()
{
delete first; delete second;
}

/*!  This routine implements "run-time" compositions. When
 there is a possible "reduction" (e.g. compositions of polynomials),
 GtransfoCompose detects it and returns a genuine Gtransfo.
 */
Gtransfo *GtransfoCompose(const Gtransfo *Left, const Gtransfo *Right)
{
  /* is Right Identity ? if Left is Identity , GtransfoIdentity::ReduceCompo does the right job */
  if (IsIdentity(Right))
    {
      return Left->Clone();
    }
  /* Try to use the ReduceCompo method from Left. If absent,
     Gtransfo::ReduceCompo return NULL. ReduceCompo is non trivial for
     polynomials */
  Gtransfo *composition = Left->ReduceCompo(Right);
  /* composition == NULL means no reduction : just build a Composition
     that pipelines "Left" and "Right" */
  if (composition == NULL) return new GtransfoComposition(Left,Right);
  else return composition;
}


// just a speed up, to avoid useless numerical derivation.
void GtransfoIdentity::Derivative(const Point &Where, 
				  GtransfoLin &Derivative, 
				  const double Step) const
{
  if (Step  || &Where) {} // warning killer
  Derivative = GtransfoLin();
}


GtransfoLin GtransfoIdentity::LinearApproximation(const Point &Where, const double Step) const
{
  if (Step  || &Where) {} // warning killer
  GtransfoLin result;
  return result; // rely on default Gtransfolin constructor;
}


void GtransfoIdentity::Write(ostream &s) const
{
  s << "GtransfoIdentity 1" << endl;
}




void GtransfoIdentity::Read(istream &s)
{
  int format;
  s >> format;
  if (format != 1)
    throw(PolokaException(" GtransfoIdentity::Read : format is not 1 " ));
}


/***************  GtransfoPoly **************************************/



//! Default transfo : identity for all degrees (>=1 )

GtransfoPoly::GtransfoPoly(const unsigned Deg)
{
  nterms = 0;
  SetDegree(Deg);
}
  

void GtransfoPoly::SetDegree(const unsigned Deg)
{
  deg = Deg;
  nterms = (deg+1)*(deg+2)/2;

  // allocate vectors, since we know their length
  monomials.reserve(nterms);
  coeffs.reserve(2*nterms);
  // fill them
  monomials.insert(monomials.begin(), nterms,0);
  coeffs.insert(coeffs.begin(), 2*nterms, 0);
  // default is supposed to be the identity
  if (deg>=1)
    {
      Coeff(1,0,0) = 1;
      Coeff(0,1,1) = 1;
    }
}

/* this is reasonably fast, when optimized */
void GtransfoPoly::apply(const double Xin, const double Yin, 
			 double &Xout, double &Yout) const
{
  /* The ordering of monomials is implemented here.
     You may not change it without updating the "mapping" routines
    Coeff(unsigned, unsigned, unsigned). 
    I (P.A.) did not find a clever way to loop over monomials. 
    Improvements welcome.
    This routine is used also by the fit to fill monomials.
    We could certainly be more elegant.

    This routine computes the monomials only once for both
    polynomials.  This is why GtransfoPoly does not use an auxilary
    class (such as PolyXY) to handle each polynomial.

    The code works even if &Xin == &Xout (or &Yin == &Yout)
  */
  vector<double> &m = (vector<double>&) monomials; //constness violation

  double xx = 1;
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      double yy = 1;
      unsigned k=ix*(ix+1)/2;
      for (unsigned iy = 0; iy<=deg-ix; ++iy)
	{
	m[k] = xx*yy;
	yy *= Yin;
	k+= ix+iy+2;
      }
    xx *= Xin;
    }
  Xout = 0;
  Yout = 0;
  const double *c = &coeffs[0];
  const double *pm = &monomials[0];
  for (int k=nterms; k--; ) Xout +=  (*(pm++))*(*(c++));
  pm = &monomials[0];
  for (int k=nterms; k--; ) Yout +=  (*(pm++))*(*(c++));
}


void GtransfoPoly::Derivative(const Point &Where,
			      GtransfoLin &Der, const double Step) const
{ /* routine checked against numerical derivatives from Gtransfo::Derivative */
  if (deg == 1)
    {
      Der = GtransfoLin(*this);
      Der.dx() = Der.dy() = 0;
      return;
    }

  double *dermx = new double [2*nterms];
  double *dermy = dermx+nterms;
  double xin = Where.x;
  double yin = Where.y;

  double xx = 1;
  double xxm1 = 1; // xx^(ix-1)
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      unsigned k=(ix)*(ix+1)/2;
      // iy = 0
      dermx[k] = ix*xxm1; 
      dermy[k] = 0;
      k+= ix+2;
      double yym1 = 1; // yy^(iy-1)
      for (unsigned iy = 1; iy<=deg-ix; ++iy)
	{
	  dermx[k] = ix*xxm1*yym1*yin;
	  dermy[k] = iy*xx*yym1;
	  yym1 *= yin;
	  k+= ix+iy+2;
	}
    xx *= xin;
    if (ix>=1) xxm1 *= xin;
    }

  Der.dx() = 0;
  Der.dy() = 0;

  const double *mx = &dermx[0];
  const double *my = &dermy[0];
  const double *c = &coeffs[0];
  // dx' 
  double a11=0, a12 = 0; 
  for (int k=nterms; k--; )
    {
      a11 += (*(mx++))*(*c);
      a12 += (*(my++))*(*(c++));
    }
  Der.a11() = a11;
  Der.a12() = a12;
  // dy'
  double a21 = 0, a22 = 0;
  mx = &dermx[0];
  my = &dermy[0];
  for (int k=nterms; k--; )
    {
      a21 += (*(mx++))*(*c);
      a22 += (*(my++))*(*(c++));
    }
  Der.a21() = a21;
  Der.a22() = a22;

  delete [] dermx;
}

void GtransfoPoly::TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const
{
  /* 
     The results from this routine were compared to what comes out from apply and 
     TransformErrors. The Derivative routine was checked against
     numerical derivatives from Gtransfo::Derivative. (P.A dec 2009).

     This routine could be made much simpler by calling apply and Derivative
     (i.e. you just suppress it, and the fallback is the generic version in Gtransfo).
     BTW, I checked that both routines provide the same result.
     This version is however faster (monomials get recycled), but we probably don't really care.
  */
  vector<double> &m = (vector<double>&) monomials; //constness violation
  FatPoint  res; // to store the result, because nothing forbids &In == &Out.
  // I would like to avoid that : 
  double *dermx = new double [2*nterms]; // monomials for derivative w.r.t. x
  double *dermy = dermx+nterms;  // same for y
  double xin = In.x;
  double yin = In.y;

  double xx = 1;
  double xxm1 = 1; // xx^(ix-1)
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      unsigned k=(ix)*(ix+1)/2;
      // iy = 0
      dermx[k] = ix*xxm1; 
      dermy[k] = 0;
      m[k] = xx;
      k+= ix+2;
      double yy = yin;
      double yym1 = 1; // yy^(iy-1)
      for (unsigned iy = 1; iy<=deg-ix; ++iy)
	{
	  m[k] = xx*yy;
	  dermx[k] = ix*xxm1*yy;
	  dermy[k] = iy*xx*yym1;
	  yym1 *= yin;
	  yy *= yin;
	  k+= ix+iy+2;
	}
    xx *= xin;
    if (ix>=1) xxm1 *= xin;
    }

  // output position 
  double xout = 0, yout=0;
  const double *c = &coeffs[0];
  const double *pm = &monomials[0];
  for (int k=nterms; k--; ) xout +=  (*(pm++))*(*(c++));
  pm = &monomials[0];
  for (int k=nterms; k--; ) yout +=  (*(pm++))*(*(c++));
  res.x = xout; res.y = yout;

  // derivatives 
  c = &coeffs[0];
  const double *mx = &dermx[0];
  const double *my = &dermy[0];
  double a11=0, a12 = 0; 
  for (int k=nterms; k--; )
    {
      a11 += (*(mx++))*(*c);
      a12 += (*(my++))*(*(c++));
    }

  double a21 = 0, a22 = 0;
  mx = &dermx[0];
  my = &dermy[0];
  for (int k=nterms; k--; )
    {
      a21 += (*(mx++))*(*c);
      a22 += (*(my++))*(*(c++));
    }

  // output co-variance
  res.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  res.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  res.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;
  Out = res;
  delete [] dermx;// cleanup  
}

/* The coefficient ordering is defined both here *AND* in the 
   GtransfoPoly::apply, GtransfoPoly::Derivative, ... routines
   Change all or none ! */

double GtransfoPoly::Coeff(const unsigned Degx, const unsigned Degy,
			   const unsigned WhichCoord) const
{
  assert((Degx+Degy<=deg) && WhichCoord<2);
  /* this assertion above is enough to ensure that the index used just
     below is within bounds since the reserved length is
     2*nterms=(deg+1)*(deg+2) */
  return coeffs[(Degx+Degy)*(Degx+Degy+1)/2+Degy+WhichCoord*nterms];
}


double& GtransfoPoly::Coeff(const unsigned Degx, const unsigned Degy,
			   const unsigned WhichCoord)
{
  assert((Degx+Degy<=deg) && WhichCoord<2);                       
  return coeffs[(Degx+Degy)*(Degx+Degy+1)/2+Degy+WhichCoord*nterms];
}

/* parameter serialization for "virtual" fits */
double GtransfoPoly::ParamRef(const int i) const
{
  assert(unsigned(i)<2*nterms);
  return coeffs[i];
}


double& GtransfoPoly::ParamRef(const int i)
{
  assert(unsigned(i)<2*nterms);
  return coeffs[i];
}

void GtransfoPoly::ParamDerivatives(const Point &Where, 
				    double *Derivatives) const
{/* first half : dxout/dpar, second half : dyout/dpar */
  double xout,yout;
  apply(Where.x, Where.y, xout,yout); // to get the monomials
  for (unsigned k=0; k<nterms; ++k)
    {
      Derivatives[k] = Derivatives[3*nterms+k] = monomials[k];
      Derivatives[nterms+k] = Derivatives[2*nterms+k] = 0;
    }
}

/*
  mapping coefficients with names (inherited from GtransfoLin,GtransfoQuad,GtransfoCub} */
typedef struct
{
  const char *name;
  const unsigned char px,py,whichCoord;
} CoeffTagStruct;
  
static const CoeffTagStruct CoeffTags[] =
  {
    {"dx",   0,0,0},
    {"dy",   0,0,1},
    {"a11",  1,0,0},
    {"a12",  0,1,0},
    {"a22",  0,1,1},
    {"a21",  1,0,1},
    {"a1x2", 2,0,0},
    {"a1xy", 1,1,0},
    {"a1y2", 0,2,0},
    {"a2y2", 0,2,1},
    {"a2xy", 1,1,1},
    {"a2x2", 2,0,1},
    {"a1x3", 3,0,0},
    {"a1x2y",2,1,0},
    {"a1xy2",1,2,0},
    {"a1y3", 0,3,0},
    {"a2y3", 0,3,1},
    {"a2xy2",1,2,1},
    {"a2x2y",2,1,1},
    {"a2x3", 3,0,1}
  };

static unsigned NCoeffTags = sizeof(CoeffTags)/sizeof(CoeffTags[0]);


static unsigned tag_pos(const char *Name)
{
  for (unsigned k=0; k<NCoeffTags; ++ k)
    {
      const CoeffTagStruct &t=CoeffTags[k];
      if (strcmp(t.name, Name) == 0)
	return k;
    }
  stringstream message;
  message << "GtransfoPoly::Coeff(const char *Name) : unknown name : \""
	  << string(Name) << '\"';
  abort();
  //throw(GtransfoException(message));
}

bool GtransfoPoly::HasCoeff(const char* Name) const
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];  
  return (t.px+t.py<=deg);
}

double GtransfoPoly::Coeff(const char *Name) const 
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];
  return Coeff(t.px, t.py, t.whichCoord);
}

double& GtransfoPoly::Coeff(const char *Name)
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];
  return Coeff(t.px, t.py, t.whichCoord);
}

/* utility for the dump(ostream&) routine */
static string monomial_string(const unsigned powx, const unsigned powy)
{
  stringstream ss;
  if (powx+powy) ss<<"*";
  if (powx>0) ss<< "x";
  if (powx>1) ss<<"^"<<powx;
  if (powy>0) ss<< "y";
  if (powy>1) ss<<"^"<<powy;
  return ss.str();
}
  
void  GtransfoPoly::dump(ostream &S) const
{
  for (unsigned ic=0; ic<2; ++ic)
    {
      if (ic==0)   S << " newx = ";
      else S << " newy = ";
      for (unsigned p = 0; p<=deg; ++p)
	for (unsigned py=0; py<=p; ++py)
	  { 
	    if (p+py != 0) S<< " + "; 
	    S << Coeff(p-py,py,ic) << monomial_string(p-py,py);
	  }
      S << endl;
    }
  S << " Linear Determinant = " << Determinant() << endl ;
}

double GtransfoPoly::Determinant() const
{
  return Coeff(1,0,0)*Coeff(0,1,1) - Coeff(0,1,0)*Coeff(1,0,1);
}



/*utility for the GtransfoPoly::fit() routine */
static GtransfoPoly my_shift_to_center(const StarMatchList &List)
{
  double xav=0;
  double yav=0;
  double count=0;
  for (StarMatchCIterator it = List.begin(); it != List.end(); ++it)
    {
      const StarMatch &a_match = *it;
      const Point &point1 = a_match.point1;
      xav += point1.x;
      yav += point1.y;
      count++;
    }
  if (count==0) count = 1;
  GtransfoPoly result(1); // default = identity
  result.Coeff(0,0,0) = -xav/count;
  result.Coeff(0,0,1) = -yav/count;
  return result;
}

static double sq(const double &x) { return x*x;}
  
double GtransfoPoly::do_the_fit(const StarMatchList &List, 
				const Gtransfo &ShiftToCenter,
				const bool UseErrors)
{
  Mat A(2*nterms,2*nterms);
  Vect B(2*nterms);
  double sumr2 = 0;
  for (StarMatchCIterator it = List.begin(); it != List.end(); ++it)
    {
      const StarMatch &a_match = *it;
      Point tmp = ShiftToCenter.apply(a_match.point1);
      FatPoint point1(tmp, a_match.point1.vx, a_match.point1.vy, 
		      a_match.point1.vxy);
      const FatPoint &point2 = a_match.point2;
      double wxx,wyy,wxy;
      FatPoint tr1;
      if (UseErrors)
	{
	  TransformPosAndErrors(point1, tr1);// it also fills the monomials
	  double vxx = (tr1.vx+point2.vx);
	  double vyy = (tr1.vy+point2.vy);
	  double vxy = (tr1.vxy+point2.vxy);
	  double det = vxx*vyy-vxy*vxy;
	  wxx = vyy/det;
	  wyy = vxx/det;
	  wxy = -vxy/det;
	}
      else 
	{
	  wxx = wyy = 1; wxy = 0;
	  apply(point1.x, point1.y ,tr1.x, tr1.y);
	}
      double resx = point2.x - tr1.x;
      double resy = point2.y - tr1.y;
      sumr2 += wxx*sq(resx) + wyy*sq(resy) 
	+2*wxy*resx*resy;

      double bxcoeff = wxx*resx + wxy*resy ;
      double bycoeff = wyy*resy + wxy*resx;
      for (unsigned j=0; j<nterms; ++j)
	{
	  for (unsigned i=0; i<=j; ++i)
	    {
	      A(i,j) += wxx*monomials[i]*monomials[j];
	      A(i+nterms,j+nterms) += wyy*monomials[i]*monomials[j];
	      A(j,i+nterms) = A(i,j+nterms) += wxy*monomials[i]*monomials[j];
	    }
	  B(j)        += bxcoeff*monomials[j];
	  B(j+nterms) += bycoeff*monomials[j];
	}
    } // end loop on points

  Vect sol(B);
  if (cholesky_solve(A,sol,"U") != 0) return false;
  for (unsigned k=0; k< 2*nterms; ++k) coeffs[k] += sol(k);
  if (List.size() == nterms) return 0;
  return (sumr2-B*sol);
}


double  GtransfoPoly::fit(const StarMatchList &List)
{
  if (List.size()< nterms)
    {
      cerr << " GtransfoPoly::fit : trying to fit a polynomial transfo of degree " << deg << " with only " << List.size() << " matches " << endl;
      return -1;
    }

  GtransfoPoly shift_to_center = my_shift_to_center(List);

  do_the_fit(List, shift_to_center, false); // get a rough solution
  do_the_fit(List, shift_to_center, true); // weight with it
  double chi2 = do_the_fit(List, shift_to_center, true); // once more
  
  (*this) = (*this)*shift_to_center;
  if (List.size() == nterms) return 0;
  return chi2;
}


Gtransfo * GtransfoPoly::ReduceCompo(const Gtransfo *Right) const
{
  const GtransfoPoly *p = dynamic_cast<const GtransfoPoly *>(Right);
  if (p) 
    {
      if (Degree() == 1 && p->Degree() == 1)
	return new GtransfoLin((*this)*(*p)); // does the composition
      else
	return new GtransfoPoly((*this)*(*p)); // does the composition
    }
  else return NULL;
}

/*  PolyXY the class used to perform polynomial algebra (and in
    particular composition) at the coefficient level. This class
    handles a single polynomial, while a GtransfoPoly is a couple of
    polynomials. This class does not have any routine to evaluate
    polynomials. Efficiency is not a concern since these routines are
    seldom used.  There is no need to expose this tool class to
    Gtransfo users.
*/


class PolyXY
{
  unsigned deg;
  unsigned nterms;
  vector<long double> coeffs;

public :

  PolyXY(const int Deg) : deg(Deg), nterms((deg+1)*(deg+2)/2)
  {
    coeffs.reserve(nterms);
    coeffs.insert(coeffs.begin(), nterms, 0L); // fill & initialize to 0.
  }

  unsigned Deg() const { return deg;}

  PolyXY(const GtransfoPoly &P, const unsigned WhichCoord) 
    : deg(P.Degree()) , nterms((deg+1)*(deg+2)/2) , coeffs(nterms, 0L)
  {
    for (unsigned px=0; px<=deg; ++px)
      for (unsigned py=0; py<=deg-px; ++py)
	Coeff(px,py) = P.Coeff(px,py,WhichCoord);
  }
	

  long double Coeff(const unsigned powx, const unsigned powy) const
  {
    assert(powx+powy<=deg);
    return coeffs.at((powx+powy)*(powx+powy+1)/2+powy);
  }

  long double &Coeff(const unsigned powx, const unsigned powy)
  {
    assert(powx+powy<=deg);
    return coeffs.at((powx+powy)*(powx+powy+1)/2+powy);
  }

};


/* =====================  PolyXY Algebra routines ================== */

static void operator += (PolyXY &Left, const PolyXY &Right)
{
  unsigned rdeg = Right.Deg();
  assert(Left.Deg()>= rdeg);
  for (unsigned i=0; i<= rdeg; ++i)
    for (unsigned j=0; j<= rdeg-i; ++j)
      Left.Coeff(i,j) += Right.Coeff(i,j);
}

/* multiplication by a scalar */
static PolyXY operator * (const long double &a, const PolyXY &P)
{
  PolyXY result(P);
  // no direct access to coefficients: do it the soft way
  unsigned deg = P.Deg();
  for (unsigned i=0; i<=deg; ++i)
    for (unsigned j=0; j<=deg-i; ++j)
      result.Coeff(i,j) *= a;
  return result;
}


/*! result(x,y) = P1(x,y)*P2(x,y) */
static PolyXY Product(const PolyXY &P1, const PolyXY &P2)
{
  unsigned deg1 = P1.Deg();
  unsigned deg2 = P2.Deg();
  PolyXY result(deg1+deg2);
  for (unsigned i1=0; i1<=deg1; ++i1)
    for (unsigned j1=0; j1<=deg1-i1; ++j1)
      for (unsigned i2=0; i2<=deg2; ++i2)
	for (unsigned j2=0; j2<=deg2-i2; ++j2)
	  result.Coeff(i1+i2,j1+j2) += P1.Coeff(i1,j1)*P2.Coeff(i2,j2);
  return result;
}
	
/* Powers[k](x,y) = P(x,y)**k, 0 <= k <= MaxP */ 
static void ComputePowers(const PolyXY &P, const unsigned MaxP, vector<PolyXY> &Powers)
{
  Powers.reserve(MaxP+1);
  Powers.push_back(PolyXY(0)); Powers[0].Coeff(0,0) = 1L;
  for (unsigned k=1; k<=MaxP; ++k) Powers.push_back(Product(Powers[k-1],P));
}

/*! result(x,y) = P(Px(x,y),Py(x,y)) */
static PolyXY Composition(const PolyXY &P, const PolyXY &Px, const PolyXY &Py)
{
  unsigned pdeg = P.Deg();
  PolyXY result(pdeg*max(Px.Deg(), Py.Deg()));
  vector<PolyXY> PxPowers;
  vector<PolyXY> PyPowers;
  ComputePowers(Px, pdeg, PxPowers);
  ComputePowers(Py, pdeg, PyPowers);
  for (unsigned px=0 ; px <= pdeg; ++px)
    for (unsigned py=0; py <= pdeg-px; ++py)
      result += P.Coeff(px,py)*Product(PxPowers.at(px), PyPowers.at(py));
  return result;
}

/* ===================== end of  PolyXY Algebra routines ============= */


/* reducing polynomial composition is the reason for PolyXY stuff : */

GtransfoPoly  GtransfoPoly::operator*(const GtransfoPoly &Right) const
{
  // split each transfo into 2d polynomials
  PolyXY plx(*this, 0);
  PolyXY ply(*this, 1);
  PolyXY prx(Right,0);
  PolyXY pry(Right,1);

  // compute the compositions
  PolyXY rx(Composition(plx,prx,pry));
  PolyXY ry(Composition(ply,prx,pry));

  //copy the results the hard way.
  GtransfoPoly result(deg*Right.deg);
  for (unsigned px=0; px<=result.deg; ++px)
    for (unsigned py=0; py<=result.deg-px; ++py)
      {
	result.Coeff(px,py,0) = rx.Coeff(px,py);
	result.Coeff(px,py,1) = ry.Coeff(px,py);
      }
  return result;
}


void GtransfoPoly::Write(ostream &s) const
{
  s << " GtransfoPoly 1"<< endl;
  s << "degree " << deg << endl;
  int oldprec=s.precision();
  s << setprecision(12);
  for (unsigned k=0;k<2*nterms; ++k)
    s << coeffs[k] << ' ';
  s << endl;
  s << setprecision(oldprec);
}

void GtransfoPoly::Read(istream &s)
{
  int format;
  s >> format;
  if (format != 1)
    throw(PolokaException(" GtransfoPoly::Read : format is not 1 " ));
  string degree;
  s >> degree >> deg;
  if (degree != "degree")
    throw(PolokaException(" GtransfoPoly::Read : expecting \"degree\" and found "+degree ));
  SetDegree(deg);
  for (unsigned k=0;k<2*nterms; ++k)
    s >> coeffs[k];
}  


GtransfoPoly *InversePolyTransfo(const Gtransfo &Direct, const Frame &F, const double Prec)
{
  StarMatchList sm;
  unsigned nx = 10;
  double stepx = F.Nx()/(nx+1);
  unsigned ny = 10;
  double stepy= F.Ny()/(ny+1);
  for (unsigned i=0 ;i<nx; ++i)
    for (unsigned j=0; j<ny; ++j)
      {
	Point in((i+0.5)*stepx, (j+0.5)*stepy);
	Point out(Direct.apply(in));
	sm.push_back(StarMatch(out, in , NULL,NULL));
      }
  unsigned npairs = sm.size();
  int maxdeg = 4;
  int degree;
  GtransfoPoly *poly = NULL;
  for (degree=1; degree<=maxdeg; ++degree)
    {
      delete poly;
      poly = new GtransfoPoly(degree);
      double chi2 = poly->fit(sm);
      if (chi2/npairs< Prec*Prec) break;
    }
  if (degree>maxdeg)
    cout << " InversePolyTransfo : Reached  max degree without reaching  requested precision = " << Prec << endl;
  return poly;
}  



/**************** GtransfoLin ***************************************/
/* GtransfoLin is a specialized constructor of GtransfoPoly 
   May be it could just disappear ??
*/

GtransfoLin::GtransfoLin(const double Dx, const double Dy , 
			 const double A11, const double A12, 
			 const double A21, const double A22) : GtransfoPoly(1)
{
  dx()  = Dx;
  a11() = A11;
  a12() = A12;
  dy()  = Dy;
  a21() = A21;
  a22() = A22;
}


GtransfoLin::GtransfoLin(const GtransfoPoly &P) : GtransfoPoly(1)
{
  if (P.Degree() !=  1) 
    {
      cout << " Trying to build a GtransfoLin from a higher order transfo. aborting " << endl;
      abort(); // should throw
    }
  (GtransfoPoly &) (*this) = P;
}
  

GtransfoLin GtransfoLin::operator*(const  GtransfoLin &Right) const
{
  // There is a general routine in GtransfoPoly that would do the job: 
  //  return GtransfoLin(GtransfoPoly::operator*(Right));
  // however, we are using this composition of linear stuff heavily in 
  // Gtransfo::LinearApproximation, itself used in InverseTransfo::apply.
  // So, for sake of efficiency, and since it is easy, we take a shortcut:
  GtransfoLin result;
  apply(Right.Dx(), Right.Dy(), result.dx(), result.dy());
  result.a11() = this->A11()*Right.A11() + this->A12()*Right.A21();
  result.a12() = this->A11()*Right.A12() + this->A12()*Right.A22();
  result.a21() = this->A21()*Right.A11() + this->A22()*Right.A21();
  result.a22() = this->A21()*Right.A12() + this->A22()*Right.A22();
  return result;
}

void GtransfoLin::Derivative(const Point &Where, GtransfoLin &Derivative, 
			     const double Step) const
{
  if (Step || &Where) {} // "unused variable" warning killer
  Derivative = *this;
  Derivative.Coeff(0,0,0) = 0;
  Derivative.Coeff(0,0,1) = 0;
}


GtransfoLin GtransfoLin::LinearApproximation(const Point &Where, 
					     const double Step) const
{
  if (Step || &Where) {} // "unused variable" warning killer
  return *this;
}  

GtransfoLin GtransfoLin::invert() const
{
  //
  //   (T1,M1) * (T2,M2) = 1 i.e (0,1) implies
  //   T1 = -M1 * T2
  //   M1 = M2^(-1)
  //

  double a11 = A11();
  double a12 = A12();
  double a21 = A21();
  double a22 = A22();
  double d = (a11*a22 - a12*a21);
  if (d == 0)
    {
      cerr << " trying to invert a singular transformation: a (nice) crash follows" << endl;
      dump(cerr);
    }

  GtransfoLin result(0,0,a22/d,-a12/d,-a21/d,a11/d);
  double rdx,rdy;
  result.apply(Dx(),Dy(),rdx,rdy);
  result.dx() = -rdx;
  result.dy() = -rdy;
  return result;
} 

Gtransfo* GtransfoLin::InverseTransfo(const double Precision,
				      const Frame& Region) const
{
  if (&Region || Precision) {} // warning killer
  return new GtransfoLin(this->invert());
}

double  GtransfoLinRot::fit(const StarMatchList &List)
{

  cout << " GTransfoLinRot::fit : not implemented !!! aborting " << endl;
  abort();
  return -1;
}


double  GtransfoLinShift::fit(const StarMatchList &List)
{
  int npairs = List.size();
  if (npairs < 3)
    {
      // cerr << " GtransfoLinShift::fit  : trying to fit a linear transfo with only " << npairs << " matches " << endl;
      return -1;
    }

  double sumr2 = 0; /* used to compute chi2 without relooping */
  /* loop on pairs  and fill */
  Vect B(2);
  Mat A(2,2);

  for (StarMatchCIterator it = List.begin(); it != List.end(); it++)
    {
      const FatPoint &point1 = it->point1;
      const FatPoint &point2 = it->point2;      
      double deltax = point2.x - point1.x;
      double deltay = point2.y - point1.y;
      double vxx = point1.vx+point2.vx;
      double vyy = point1.vy+point2.vy;
      double vxy = point1.vxy+point2.vxy;
      double det = vxx*vyy-vxy*vxy;      
      double wxx = vyy/det;
      double wyy = vxx/det;
      double wxy = -vxy/det;
      B(0) +=  deltax*wxx + wxy*deltay;
      B(1) +=  deltay*wyy + wxy*deltax;
      A(0,0) += wxx;
      A(1,1) += wyy;
      A(0,1) += wxy;
      sumr2 += deltax*deltax*wxx + deltay*deltay*wyy 
	+2.*wxy*deltax*deltay;
    }
  A(1,0) = A(0,1); 
  Vect sol(B);
  if (cholesky_solve(A,sol,"U") != 0) return -1;
  (*this) = GtransfoLinShift(sol(0), sol(1));
  return (sumr2  - sol*B); // chi2 compact form
}


GtransfoLinRot::GtransfoLinRot(const double AngleRad, const Point *Center, 
			       const double ScaleFactor)
{
  double c = ScaleFactor*cos(AngleRad);
  double s = ScaleFactor*sin(AngleRad);
  a11() = a22() = c;
  a21() = s;
  a12() = -s;
  
  // we want that the center does not move : T+M*C = C ==> T = C - M*C
  Point a_point(0.,0.);
  if (Center) a_point = *Center;

  dx() = dy() = 0;
  apply(a_point.x, a_point.y, dx(), dy()); // compute M*C
  dx() = a_point.x - Dx(); dy() = a_point.y - dy();
}





static double deg2rad(const double &deg)
{
  return deg*M_PI/180.;
}


static double rad2deg(const double &rad)
{
  return rad*180./M_PI;
}

/*************  WCS transfo ******************/
/************** LinPix2Tan *******************/

/* Implementation note : it seemed wise to incorporate
   the radians to degrees convertion into the linPix2Tan
   part (right in the constructor), and to do the 
   opposite operation in the LinPart routine.
   When I was coding the fit, I realized that it was a 
   bad idea. then I realized that the fitting routine 
   itself was probably useless. Finaly, for a persistor,
   it seems bizarre that the stored data is different
   from what convention (angles in degrees for WCS's)
   would expect.
   So, no more "automatic" degrees to radians and 
   radians to degrees conversion. They are explicitely 
   done in apply (for TanPix2RaDec and TanRaDec2Pix).
   This is a minor concern though....
*/

TanPix2RaDec::TanPix2RaDec(const GtransfoLin &Pix2Tan, 
			   const Point &TangentPoint, 
			   const GtransfoPoly* Corrections)
{
  /* the angles returned by linPix2Tan should be in 
     degrees. */
  linPix2Tan = Pix2Tan;
  ra0  = deg2rad(TangentPoint.x);
  dec0 = deg2rad(TangentPoint.y);
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = NULL;
  if (Corrections) corr = new GtransfoPoly(*Corrections);
}

// ": Gtransfo" suppresses a warning
TanPix2RaDec::TanPix2RaDec(const TanPix2RaDec &Original) : Gtransfo()
{
  corr = NULL;
  *this = Original;
}

void TanPix2RaDec::operator = (const TanPix2RaDec &Original)
{
  linPix2Tan = Original.linPix2Tan;
  ra0 = Original.ra0;
  dec0 = Original.dec0;
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = NULL;
  if (Original.corr) corr = new GtransfoPoly(*Original.corr);
}  



// ": Gtransfo" suppresses a warning
TanPix2RaDec::TanPix2RaDec() : Gtransfo(), linPix2Tan()
{
  ra0 = 0; dec0 = 0;
  cos0 = 1; sin0 = 0;
  corr = NULL;
}


void TanPix2RaDec::apply(const double Xin, const double Yin, 
			 double &Xout, double &Yout) const
{
  double l,m; // radians in the tangent plane
  linPix2Tan.apply(Xin,Yin,l,m); // l, m in degrees.
  if (corr) corr->apply(l,m,l,m); // still in degrees.
  l = deg2rad(l); 
  m = deg2rad(m); // now in radians
  double dect = cos0 - m * sin0;
  if (dect == 0)
    {
      cerr << " no sideral coordinates at pole ! " << endl;
      Xout = 0;
      Yout = 0;
      return;
    }
  double rat = ra0 + atan2(l, dect);
  dect = atan(cos(rat-ra0) * (m * cos0 + sin0) / dect);
  if (rat - ra0 >  M_PI) rat -= (2.*M_PI);
  if (rat - ra0 < -M_PI) rat += (2.*M_PI);
  if (rat < 0.0) rat += (2.*M_PI);
  // convert to deg
  Xout = rad2deg(rat);
  Yout = rad2deg(dect);
}


Gtransfo * TanPix2RaDec::ReduceCompo(const Gtransfo *Right) const
{
  const GtransfoLin *lin = dynamic_cast<const GtransfoLin *>(Right);
  if (lin && lin->Degree() == 1) return new TanPix2RaDec((*this)*(*lin));
  return NULL;
}


TanPix2RaDec TanPix2RaDec::operator *(const GtransfoLin &Right) const
{
  TanPix2RaDec result(*this);
  result.linPix2Tan = result.linPix2Tan * Right;
  return result;
}

TanRaDec2Pix TanPix2RaDec::invert() const
{
  if (corr != NULL)
    {
      cerr << " You are inverting a TanPix2RaDec with corrections " << endl;
      cerr << " The inverse you get ignores the corrections !!!!!!" << endl;
    }
  return TanRaDec2Pix(LinPart().invert(),TangentPoint());
}

Gtransfo* TanPix2RaDec::RoughInverse(const Frame &Region) const
{
  if (&Region) {} 
  return new TanRaDec2Pix(LinPart().invert(),TangentPoint());
}

Gtransfo*  TanPix2RaDec::InverseTransfo(const double Precision,
					const Frame& Region) const
{
  if (!corr) return new TanRaDec2Pix(LinPart().invert(),TangentPoint());
  else return new GtransfoInverse(this, Precision, Region);
}

void TanPix2RaDec::SetCorrections(const GtransfoPoly *Corrections)
{
  if (corr) delete corr;
  corr = NULL;
  if (Corrections)
    {
      corr = (GtransfoPoly*) Corrections->Clone();
    }
}

Point TanPix2RaDec::TangentPoint() const
{
  return Point(rad2deg(ra0),rad2deg(dec0));
}

GtransfoLin TanPix2RaDec::LinPart() const
{
  return linPix2Tan;
}


const GtransfoPoly* TanPix2RaDec::Corr() const
{
  return corr;
}

GtransfoPoly TanPix2RaDec::Pix2TangentPlane() const
{
  if (corr) return (*corr)*linPix2Tan;
  else return linPix2Tan;
}


Point TanPix2RaDec::CrPix() const
{
  /* CRPIX's are defined by:
                    ( CD1_1  CD1_2 )   (X - crpix1)
     transformed =  (              ) * (          )
                    ( CD2_1  CD2_2 )   (Y - crpix2)

     so that CrPix is the point which transforms to (0,0) 
  */
  const GtransfoLin inverse = linPix2Tan.invert();
  return Point(inverse.Dx(),inverse.Dy());
}
  

Gtransfo *TanPix2RaDec::Clone() const
{
  return new TanPix2RaDec(LinPart(), TangentPoint(), corr);
}

void TanPix2RaDec::dump(ostream &stream) const
{
  stream << " lin part :" << endl << linPix2Tan;
  Point tp = TangentPoint();
  stream << " tangent point " << tp.x << ' ' << tp.y << endl;
  Point crpix = CrPix();
  stream << " crpix " << crpix.x << ' ' << crpix.y << endl;
  if (corr) stream << " correction: " << endl << *corr;
}


double  TanPix2RaDec::fit(const StarMatchList &List)
{ 
  /* OK we could implement this routine, but it is
     probably useless since to do the match, we have to
     project from sky to tangent plane. When a match is
     found, it is easier to carry out the fit in the 
     tangent plane, rather than going back to the celestial
     sphere (and reproject to fit...). Anyway if this
     message shows up, we'll think about it.
  */
  if (&List) {} // warning killer
  throw(PolokaException("TanPix2RaDec::fit is NOT implemented (although it is doable)) "));
  return -1;
}



TanPix2RaDec::~TanPix2RaDec()
{ 
  if (corr) delete corr;
}

/***************  reverse transfo TanRaDec2Pix ********/


TanRaDec2Pix::TanRaDec2Pix(const GtransfoLin &Tan2Pix, const Point &TangentPoint) : linTan2Pix(Tan2Pix)
{/* the radian to degrees convertion after projection 
    is handled in apply */
  ra0  = deg2rad(TangentPoint.x);
  dec0 = deg2rad(TangentPoint.y);
  cos0 = cos(dec0);
  sin0 = sin(dec0);
}

TanRaDec2Pix::TanRaDec2Pix() : linTan2Pix()
{
  ra0 = dec0 = 0;
  cos0 = 1;
  sin0 = 0;
}
  


Point TanRaDec2Pix::TangentPoint() const
{
  return Point(rad2deg(ra0),rad2deg(dec0));
}

GtransfoLin TanRaDec2Pix::LinPart() const
{
  return linTan2Pix;
}

// Use analytic derivatives, computed at the same time as the transform itself
void TanRaDec2Pix::TransformPosAndErrors(const FatPoint &In, 
					 FatPoint &Out) const
{
  /* this routine is very similar to apply, but also propagates errors.
     The deg2rad and rad2deg are ignored for errors because they act as
     2 global scalings that cancel each other. 
     Derivatives were computed using maple:

     l1 := sin(a - a0) cos(d)
     m1 := sin(d)*sin(d0)+cos(d)*cos(d0)*cos(a-a0);
     l2 := sin(d)*cos(d0)-cos(d)*sin(d0)*cos(a-a0);
     simplify(diff(l1/m1,a);
     simplify(diff(l1/m1,d);
     simplify(diff(l2/m1,a);
     simplify(diff(l2/m1,d);

     Checked against Gtransfo::TransformPosAndErrors (dec 09)
  */
  double ra = deg2rad(In.x);
  double dec = deg2rad(In.y);
  if (ra-ra0 >  M_PI ) ra -= (2.* M_PI);
  if (ra-ra0 < -M_PI ) ra += (2.* M_PI);
  
  double coss = cos(dec);
  double sins = sin(dec);
  double sinda = sin(ra-ra0);
  double cosda = cos(ra -ra0);
  double l = sinda * coss;
  double m = sins * sin0 + coss * cos0 * cosda;
  l = l/m;
  m = (sins * cos0 - coss * sin0 * cosda)/ m;

  // derivatives
  double deno = sq(sin0)-sq(coss)+sq(coss*cos0)*(1+sq(cosda))+2*sins*sin0*coss*cos0*cosda;
  double a11 = coss*(cosda*sins*sin0+coss*cos0)/deno;
  double a12 = sinda*sin0/deno;
  double a21 = coss*sinda*sins/deno;
  double a22 = cosda/deno;
  
  FatPoint tmp;
  tmp.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  tmp.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  tmp.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;

  // l and m are now coordinates in the tangent plane, in radians.
  tmp.x = rad2deg(l);
  tmp.y = rad2deg(m);
  
  linTan2Pix.TransformPosAndErrors(tmp, Out);
}  


void TanRaDec2Pix::apply(const double Xin, const double Yin, double &Xout, double &Yout) const
{
  double ra = deg2rad(Xin);
  double dec = deg2rad(Yin);
  if (ra-ra0 >  M_PI ) ra -= (2.* M_PI);
  if (ra-ra0 < -M_PI ) ra += (2.* M_PI);
  
  double coss = cos(dec);
  double sins = sin(dec);
  double l = sin(ra-ra0) * coss;
  double m = sins * sin0 + coss * cos0 * cos (ra -ra0);
  l = l/m;
  m = (sins * cos0 - coss * sin0 * cos( ra -ra0))/ m;
  // l and m are now coordinates in the tangent plane, in radians.
  l = rad2deg(l);
  m = rad2deg(m);
  linTan2Pix. apply(l,m, Xout, Yout);
}


TanPix2RaDec TanRaDec2Pix::invert() const
{
  return TanPix2RaDec(LinPart().invert(),TangentPoint());
}


void TanRaDec2Pix::dump(ostream &stream) const
{
  Point tp = TangentPoint();
  stream << " tan2pix " << linTan2Pix << " tangent point " << tp.x << ' ' << tp.y << endl;
}

Gtransfo* TanRaDec2Pix::RoughInverse(const Frame &Region) const
{
  if (&Region) {} //
  return new TanPix2RaDec(LinPart().invert(),TangentPoint());
}

Gtransfo* TanRaDec2Pix::InverseTransfo(const double Precision,
					const Frame& Region) const
{
  if (&Region || Precision) {}
  return new TanPix2RaDec(LinPart().invert(),TangentPoint());
}


Gtransfo *TanRaDec2Pix::Clone() const
{
  return new TanRaDec2Pix(*this);
}

double TanRaDec2Pix::fit(const StarMatchList &List)
{
  if (&List) {} // warning killer
  throw(PolokaException("TanRaDec2Pix::fit is NOT implemented (although it is doable)) "));
  return -1;
}

/*************  a "run-time" transfo, that does not require to 
modify this file */

UserTransfo::UserTransfo(GtransfoFun &Fun, const void *UserData):
  userFun(Fun), userData(UserData) 
{}

void UserTransfo::apply(const double Xin, const double Yin, 
			double &Xout, double &Yout) const
{
  userFun(Xin,Yin,Xout,Yout, userData);
}

void UserTransfo::dump(ostream &stream) const
{
  stream << "UserTransfo with user function @ " << userFun  
	 << "and userData@ " << userData << endl;
}

double UserTransfo::fit(const StarMatchList &List)
{
  if (&List) {} // warning killer
  throw(PolokaException("UserTransfo::fit is NOT implemented (and will never be)) "));
  return -1;
}

Gtransfo *UserTransfo::Clone() const
{
  return new UserTransfo(*this);
}


/*************************************************************/


Gtransfo* GtransfoRead(const std::string &FileName)
{
  ifstream s(FileName.c_str());
  if (!s) throw(PolokaException(" GtransfoRead : cannot open " + FileName));
  try
    {
      Gtransfo *res = GtransfoRead(s);
      s.close();
      return res;
    }
  catch (PolokaException e) {
    throw(e.message()+" in file "+FileName);}
}

Gtransfo* GtransfoRead(istream &s)
{
  std::string type;
  s >> type;
  if (s.fail()) throw(PolokaException("GtransfoRead : could not find a Gtransfotype"));
  if (type == "GtransfoIdentity")
    {GtransfoIdentity* res = new GtransfoIdentity(); res->Read(s); return res;}
  else if (type == "GtransfoPoly")
    {GtransfoPoly* res = new GtransfoPoly(); res->Read(s); return res;}
  else
    throw(PolokaException(" No reader for Gtransfo type "+ type));
}


/* clearly this assumes that the Transfo is essentially a translation or a simple rotation (n*pi/2) could probably be improved */
Frame ApplyTransfo(const Frame& inputframe, const Gtransfo &T, const WhichTransformed W) 
{
  // 2 opposite corners
  double xtmin1, xtmax1, ytmin1, ytmax1;
  T.apply(inputframe.xMin,inputframe.yMin,xtmin1,ytmin1);
  T.apply(inputframe.xMax,inputframe.yMax,xtmax1,ytmax1);
  Frame fr1(min(xtmin1,xtmax1), min(ytmin1,ytmax1), 
	    max(xtmin1,xtmax1), max(ytmin1,ytmax1));
  // 2 other corners
  double xtmin2, xtmax2, ytmin2, ytmax2;
  T.apply(inputframe.xMin, inputframe.yMax, xtmin2, ytmax2);
  T.apply(inputframe.xMax, inputframe.yMin, xtmax2, ytmin2);
  Frame fr2(min(xtmin2,xtmax2), min(ytmin2,ytmax2), 
	    max(xtmin2,xtmax2), max(ytmin2,ytmax2));

  if (W == SmallFrame) return fr1*fr2;
  return fr1+fr2;
}

