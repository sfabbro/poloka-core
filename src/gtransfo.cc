#include "gtransfo.h"
#include <iostream>
#include <iomanip>
#include <iterator> /* for ostream_iterator */
#include "frame.h"

bool IsIdentity(const Gtransfo *a_transfo)
{ return (dynamic_cast<const GtransfoIdentity*>(a_transfo) != NULL);}



/********* Gtransfo ***********************/

Gtransfo* Gtransfo::ReduceCompo(const Gtransfo *Right) const 
{// by default no way to compose
  if (Right) {} // avoid a warning
  return NULL;
}

#ifdef USE_ROOT
#include <TFile.h>
#include <TKey.h> 


Gtransfo *GtransfoRead(const char* FileName, const char* ObjectName)
{
  TFile tfile(FileName);
  TKey *tkey  = NULL;
  if (ObjectName)
    {
      tkey = tfile.FindKey(ObjectName);
      if (!tkey)
	cerr << " Cannot find a Gtransfo named " << ObjectName << " in " <<  FileName  << endl;    
    }
  else
    {
      TIter nextkey(tfile.GetListOfKeys());
      tkey = (TKey*) nextkey();
      if (!tkey)
      cerr << " Cannot find anything in " <<  FileName  << endl;    
    }
  if (!tkey)
    {
      return NULL;
    }
  Gtransfo *t = dynamic_cast<Gtransfo*>(tkey->ReadObj());
  return t;
  // tfile.Close() called by destructor
}

int Gtransfo::root_write(const char*FileName, const char *FileStatus, const char* ObjectName)
{
  TFile tfile(FileName,FileStatus); 
  int count = Write(ObjectName);
  tfile.Close();
  return count;
}

#endif /* USE_ROOT */

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
  Der.a11 = (xp-xp0)/Step;
  Der.a21 = (yp-yp0)/Step;
  apply(x, y + Step, xp, yp);
  Der.a12 = (xp-xp0)/Step;
  Der.a22 = (yp-yp0)/Step;
  Der.dx = 0;
  Der.dy = 0;
}  

GtransfoLin Gtransfo::LinearApproximation(const Point &Where, 
					  const double Step) const
{
  Point outWhere = apply(Where);
  GtransfoLin der;
  Derivative(Where, der, Step);
  return GtransfoLinShift(outWhere.x, outWhere.y)*der*GtransfoLinShift(-Where.x, -Where.y);
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
     we wnat to compute M*V*tp(M)
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

// dummy implementation for derived classes which have Npar ==0;
void Gtransfo::ParamDerivatives(const Point &Where, double *Derivatives) const
{
  if (&Where && Derivatives) {}; // warning killer
  if (Npar() == 0)
    std::cerr << "Gtransfo::ParamDerivatives should never be called" 
	      << std::endl; 
}


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
  return 0; 
}

double &Gtransfo::ParamRef(const int i)
{
  if (i) {} // warning killer;
  return *(double *)(NULL); 
}



ostream & operator << (ostream &stream, const Gtransfo & T)
           {T.dump(stream); return stream;}


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
    if (&Region || &Precision ) {} //
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
// does not need to be public. Invoked  by Compose(Left,Right)

#include "gtransfocomposition.h"



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

Gtransfo *GtransfoCompose(const Gtransfo *Left, const Gtransfo *Right)
{
  /* is Right Identity ? if Left is Identity , GtransfoIdentity::ReduceCompo does the right job */
  if (IsIdentity(Right))
    {
      return Left->Clone();
    }
  /* try to use the ReduceCompo method of "this" */
  Gtransfo *composition = Left->ReduceCompo(Right);
  /* failure : just build a Composition that pipelines "this" and "Right" */
  if (composition == NULL) return new GtransfoComposition(Left,Right);
  else return composition;
}


// just a speed up, to avoid useless numerical derivation.
void GtransfoIdentity::Derivative(const Point &Where, 
				  GtransfoLin &Derivative, 
				  const double Step) const
{
  if (Step  || &Where) {} // warning killer
  Derivative.a11 = 1; Derivative.a12 = 0;
  Derivative.a21 = 0; Derivative.a22 = 1;
  Derivative.dx = Derivative.dy = 0;
}


GtransfoLin GtransfoIdentity::LinearApproximation(const Point &Where, const double Step) const
{
  if (Step  || &Where) {} // warning killer
  GtransfoLin result;
  return result; // rely on default Gtransfolin constructor;
}

/**************** GtransfoLin ***************************************/

GtransfoLin GtransfoLin::operator*(const  GtransfoLin &T2) const
{
  //   X' = Trans(dx,dy) * M(a11,a12,a21,a22) * X
  //   (T1,M1)*(T2,M2)*X = (T1,M1)*(T2+M2*X) = T1+M1*T2 + M1*M2*X = ((T1,M1)*T2, M1*M2) * X
GtransfoLin result;
apply(T2.dx, T2.dy, result.dx, result.dy);
result.a11 = this->a11*T2.a11 + this->a12*T2.a21;
result.a12 = this->a11*T2.a12 + this->a12*T2.a22;
result.a21 = this->a21*T2.a11 + this->a22*T2.a21;
result.a22 = this->a21*T2.a12 + this->a22*T2.a22;
return result;
}

void GtransfoLin::Derivative(const Point &Where, GtransfoLin &Derivative, 
			     const double Step) const
{
  if (Step || &Where) {} // "unused variable" warning killer
  Derivative = *this;
  Derivative.dx = 0;
  Derivative.dy = 0;
}


GtransfoLin GtransfoLin::LinearApproximation(const Point &Where, 
					     const double Step) const
{
  if (Step || &Where) {} // "unused variable" warning killer
  return *this;
}  

GtransfoLin GtransfoLin::invert() const
{
  //   following the above routine :
  //   (T1,M1) * (T2,M2) = 1 i.e (0,1) implies
  //   T1 = -M1 * T2
  //   M1 = M2^(-1)
  //

  /* this invert() routine does not apply to derived classes 
     (Quad and Cub transfos). With those derived classes,
     we can only trigger a run time error */

  if (dynamic_cast<const GtransfoQuad *>(this))
    {
      cerr << " GtransfoLin::invert is called with a non linear transfo"
	   << " no invertion possible, this is a very serious error " << endl;
      return GtransfoLin();
    }
  GtransfoLin result;
  double d = (a11*a22 - a12*a21);
  if (d == 0)
    {
      cerr << " trying to invert a singular transformation: a (nice) crash follows" << endl;
      dump(cerr);
    }
  result.a11 = a22/d;
  result.a12 = -a12/d;
  result.a22 = a11/d;
  result.a21 = -a21/d;
  result.dx  = -(result.a11*dx + result.a12*dy);
  result.dy  = -(result.a21*dx + result.a22*dy);
  return result;
} 

Gtransfo* GtransfoLin::InverseTransfo(const double Precision,
				      const Frame& Region) const
{
  if (&Precision || &Region) {} // warning killer
  return new GtransfoLin(this->invert());
}

void GtransfoLin::dump(ostream &stream) const
{
  ios::fmtflags  old_flags =  stream.flags(); 
  int oldprec = stream.precision();
  stream << resetiosflags(ios::scientific) ;
  stream << setiosflags(ios::fixed) ;
  stream << setprecision(8);
  stream << " newx = " << dx ;
  stream << " + " << a11 << "*x + " << a12 << "*y" << endl ;
  stream << " newy = " << dy ;
  stream << " + " << a21 << "*x + " << a22 << "*y" << endl ;
  stream << " Determinant = " << Determinant() << endl;
  stream.flags(old_flags);
  stream << setprecision(oldprec) ;
}



double GtransfoLin::ParamRef(const int i) const
{
  return this->*(Params[i].Item);
}

double &GtransfoLin::ParamRef(const int i)
{
  return this->*(Params[i].Item);
}

void GtransfoLin::ParamDerivatives(const Point &Where, double *Derivatives) const
{
  Derivatives[0] = Derivatives[9]  = 1;
  Derivatives[1] = Derivatives[10] = Where.x;
  Derivatives[2] = Derivatives[11] = Where.y;
  for (int i=3; i<9; ++i) Derivatives[i] = 0;
}

/* the ordering in the following array and the above routine are 
   linked together! */

GtransfoLin::LinParams GtransfoLin::Params[6] =
{{&GtransfoLin::dx}, {&GtransfoLin::a11}, {&GtransfoLin::a12},
 {&GtransfoLin::dy}, {&GtransfoLin::a21}, {&GtransfoLin::a22}};

  
#include "starmatch.h"


//used by fitting routines to set <x>=<y>=0, to limit numerical troubles.
static GtransfoLinShift ShiftToCenter(const StarMatchList &List)
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
  return GtransfoLinShift(-xav/count, -yav/count);
}
  

#include "vutils.h" /* for SymMatInv */

double GtransfoLin::fit(const StarMatchList &List)
{
#define NPAR 6
double a[NPAR][NPAR]; double b[NPAR]; double hx[NPAR]; double hy[NPAR];
// fit T defined by
//  x' = dx + a11*x + a12*y
//  y' = dy + a21*x + a22*y
// the ordering of the parameters is (a11,a12,a21,a22,dx,dy))
int i,j;
for (i=0; i < NPAR; i++) 
  { b[i] =0; hx[i] = 0. ; hy[i] = 0. ; for (j=0; j<NPAR; j++) a[i][j] = 0.;}

double sumr2 = 0; /* used to compute chi2 without relooping */
/* loop on pairs  and fill a and b */

GtransfoLinShift shift_to_center = ShiftToCenter(List);

StarMatchCIterator it;

for (it = List.begin(); it != List.end(); ++it)
  {
  const StarMatch &a_match = *it;
  Point point1 = shift_to_center.apply(a_match.point1);
  const Point &point2 = a_match.point2;
  hx[0] = point1.x;
  hx[1] = point1.y;
  hx[4] = 1.0;
  //
  hy[2] = point1.x;
  hy[3] = point1.y;
  hy[5] = 1.0;
  for (i=0; i< NPAR; i++)
    {
    b[i] += hx[i]*point2.x + hy[i] * point2.y;
    for (j=0; j<NPAR; j++) a[i][j] += hx[i]*hx[j] + hy[i]*hy[j];
    }
  sumr2 += point2.x*point2.x + point2.y*point2.y;
  }

/* symetrize the matrix */
//for (i=0; i<NPAR; i++) for (j=0; j<i; j++ ) a[i][j] = a[j][i];

/* solve */
if (!SymMatInv(&a[0][0],NPAR))   return -1;

double res[NPAR];
MatVec(&a[0][0], NPAR, NPAR, b, res);

a11 = res[0];
a12 = res[1];
a21 = res[2];
a22 = res[3]; 
dx  = res[4];
dy  = res[5];
 (*this) =  (*this) * shift_to_center;
/* compute chi2 */
return (sumr2 - ScalProd(res,b,NPAR));
}

#undef NPAR


#ifdef OLD_STUFF
double  GtransfoLin::fit(const StarMatchList &List)
{
int npairs = List.size();
if (List.size()< 3)
  {
  cerr << " GTransfoLin::fit  : trying to fit a linear transfo with only " << npairs << " matches " << endl;
  return -1;
  }


/* do the job twice in case of numerical problem the first time */
GtransfoLin firstLoop, secondLoop;
double chi2;
firstLoop.do_the_fit(chi2, List);

Gtransfo *globalTransfo = GtransfoCompose(&firstLoop,PriorTransfo);

/* second turn : secondLoop will be in principle Identity (if no serious numerical problem showed up)  */
secondLoop.do_the_fit( chi2, List, globalTransfo, PosteriorTransfo);
*this = secondLoop*firstLoop;
delete globalTransfo;
return chi2;
}
#endif /* OLD_STUFF */


double  GtransfoLinRot::fit(const StarMatchList &List)
{
int npairs = List.size();
if (List.size()< 3)
  {
  cerr << " GTransfoLinRot::fit : trying to fit a linear transfo with only " << npairs << " matches " << endl;
  return -1;
  }


#define NPAR 4
double a[NPAR][NPAR]; double b[NPAR]; double hx[NPAR]; double hy[NPAR];
// fit T defined by
//  x' = ax - by + c 
//  y' = bx + ay + d 
// the ordering of the parameters is (a,b,c,d)
int i,j;
for (i=0; i < NPAR; i++) { b[i] =0; for (j=0; j<NPAR; j++) a[i][j] = 0.;}

GtransfoLinShift shift_to_center = ShiftToCenter(List);

double sumr2 = 0; /* used to compute chi2 without relooping */
/* loop on pairs  and fill a and b */
StarMatchCIterator it;
for (it = List.begin(); it != List.end(); it++)
  {
  const StarMatch &a_match = *it;
  Point point1 = shift_to_center.apply(a_match.point1);
  const Point &point2 = a_match.point2;
  hx[0] = point1.x;
  hx[1] = -point1.y;
  hx[2] = 1.0;
  hx[3] = 0;

  hy[0] = point1.y;
  hy[1] = point1.x;
  hy[2] = 0.;
  hy[3] = 1.0;
  for (i=0; i< NPAR; i++)
    {
    b[i] += hx[i]*point2.x + hy[i] * point2.y;
    for (j=0; j<NPAR; j++) a[i][j] += hx[i]*hx[j] + hy[i]*hy[j];
    }
  sumr2 += point2.x*point2.x + point2.y*point2.y;
  }
/* symetrize the matrix */
//for (i=0; i<NPAR; i++) for (j=0; j<i; j++ ) a[i][j] = a[j][i];

/* solve */
if (!SymMatInv(&a[0][0],NPAR))
  cout << " could not invert in GtransfoLin::fit " << endl;

double res[NPAR];
 for (i=0; i<NPAR; i++) {res[i] = 0; for (j=0;j<NPAR;j++) res[i] += a[i][j]*b[j];}

a11 = res[0];
a21 = res[1];
dx = res[2];
dy = res[3];
a12 = -a21;
a22 = a11;
/* compute chi2 */
return (sumr2 - (res[0]*b[0]+res[1]*b[1]+res[2]*b[2]+res[3]*b[3]));
}
#undef NPAR


double  GtransfoLinShift::fit(const StarMatchList &List)
{
int npairs = List.size();
if (npairs < 3)
  {
    // cerr << " GTransfoLinShift::fit  : trying to fit a linear transfo with only " << npairs << " matches " << endl;
  return -1;
  }

// just fitting offsets: they are defined by the average value of delta_x and delta_y

double sumr2 = 0; /* used to compute chi2 without relooping */
/* loop on pairs  and fill */
double sDeltaX = 0;
double sDeltaY = 0;
double s1 =0;
StarMatchCIterator it;
for (it = List.begin(); it != List.end(); it++)
  {
  double deltax = it->point2.x - it->point1.x;
  double deltay = it->point2.y - it->point1.y;
  sDeltaX += deltax;
  sDeltaY += deltay;
  s1 += 1;
  sumr2 += deltax*deltax + deltay*deltay;
  }

dx = sDeltaX /s1;
dy = sDeltaY /s1;
return (sumr2 + sDeltaX*dx + sDeltaY*dy);
}

Gtransfo *GtransfoLin::ReduceCompo(const Gtransfo *Right) const
{
  // if (IsIdentity(Right))  return this->Clone(); done in Compose(Left,Right)
  // try cub:
  const GtransfoCub *pc = dynamic_cast<const GtransfoCub*>(Right);
  if (pc != NULL) return new GtransfoCub( (*this) * (*pc) );
  // try quad
  const GtransfoQuad *pq = dynamic_cast<const GtransfoQuad*>(Right);
  if (pq != NULL) return new GtransfoQuad((*this) * (*pq));
  //! try lin
  const GtransfoLin *pl = dynamic_cast<const GtransfoLin*>(Right);
  if (pl) return new GtransfoLin((*this) * (*pl));
  return NULL;
}


#include <math.h>

GtransfoLinRot::GtransfoLinRot(const double AngleRad, const Point *Center, 
			       const double ScaleFactor)
{
double c = ScaleFactor*cos(AngleRad);
double s = ScaleFactor*sin(AngleRad);
a11= a22 = c;
a21 = s;
a12 = -s;

// we want that the center does not move : T+M*C = C ==> T = C - M*C
Point a_point(0.,0.);
if (Center) a_point = *Center;

dx = dy = 0;
apply(a_point.x, a_point.y, dx, dy); // compute M*C
dx = a_point.x - dx; dy = a_point.y - dy;
}




//**********************************************************************************
//***************************** GtransfoQuad ***************************************


#include "vutils.h" /* for SymMatInv */
#include "starmatch.h"

GtransfoQuad::GtransfoQuad(const GtransfoLin & Lin)
{
  identity();
  const GtransfoCub *c = dynamic_cast<const GtransfoCub *>(&Lin);
  if (c) {*this = *c; return;}
  const GtransfoQuad *q = dynamic_cast<const GtransfoQuad *>(&Lin);
  if (q) { GtransfoQuad &tq = *this;  tq = *q; return;}
  GtransfoLin &tl = *this;
  tl = Lin;
}


Gtransfo *GtransfoQuad::ReduceCompo(const Gtransfo *Right) const
{
  //! try lin
  const GtransfoLin *pl = dynamic_cast<const GtransfoLin*>(Right);
  if (pl && pl->Degree() == 1) return new GtransfoQuad((*this) * (*pl));
  return NULL;
}


void GtransfoQuad::dump(ostream &stream) const
{
  ios::fmtflags  old_flags =  stream.flags(); 
  int oldprec = stream.precision();
  stream << resetiosflags(ios::scientific) ;
  stream << setiosflags(ios::fixed) ;
  stream << setprecision(8);
  stream << " newx = " << dx ;
  stream << " + " << a11 << "*x + " << a12 
	 << "*y +" << a1x2 << "*x^2 +" 
	 << a1xy << "*x*y +" << a1y2 
	 << "*y^2" << endl ;
  stream << " newy = " << dy ;
  stream << " + " << a21 << "*x + " << a22 << "*y +" << 
    a2x2 << "*x^2 +" << a2xy << "*x*y +" << a2y2 << "*y^2" << endl ;
  stream << " Linear Determinant = " << a11*a22-a21*a12 << endl ;
  stream.flags(old_flags);
  stream << setprecision(oldprec);
}


double GtransfoQuad::ParamRef(const int i) const
{
  return this->*(Params[i].Item);
}

double &GtransfoQuad::ParamRef(const int i)
{
  return this->*(Params[i].Item);
}

void GtransfoQuad::ParamDerivatives(const Point &Where, double *Derivatives) const
{
  Derivatives[0] = Derivatives[18]  = 1;
  Derivatives[1] = Derivatives[19] = Where.x;
  Derivatives[2] = Derivatives[20] = Where.y;
  Derivatives[3] = Derivatives[21] = Where.x * Where.x;
  Derivatives[4] = Derivatives[22] = Where.x * Where.y;
  Derivatives[5] = Derivatives[23] = Where.y * Where.y;
  for (int i=6; i<24; ++i) Derivatives[i] = 0;
}

/* the ordering in the following array and the above routine are 
   linked together! */

GtransfoQuad::QuadParams GtransfoQuad::Params[12] =
{{&GtransfoQuad::dx},   {&GtransfoQuad::a11},  {&GtransfoQuad::a12},
 {&GtransfoQuad::a1x2}, {&GtransfoQuad::a1xy}, {&GtransfoQuad::a1y2},
 {&GtransfoQuad::dy},   {&GtransfoQuad::a21},  {&GtransfoQuad::a22},
 {&GtransfoQuad::a2x2}, {&GtransfoQuad::a2xy}, {&GtransfoQuad::a2y2},};


double  GtransfoQuad::fit(const StarMatchList &List)
{
  int npairs = List.size();
  if (List.size()< 3)
    {
    cerr << " GtransfoQuad::fit : trying to fit a quadratic transfo with only " << npairs << " matches " << endl;
    return -1;
    }
  
#define NPAR 6
  double a[NPAR][NPAR];  double bx[NPAR]; double by[NPAR]; 
  
  // fit T defined by
  //  x' = dx + a11*x + a12*y + a1x2*x^2 + a1xy*x*y + a1y2*y^2
  //  y' = dy + a21*x + a22*y + a2x2*x^2 + a2xy*x*y + a2y2*y^2
  
  double sumr2 = 0; /* used to compute chi2 without relooping */
  
  // loop on pairs and fill in the 2 Gram matrices ax and ay 
  // hx = (1 x y x^2 xy y^2) same for hy 
  // b is hx * point2
  
  memset(a,0,sizeof(double)*NPAR*NPAR);
  memset(bx,0,sizeof(double)*NPAR);
  memset(by,0,sizeof(double)*NPAR);
  
  GtransfoLinShift shift_to_center = ShiftToCenter(List);

  for (StarMatchCIterator it = List.begin(); it != List.end(); ++it)
    {
      const StarMatch &a_match = *it;
      Point point1 = shift_to_center.apply(a_match.point1);
      const Point &point2 = a_match.point2;
      double h[NPAR];
      h[0] = 1.0;
      h[1] = point1.x;
      h[2] = point1.y;
      h[3] = point1.x * point1.x;
      h[4] = point1.x * point1.y;
      h[5] = point1.y * point1.y;

      for (int i=0; i< NPAR; ++i)
	{
	  bx[i] += h[i] * point2.x;
	  by[i] += h[i] * point2.y;
	  for (int j=0; j<NPAR; ++j) a[i][j] += h[i]*h[j];
	}
      sumr2 += point2.x*point2.x + point2.y*point2.y;
    }

  // now invert the symetric matrices
  if (!SymMatInv(&a[0][0],NPAR))
    {
      cout << " could not invert the matrix  in GtransfoQuad::fit " << endl;
     return -1;
    }

  double paramx[NPAR];
  double paramy[NPAR];
  for (int i=0; i<NPAR; ++i) 
    {
      paramx[i] = 0; 
      paramy[i] = 0;
      for (int j=0; j<NPAR; ++j) 
	{
	  paramx[i] += a[i][j]*bx[j];
	  paramy[i] += a[i][j]*by[j];
	}
    }
  
  // fill in the estimated parameters
  dx  = paramx[0];
  dy  = paramy[0];
  a11 = paramx[1];
  a21 = paramy[1];
  a12 = paramx[2];
  a22 = paramy[2]; 
  a1x2 = paramx[3];
  a2x2 = paramy[3];
  a1xy = paramx[4];
  a2xy = paramy[4];
  a1y2 = paramx[5];
  a2y2 = paramy[5];

  // compute and return chi2 for the two fits
  *this = *this * shift_to_center;

  return (sumr2 - ScalProd(paramx,bx,NPAR) - ScalProd(paramy,by,NPAR));

#undef NPAR
 
}


// need to overload GtransfoLin::InverseTransfo !
Gtransfo* GtransfoQuad::InverseTransfo(const double Precision,
				      const Frame& Region) const
{
  return Gtransfo::InverseTransfo(Precision, Region);
}

// same as above
void GtransfoQuad::Derivative(const Point &Where, 
			      GtransfoLin &Derivative, 
			      const double Step) const
{
  Gtransfo::Derivative(Where, Derivative, Step);
}

// overload GtransfoLin routine
GtransfoLin GtransfoQuad::LinearApproximation(const Point &Where, 
					      const double Step) const
{
  return Gtransfo::LinearApproximation(Where, Step);
}


/* 
we use the same routine to perform GtransfoLin*GtransfoQuad
and GtransfoQuad*GtransfoLin: the product of 2 GtransfoQuad 
truncated to 2nd order. no need to take temporary copies even for Q=Q*L; */

GtransfoQuad operator*(const GtransfoLin &L, const GtransfoQuad &R)
{
  GtransfoQuad temp(L);
  return temp.truncated_product(R);
}

GtransfoQuad GtransfoQuad::operator*(const  GtransfoLin &R) const
{
  return truncated_product(R);
}


GtransfoQuad GtransfoQuad::truncated_product(const GtransfoQuad &R) const
{
  //computed using Maple...
  GtransfoQuad result;
  result.dx = dx + a11*R.dx + a12*R.dy + a1x2*R.dx*R.dx + a1xy*R.dx*R.dy + a1y2*R.dy*R.dy;

  result.a11 =  a11*R.a11 + a12*R.a21 + 2.0*a1x2*R.dx*R.a11 + a1xy*R.a11*R.dy + 
    a1xy*R.dx*R.a21 + 2.0*a1y2*R.dy*R.a21;

  result.a12 =  a11*R.a12 + a12*R.a22 + 2.0*a1x2*R.dx*R.a12 + 
    a1xy*R.a12*R.dy + a1xy*R.dx*R.a22 + 2.0*a1y2*R.dy*R.a22;

  result.a1x2 = a11*R.a1x2 + a12*R.a2x2 + a1x2*R.a11*R.a11 + 2.0*a1x2*R.dx*R.a1x2 + 
    a1xy*R.a1x2*R.dy + a1xy*R.a11*R.a21 + a1xy*R.dx*R.a2x2 + a1y2*R.a21*R.a21 + 
    2.0*a1y2*R.dy*R.a2x2;

  result.a1xy = a11*R.a1xy + a12*R.a2xy + 2.0*a1x2*R.a12*R.a11 + 2.0*a1x2*R.dx*R.a1xy + 
    a1xy*R.a1xy*R.dy + a1xy*R.a11*R.a22 + a1xy*R.a12*R.a21 + a1xy*R.dx*R.a2xy + 
    2.0*a1y2*R.a22*R.a21 + 2.0*a1y2*R.dy*R.a2xy;


  result.a1y2 = a11*R.a1y2 + a12*R.a2y2 + a1x2*R.a12*R.a12 + 2.0*a1x2*R.dx*R.a1y2 + 
    a1xy*R.a1y2*R.dy + a1xy*R.a12*R.a22 + a1xy*R.dx*R.a2y2 + a1y2*R.a22*R.a22 + 
    2.0*a1y2*R.dy*R.a2y2;


  result.dy = dy + a21*R.dx + a22*R.dy + a2x2*R.dx*R.dx + a2xy*R.dx*R.dy + a2y2*R.dy*R.dy;

  result.a21 =  a21*R.a11 + a22*R.a21 + 2.0*a2x2*R.dx*R.a11 + a2xy*R.a11*R.dy + 
    a2xy*R.dx*R.a21 + 2.0*a2y2*R.dy*R.a21;

   result.a22 = a21*R.a12 + a22*R.a22 + 2.0*a2x2*R.dx*R.a12 + a2xy*R.a12*R.dy + 
     a2xy*R.dx*R.a22 + 2.0*a2y2*R.dy*R.a22;


   result.a2x2 = a21*R.a1x2 + a22*R.a2x2 + a2x2*R.a11*R.a11 + 2.0*a2x2*R.dx*R.a1x2 
     + a2xy*R.a1x2*R.dy + a2xy*R.a11*R.a21 + a2xy*R.dx*R.a2x2 + a2y2*R.a21*R.a21 + 2.0*a2y2*R.dy*R.a2x2;

   result.a2xy = a21*R.a1xy + a22*R.a2xy + 2.0*a2x2*R.a12*R.a11 + 2.0*a2x2*R.dx*R.a1xy + 
     a2xy*R.a1xy*R.dy + a2xy*R.a11*R.a22 + a2xy*R.a12*R.a21 + a2xy*R.dx*R.a2xy + 
     2.0*a2y2*R.a22*R.a21 + 2.0*a2y2*R.dy*R.a2xy;

   result.a2y2 = a21*R.a1y2 + a22*R.a2y2 + a2x2*R.a12*R.a12 + 2.0*a2x2*R.dx*R.a1y2 + 
     a2xy*R.a1y2*R.dy + a2xy*R.a12*R.a22 + a2xy*R.dx*R.a2y2 + a2y2*R.a22*R.a22 + 
     2.0*a2y2*R.dy*R.a2y2;

  return result;
}



#ifdef STORAGE
GtransfoQuad GtransfoQuad::truncated_product(const GtransfoQuad &R) const
{
  //computed using Maple...
  GtransfoQuad result;
  result.dx = dx+a11*R.dx+a1y2*R.dy*R.dy+a1x2*R.dx*R.dx+a12*R.dy +a1xy*R.dx*R.dy;
  result.a11 = (a12*R.a21+2.0*a1x2*R.dx*R.a11+2.0*a1y2*R.dy*R.a21
		+a1xy*R.dx*R.a21+a1xy*R.a11*R.dy+a11*R.a11);
  result.a12 = (a12*R.a22+2.0*a1x2*R.dx*R.a12+2.0*a1y2*R.dy*R.a22+a1xy*R.dx*R.a22
		+a1xy*R.a12*R.dy+a11*R.a12);
  result.a1x2 = (a1xy*R.a1x2*R.dy+a1xy*R.a11*R.a21+a1xy*R.dx*R.a2x2+a11*R.a1x2+
		 a1x2*(2.0*R.dx*R.a1x2+R.a11*R.a11)+a12*R.a2x2+
		 a1y2*(2.0*R.dy*R.a2x2+R.a21*R.a21));

  result.a1xy = (a1xy*R.dx*R.a2xy+a1x2*(2.0*R.dx*R.a1xy+2.0*R.a12*R.a11)
		 +a12*R.a2xy+a1xy*R.a12*R.a21+a1xy*R.a11*R.a22
		 +a1y2*(2.0*R.dy*R.a2xy+2.0*R.a22*R.a21)+a11*R.a1xy+a1xy*R.a1xy*R.dy);

  result.a1y2 = (a12*R.a2y2+a1y2*(2.0*R.dy*R.a2y2+R.a22*R.a22)
		 +a11*R.a1y2+a1x2*(2.0*R.dx*R.a1y2+R.a12*R.a12)+
		 a1xy*R.a1y2*R.dy+a1xy*R.a12*R.a22+a1xy*R.dx*R.a2y2);

  result.dy = dy+a21*R.dx+a2y2*R.dy*R.dy+a2x2*R.dx*R.dx+a22*R.dy+a2xy*R.dx*R.dy;
  result.a21 = (a22*R.a21+2.0*a2x2*R.dx*R.a11+2.0*a2y2*R.dy*R.a21
		+a2xy*R.dx*R.a21+a2xy*R.a11*R.dy+a21*R.a11);
  result.a22 = (a22*R.a22+2.0*a2x2*R.dx*R.a12+2.0*a2y2*R.dy*R.a22
		+a2xy*R.dx*R.a22+a2xy*R.a12*R.dy+a21*R.a12);
  result.a2x2 = (a2xy*R.a1x2*R.dy+a2xy*R.a11*R.a21+a2xy*R.dx*R.a2x2+a21*R.a1x2
		 +a2x2*(2.0*R.dx*R.a1x2+R.a11*R.a11)+
		 a22*R.a2x2+a2y2*(2.0*R.dy*R.a2x2+R.a21*R.a21));
  result.a2xy = (a2xy*R.dx*R.a2xy+a2x2*(2.0*R.dx*R.a1xy+2.0*R.a12*R.a11)+
		 a22*R.a2xy+a2xy*R.a12*R.a21+a2xy*R.a11*R.a22+
		 a2y2*(2.0*R.dy*R.a2xy+2.0*R.a22*R.a21)+a21*R.a1xy+a2xy*R.a1xy*R.dy);
  result.a2y2 = (a22*R.a2y2+a2y2*(2.0*R.dy*R.a2y2+R.a22*R.a22)+
		 a21*R.a1y2+a2x2*(2.0*R.dx*R.a1y2+R.a12*R.a12)+
		 a2xy*R.a1y2*R.dy+a2xy*R.a12*R.a22+a2xy*R.dx*R.a2y2);
  return result;
}
#endif

//**********************************************************************************
//***************************** GtransfoCub ***************************************


#include "vutils.h" /* for SymMatInv */
#include "starmatch.h"

GtransfoCub::GtransfoCub(const GtransfoLin & Lin)
{
  identity(); // sets to zero almost everything;
  const GtransfoCub *c = dynamic_cast<const GtransfoCub *>(&Lin);
  if (c) {*this = *c; return;}
  const GtransfoQuad *q = dynamic_cast<const GtransfoQuad *>(&Lin);
  if (q) { GtransfoQuad &tq = *this;  tq = *q; return;}
  GtransfoLin &tl = *this;
  tl = Lin;
}

// isn't this routine useless?, the one above does the right job.
GtransfoCub::GtransfoCub(const GtransfoQuad & Quad)
{
  identity(); // sets almost everything to zero
  const GtransfoCub *c = dynamic_cast<const GtransfoCub *>(&Quad);
  if (c) { *this = *c; return;}
  GtransfoQuad &tq = *this; tq = Quad;
}

Gtransfo *GtransfoCub::ReduceCompo(const Gtransfo *Right) const
{
  //! try lin
  const GtransfoLin *pl = dynamic_cast<const GtransfoLin*>(Right);
  if (pl && pl->Degree() == 1) return new GtransfoCub((*this) * (*pl));
  return NULL;
}

GtransfoCub operator*(const GtransfoLin &L, const GtransfoCub &R)
{
  //computed by Maple...
  GtransfoCub result;

  result.dx =  L.dx + L.a11*R.dx + L.a12*R.dy;
  result.a11 = L.a11*R.a11 + L.a12*R.a21;
  result.a12 = L.a11*R.a12 + L.a12*R.a22;
  result.a1x2 = L.a11*R.a1x2 + L.a12*R.a2x2;
  result.a1xy = L.a11*R.a1xy + L.a12*R.a2xy;
  result.a1y2 = L.a11*R.a1y2 + L.a12*R.a2y2;
  result.a1x3 = L.a11*R.a1x3 + L.a12*R.a2x3;
  result.a1x2y = L.a11*R.a1x2y + L.a12*R.a2x2y;
  result.a1xy2  = L.a11*R.a1xy2 + L.a12*R.a2xy2;
  result.a1y3 = L.a11*R.a1y3 + L.a12*R.a2y3;
  
  result.dy = L.dy + L.a21*R.dx + L.a22*R.dy;
  result.a21 = L.a21*R.a11 + L.a22*R.a21;
  result.a22 = L.a21*R.a12 + L.a22*R.a22;
  result.a2x2 = L.a21*R.a1x2 + L.a22*R.a2x2;
  result.a2xy = L.a21*R.a1xy + L.a22*R.a2xy;
  result.a2y2 = L.a21*R.a1y2 + L.a22*R.a2y2;
  result.a2x3 = L.a21*R.a1x3 + L.a22*R.a2x3;
  result.a2x2y = L.a21*R.a1x2y + L.a22*R.a2x2y;
  result.a2xy2 = L.a21*R.a1xy2 + L.a22*R.a2xy2;
  result.a2y3 = L.a21*R.a1y3 + L.a22*R.a2y3;

  return result;
}

GtransfoCub operator*(const GtransfoCub &L, const GtransfoLin &R)
{
  //computed by Maple...
  GtransfoCub result;
  result.dx = L.dx + L.a11*R.dx + L.a12*R.dy + L.a1x2*R.dx*R.dx + L.a1xy*R.dx*R.dy + L.a1y2*R.dy*R.dy + 
    L.a1x3*R.dx*R.dx*R.dx + L.a1x2y*R.dx*R.dx*R.dy + L.a1xy2*R.dx*R.dy*R.dy + L.a1y3*R.dy*R.dy*R.dy;

  result.a11 = L.a11*R.a11 + L.a12*R.a21 + 2.0*L.a1x2*R.dx*R.a11 + L.a1xy*R.a11*R.dy + L.a1xy*R.dx*R.a21 + 2.0*
L.a1y2*R.dy*R.a21 + 3.0*L.a1x3*R.dx*R.dx*R.a11 + 2.0*L.a1x2y*R.dx*R.dy*R.a11 + L.a1x2y*R.dx*R.dx*R.a21 + 
L.a1xy2*R.a11*R.dy*R.dy + 2.0*L.a1xy2*R.dx*R.dy*R.a21 + 3.0*L.a1y3*R.dy*R.dy*R.a21;

  result.a12 = L.a11*R.a12 + L.a12*R.a22 + 2.0*L.a1x2*R.dx*R.a12 + L.a1xy*R.a12*R.dy + L.a1xy*R.dx*R.a22 + 2.0*
L.a1y2*R.dy*R.a22 + 3.0*L.a1x3*R.dx*R.dx*R.a12 + 2.0*L.a1x2y*R.dx*R.dy*R.a12 + L.a1x2y*R.dx*R.dx*R.a22 + 
L.a1xy2*R.a12*R.dy*R.dy + 2.0*L.a1xy2*R.dx*R.dy*R.a22 + 3.0*L.a1y3*R.dy*R.dy*R.a22;

  result.a1x2 = L.a1x2*R.a11*R.a11 + L.a1xy*R.a11*R.a21 + L.a1y2*R.a21*R.a21 + 3.0*L.a1x3*R.dx*R.a11*R.a11 + 
L.a1x2y*R.a11*R.a11*R.dy + 2.0*L.a1x2y*R.dx*R.a21*R.a11 + 2.0*L.a1xy2*R.a11*R.dy*R.a21 + L.a1xy2*R.dx*
R.a21*R.a21 + 3.0*L.a1y3*R.dy*R.a21*R.a21;

  result.a1xy = 2.0*L.a1x2*R.a12*R.a11 + L.a1xy*R.a11*R.a22 + L.a1xy*R.a12*R.a21 + 2.0*L.a1y2*R.a22*R.a21 + 
6.0*L.a1x3*R.dx*R.a11*R.a12 + 2.0*L.a1x2y*R.a12*R.dy*R.a11 + 2.0*L.a1x2y*R.dx*R.a22*R.a11 + 2.0*
L.a1x2y*R.dx*R.a21*R.a12 + 2.0*L.a1xy2*R.a11*R.dy*R.a22 + 2.0*L.a1xy2*R.a12*R.dy*R.a21 + 2.0*L.a1xy2*
R.dx*R.a22*R.a21 + 6.0*L.a1y3*R.dy*R.a21*R.a22;

  result.a1y2 = L.a1x2*R.a12*R.a12 + L.a1xy*R.a12*R.a22 + L.a1y2*R.a22*R.a22 + 3.0*L.a1x3*R.dx*R.a12*R.a12 + 
L.a1x2y*R.a12*R.a12*R.dy + 2.0*L.a1x2y*R.dx*R.a22*R.a12 + 2.0*L.a1xy2*R.a12*R.dy*R.a22 + L.a1xy2*R.dx*
R.a22*R.a22 + 3.0*L.a1y3*R.dy*R.a22*R.a22;

  result.a1x3 = L.a1x3*R.a11*R.a11*R.a11 + L.a1x2y*R.a11*R.a11*R.a21 + L.a1xy2*R.a11*R.a21*R.a21 + L.a1y3*
R.a21*R.a21*R.a21;

  result.a1x2y = 3.0*L.a1x3*R.a12*R.a11*R.a11 + L.a1x2y*R.a11*R.a11*R.a22 + 2.0*L.a1x2y*R.a12*R.a21*
R.a11 + 2.0*L.a1xy2*R.a11*R.a22*R.a21 + L.a1xy2*R.a12*R.a21*R.a21 + 3.0*L.a1y3*R.a22*R.a21*R.a21;

  result.a1xy2 = 3.0*L.a1x3*R.a12*R.a12*R.a11 + 2.0*L.a1x2y*R.a12*R.a22*R.a11 + L.a1x2y*R.a12*R.a12*
R.a21 + L.a1xy2*R.a11*R.a22*R.a22 + 2.0*L.a1xy2*R.a12*R.a22*R.a21 + 3.0*L.a1y3*R.a22*R.a22*R.a21;

  result.a1y3 = L.a1x3*R.a12*R.a12*R.a12 + L.a1x2y*R.a12*R.a12*R.a22 + L.a1xy2*R.a12*R.a22*R.a22 + L.a1y3*
R.a22*R.a22*R.a22;

  result.dy = L.dy + L.a21*R.dx + L.a22*R.dy + L.a2x2*R.dx*R.dx + L.a2xy*R.dx*R.dy + L.a2y2*R.dy*R.dy + L.a2x3*R.dx*
R.dx*R.dx + L.a2x2y*R.dx*R.dx*R.dy + L.a2xy2*R.dx*R.dy*R.dy + L.a2y3*R.dy*R.dy*R.dy;

  result.a21 = L.a21*R.a11 + L.a22*R.a21 + 2.0*L.a2x2*R.dx*R.a11 + L.a2xy*R.a11*R.dy + L.a2xy*R.dx*R.a21 + 2.0*
L.a2y2*R.dy*R.a21 + 3.0*L.a2x3*R.dx*R.dx*R.a11 + 2.0*L.a2x2y*R.dx*R.dy*R.a11 + L.a2x2y*R.dx*R.dx*R.a21 + 
L.a2xy2*R.a11*R.dy*R.dy + 2.0*L.a2xy2*R.dx*R.dy*R.a21 + 3.0*L.a2y3*R.dy*R.dy*R.a21;

  result.a22 = L.a21*R.a12 + L.a22*R.a22 + 2.0*L.a2x2*R.dx*R.a12 + L.a2xy*R.a12*R.dy + L.a2xy*R.dx*R.a22 + 2.0*
L.a2y2*R.dy*R.a22 + 3.0*L.a2x3*R.dx*R.dx*R.a12 + 2.0*L.a2x2y*R.dx*R.dy*R.a12 + L.a2x2y*R.dx*R.dx*R.a22 + 
L.a2xy2*R.a12*R.dy*R.dy + 2.0*L.a2xy2*R.dx*R.dy*R.a22 + 3.0*L.a2y3*R.dy*R.dy*R.a22;

  result.a2x2 = L.a2x2*R.a11*R.a11 + L.a2xy*R.a11*R.a21 + L.a2y2*R.a21*R.a21 + 3.0*L.a2x3*R.dx*R.a11*R.a11 + 
L.a2x2y*R.a11*R.a11*R.dy + 2.0*L.a2x2y*R.dx*R.a21*R.a11 + 2.0*L.a2xy2*R.a11*R.dy*R.a21 + L.a2xy2*R.dx*
R.a21*R.a21 + 3.0*L.a2y3*R.dy*R.a21*R.a21;

  result.a2xy = 2.0*L.a2x2*R.a12*R.a11 + L.a2xy*R.a11*R.a22 + L.a2xy*R.a12*R.a21 + 2.0*L.a2y2*R.a22*R.a21 + 
6.0*L.a2x3*R.dx*R.a11*R.a12 + 2.0*L.a2x2y*R.a12*R.dy*R.a11 + 2.0*L.a2x2y*R.dx*R.a22*R.a11 + 2.0*
L.a2x2y*R.dx*R.a21*R.a12 + 2.0*L.a2xy2*R.a11*R.dy*R.a22 + 2.0*L.a2xy2*R.a12*R.dy*R.a21 + 2.0*L.a2xy2*
R.dx*R.a22*R.a21 + 6.0*L.a2y3*R.dy*R.a21*R.a22;

  result.a2y2 = L.a2x2*R.a12*R.a12 + L.a2xy*R.a12*R.a22 + L.a2y2*R.a22*R.a22 + 3.0*L.a2x3*R.dx*R.a12*R.a12 + 
L.a2x2y*R.a12*R.a12*R.dy + 2.0*L.a2x2y*R.dx*R.a22*R.a12 + 2.0*L.a2xy2*R.a12*R.dy*R.a22 + L.a2xy2*R.dx*
R.a22*R.a22 + 3.0*L.a2y3*R.dy*R.a22*R.a22;

  result.a2x3 = L.a2x3*R.a11*R.a11*R.a11 + L.a2x2y*R.a11*R.a11*R.a21 + L.a2xy2*R.a11*R.a21*R.a21 + L.a2y3*
R.a21*R.a21*R.a21;

  result.a2x2y = 3.0*L.a2x3*R.a12*R.a11*R.a11 + L.a2x2y*R.a11*R.a11*R.a22 + 2.0*L.a2x2y*R.a12*R.a21*
R.a11 + 2.0*L.a2xy2*R.a11*R.a22*R.a21 + L.a2xy2*R.a12*R.a21*R.a21 + 3.0*L.a2y3*R.a22*R.a21*R.a21;

  result.a2xy2 = 3.0*L.a2x3*R.a12*R.a12*R.a11 + 2.0*L.a2x2y*R.a12*R.a22*R.a11 + L.a2x2y*R.a12*R.a12*
R.a21 + L.a2xy2*R.a11*R.a22*R.a22 + 2.0*L.a2xy2*R.a12*R.a22*R.a21 + 3.0*L.a2y3*R.a22*R.a22*R.a21;

  result.a2y3 =  L.a2x3*R.a12*R.a12*R.a12 + L.a2x2y*R.a12*R.a12*R.a22 + L.a2xy2*R.a12*R.a22*R.a22 + L.a2y3*
R.a22*R.a22*R.a22;

  return result;
}


void GtransfoCub::dump(ostream &stream) const
{
  ios::fmtflags  old_flags =  stream.flags(); 
  int oldprec = stream.precision();
  stream << resetiosflags(ios::scientific) ;
  stream << setiosflags(ios::fixed) ;
  stream << setprecision(8);
  stream << " newx = " << dx ;
  /*
    stream << resetiosflags(ios::fixed) ;
    stream << setiosflags(ios::scientific) ;
    stream << setprecision(8);
  */
  stream << " + " << a11 << "*x + " << a12 
	 << "*y +" << a1x2 << "*x^2 +" << a1xy 
	 << "*x*y +" << a1y2 << "*y^2 +" << a1x3 
	 << "*x^3 +" << a1x2y << "*x^2*y +" 
	 << a1xy2 << "*x*y^2 + "<< a1y3 
	 << "*y^3"<<endl ;
  /* 
     stream << resetiosflags(ios::scientific) ;
     stream << setiosflags(ios::fixed) ;
     stream << setprecision(8);	   
  */
  stream << " newy = " << dy ;
  /* 
     stream << resetiosflags(ios::fixed) ;
     stream << setiosflags(ios::scientific) ;
     stream << setprecision(8);
  */
  stream << " + " << a21 << "*x + " << a22 
	 << "*y +" << a2x2 << "*x^2 +" << a2xy 
	 << "*x*y +" << a2y2 << "*y^2 +" << a2x3 
	 << "*x^3 +" << a2x2y << "*x^2*y +" 
	 << a2xy2 << "*x*y^2 + "<< a2y3 << "*y^3"
	 <<endl ;
  stream << " Linear Determinant = " << a11*a22-a21*a12 << endl;
  stream.flags(old_flags);
  stream << setprecision(oldprec) ;
}


double GtransfoCub::ParamRef(const int i) const
{
  return this->*(Params[i].Item);
}

double &GtransfoCub::ParamRef(const int i)
{
  return this->*(Params[i].Item);
}

void GtransfoCub::ParamDerivatives(const Point &Where, double *Derivatives) const
{
  Derivatives[0] = Derivatives[30]  = 1;
  Derivatives[1] = Derivatives[31] = Where.x;
  Derivatives[2] = Derivatives[32] = Where.y;
  double x2 = Where.x *Where.x;
  Derivatives[3] = Derivatives[33] = x2;
  Derivatives[4] = Derivatives[34] = Where.x * Where.y;
  double y2 = Where.y * Where.y;
  Derivatives[5] = Derivatives[35] = y2;
  Derivatives[6] = Derivatives[36] = x2 * Where.x;
  Derivatives[7] = Derivatives[37] = x2 * Where.y;
  Derivatives[8] = Derivatives[38] = Where.x * y2;
  Derivatives[9] = Derivatives[39] = Where.y * y2;
  for (int i=10; i<30; ++i) Derivatives[i] = 0;
}

/* the ordering in the following array and the above routine are 
   linked together! */

GtransfoCub::CubParams GtransfoCub::Params[20] =
{{&GtransfoCub::dx},   {&GtransfoCub::a11},  {&GtransfoCub::a12},
 {&GtransfoCub::a1x2}, {&GtransfoCub::a1xy}, {&GtransfoCub::a1y2},
 {&GtransfoCub::a1x3 },{&GtransfoCub::a1x2y},{&GtransfoCub::a1xy2},
 {&GtransfoCub::a1y3 },

 {&GtransfoCub::dy},   {&GtransfoCub::a21},  {&GtransfoCub::a22},
 {&GtransfoCub::a2x2}, {&GtransfoCub::a2xy}, {&GtransfoCub::a2y2},
 {&GtransfoCub::a2x3 },{&GtransfoCub::a2x2y},{&GtransfoCub::a2xy2},
 {&GtransfoCub::a2y3 }
};


static double  ComputeChi2(const StarMatchList &List,
			   const Gtransfo *JustFitted)
{
  double chi2 = 0;
  
  for (StarMatchCIterator it = List.begin(); it != List.end(); ++it)
    {
    const StarMatch &a_match = *it;
    Point point1 = JustFitted->apply(a_match.point1);
    const Point &point2 = a_match.point2;
    chi2 += point1.Dist2(point2);
    }
  return chi2;
}
    
      
double  GtransfoCub::fit(const StarMatchList &List)
{
int npairs = List.size();
if (List.size()< 10)
  {
  cerr << " GtransfoCub::fit : trying to fit a cubic transfo with only " << npairs << " matches " << endl;
  return -1;
  }

#define NPAR 10
double a[NPAR][NPAR];  double bx[NPAR]; double by[NPAR]; double h[NPAR];

// fit T defined by
//  x' = dx + a11*x + a12*y + a1x2*x^2 + a1xy*x*y + a1y2*y^2 + a1x3*x^3 + a1x2y*x^2*y + a1xy2*x*y^2 + a1y3*y^3
//  y' = dy + a21*x + a22*y + a2x2*x^2 + a2xy*x*y + a2y2*y^2 + a2x3*x^3 + a2x2y*x^2*y + a2xy2*x*y^2 + a2y3*y^3

int i,j;
double sumr2 = 0; /* used to compute chi2 without relooping */

// loop on pairs and fill in the 2 Gram matrices ax and ay 
// h = (1 x y x^2 xy y^2 x^3 x^2y xy^2 y^3) 
// b is h * point2

StarMatchCIterator it;
memset(a,0,sizeof(double)*NPAR*NPAR);
memset(bx,0,sizeof(double)*NPAR);
memset(by,0,sizeof(double)*NPAR);

 // trick to overcome numerical problems (in case the averaage of x and y are way off 0)
 // we fit in an offset frame and transform back the result at the end.
 GtransfoLinShift shift_to_center= ShiftToCenter(List);

 for (it = List.begin(); it != List.end(); ++it)
    {
    const StarMatch &a_match = *it;
    Point point1 = shift_to_center.apply(a_match.point1);
    const Point &point2 = a_match.point2;
    h[0] = 1.0;
    h[1] = point1.x;
    h[2] = point1.y;
    h[3] = point1.x * point1.x;
    h[4] = point1.x * point1.y;
    h[5] = point1.y * point1.y;
    h[6] = point1.x * point1.x * point1.x;
    h[7] = point1.x * point1.x * point1.y;
    h[8] = point1.x * point1.y * point1.y;
    h[9] = point1.y * point1.y * point1.y;

    for (i=0; i< NPAR; i++)
      {
	bx[i] += h[i] * point2.x;
	by[i] += h[i] * point2.y;
	for (j=0; j<NPAR; j++) a[i][j] += h[i]*h[j];
      }
    sumr2 += point2.x*point2.x + point2.y*point2.y;
  }
// now invert the symetric matrices
 if (!SymMatInv(&a[0][0],NPAR))
   {
   cout << " could not invert the matrix  in GtransfoCub::fit " << endl;
   return -1;
   }
 
 double paramx[NPAR]; double paramy[NPAR];

 for (i=0; i<NPAR; i++) 
   {
     paramx[i] = 0; 
     paramy[i] = 0; 
     for (j=0;j<NPAR;j++) 
       {
	 paramx[i] += a[i][j]*bx[j];
	 paramy[i] += a[i][j]*by[j];
       }
   }

 
 // fill in the estimated parameters
 dx  = paramx[0];
 dy  = paramy[0];
 a11 = paramx[1];
 a21 = paramy[1];
 a12 = paramx[2];
 a22 = paramy[2]; 
 a1x2 = paramx[3];
 a2x2 = paramy[3];
 a1xy = paramx[4];
 a2xy = paramy[4];
 a1y2 = paramx[5];
 a2y2 = paramy[5];

 a1x3 = paramx[6];
 a2x3 = paramy[6];
 a1x2y = paramx[7];
 a2x2y = paramy[7];
 a1xy2 = paramx[8];
 a2xy2 = paramy[8];
 a1y3 = paramx[9];
 a2y3 = paramy[9];
 

 // go back to the original frame.
 *this = *this * shift_to_center;
 double chi2 = ComputeChi2(List, this );

 // compute and return chi2 for the two fits
 return chi2; // (sumr2 - ScalProd(paramx,bx,NPAR) - ScalProd(paramy,by,NPAR));

#undef NPAR
 
}

/*
typedef double GtransfoCub:: *GtransfoCubItem;
static struct { const char *name; GtransfoCubItem value; }  
*/

#ifndef __CINT__

GtransfoCub::CubAssoc GtransfoCub::WCS3Names [] = 
  {

    {"DX",    &GtransfoCub::dx   },
    {"DY",    &GtransfoCub::dy   },
    {"A11",   &GtransfoCub::a11  },
    {"A12",   &GtransfoCub::a12  },
    {"A21",   &GtransfoCub::a21  },
    {"A22",   &GtransfoCub::a22  },
    {"A1X2",  &GtransfoCub::a1x2 },
    {"A1XY",  &GtransfoCub::a1xy },
    {"A1Y2",  &GtransfoCub::a1y2 },
    {"A2X2",  &GtransfoCub::a2x2 },
    {"A2XY",  &GtransfoCub::a2xy },
    {"A2Y2",  &GtransfoCub::a2y2 },
    {"A1X3",  &GtransfoCub::a1x3 },
    {"A1X2Y", &GtransfoCub::a1x2y},
    {"A1XY2", &GtransfoCub::a1xy2},
    {"A1Y3",   &GtransfoCub::a1y3 },
    {"A2X3",  &GtransfoCub::a2x3 },
    {"A2X2Y", &GtransfoCub::a2x2y},
    {"A2XY2", &GtransfoCub::a2xy2},
    {"A2Y3",  &GtransfoCub::a2y3 },
    {"end",NULL}
};


   GtransfoCub::CubAssoc GtransfoCub::DJNames [] = {

    {"1_0",    &GtransfoCub::dx   },
    {"2_0",    &GtransfoCub::dy   },
    {"1_1",   &GtransfoCub::a11  },
    {"1_2",   &GtransfoCub::a12  },
    {"2_1",   &GtransfoCub::a21  },
    {"2_2",   &GtransfoCub::a22  },
    {"1_4",  &GtransfoCub::a1x2 },  /* 1_3 and 2_3 are coeffs for r */
    {"1_5",  &GtransfoCub::a1xy },
    {"1_6",  &GtransfoCub::a1y2 },
    {"2_4",  &GtransfoCub::a2x2 },
    {"2_5",  &GtransfoCub::a2xy },
    {"2_6",  &GtransfoCub::a2y2 },
    {"1_7",  &GtransfoCub::a1x3 },
    {"1_8", &GtransfoCub::a1x2y},
    {"1_9", &GtransfoCub::a1xy2},
    {"1_10",   &GtransfoCub::a1y3 },
    {"2_7",  &GtransfoCub::a2x3 },
    {"2_8", &GtransfoCub::a2x2y},
    {"2_9", &GtransfoCub::a2xy2},
    {"2_10",  &GtransfoCub::a2y3 },
    {"end",NULL}
};



/* 
   This mapping of the polynomial coefficients was proposed in a draft
   paper about distorions handling in WCS's (Calbretta and Greisen),
   but it disappeared in a later version.  It is no longer in the
   current version (by May 03) of the paper (which is really a draft).
   The actual mapping adopted is the one from Astrometrix (described
   in the documentation of the package), see
   http://terapix.iap.fr/soft/.

   The difference with the "DJ" stuff is that for the second polynomial
   (the DV2_# keys), the role of x and y are swapped with respect to the 
   first polynomial.

*/


   GtransfoCub::CubAssoc GtransfoCub::DVNames [] = {

    {"1_0",    &GtransfoCub::dx   },
    {"2_0",    &GtransfoCub::dy   },
    {"1_1",   &GtransfoCub::a11  },
    {"1_2",   &GtransfoCub::a12  },
    {"2_1",   &GtransfoCub::a22  },
    {"2_2",   &GtransfoCub::a21  },
    {"1_4",  &GtransfoCub::a1x2 },  /* 1_3 and 2_3 are coeffs for r */
    {"1_5",  &GtransfoCub::a1xy },
    {"1_6",  &GtransfoCub::a1y2 },
    {"2_4",  &GtransfoCub::a2y2 },
    {"2_5",  &GtransfoCub::a2xy },
    {"2_6",  &GtransfoCub::a2x2 },
    {"1_7",  &GtransfoCub::a1x3 },
    {"1_8", &GtransfoCub::a1x2y},
    {"1_9", &GtransfoCub::a1xy2},
    {"1_10",   &GtransfoCub::a1y3 },
    {"2_7",  &GtransfoCub::a2y3 },
    {"2_8", &GtransfoCub::a2xy2},
    {"2_9", &GtransfoCub::a2x2y},
    {"2_10",  &GtransfoCub::a2x3 },
    {"end",NULL}
};


void GtransfoCub::GetValues(vector<NamedValue> &Values, const string &KeyHeader) const
{
  Values.clear();
  GtransfoCub::CubAssoc *assoc = NULL;
  if (KeyHeader == "DJ") assoc = DJNames;
  else if (KeyHeader == "DV" || KeyHeader == "PV") assoc = DVNames;
  else if (KeyHeader == "WCS3") assoc = WCS3Names;
  else
    {
      cerr << "GtransfoCub::GetValues : unknown key Mapping scheme : "<< KeyHeader  << endl;
      return;
    }
  for (unsigned int i=0; assoc[i].value != NULL ; ++i)
    {
      Values.push_back(NamedValue(assoc[i].name, this->*(assoc[i].value)));
    }
}


void GtransfoCub::SetValues(vector<NamedValue> &Values, const string &KeyHeader)
{
  GtransfoCub::CubAssoc *assoc = NULL;
  if (KeyHeader == "DJ") assoc = DJNames;
  else if (KeyHeader == "DV" || KeyHeader == "PV") assoc = DVNames;
  else if (KeyHeader == "WCS3") assoc = WCS3Names;
  else
    {
      cerr << "GtransfoCub::SetValues : unknown key Mapping scheme : "<< KeyHeader  << endl;
      return;
    }
  for (unsigned int i=0; i < Values.size(); ++i)
    {
      string name = Values[i].name;
      unsigned int j;
      for (j=0; assoc[j].value != NULL ; ++j)
        {
	  if (assoc[j].name == name) break;
	}
      if (assoc[j].value == NULL)
	{
	  cerr << " GtransfoCub::SetValues : no item named " << name << endl;
	  continue;
	}
      this->*(assoc[j].value) = Values[i].value;
    }
}



GtransfoLin *GtransfoToLin(const Gtransfo* transfo)
{

  const GtransfoLin *lin = dynamic_cast<const GtransfoLin*> (transfo);
  if (lin != NULL) 
    {
      //cout << *lin << endl;
      return new GtransfoLin(*lin);
    }
  
  const GtransfoQuad *quad = dynamic_cast<const GtransfoQuad*> (transfo);
  if (quad != NULL) 
    {
      //cout << *quad << endl;
      return  new GtransfoLin(quad->dX(),  quad->dY(),
			      quad->A11(), quad->A12(), 
			      quad->A21(), quad->A22());
    }

  const GtransfoCub *cub = dynamic_cast<const GtransfoCub*> (transfo);
  if (cub != NULL) 
    {
      //cout << *cub << endl;
      return  new GtransfoLin(cub->dX(),  cub->dY(),
			      cub->A11(), cub->A12(), 
			      cub->A21(), cub->A22());
    }
  const GtransfoIdentity *id = dynamic_cast<const GtransfoIdentity*> (transfo);
  if (id != NULL) 
    {
      return new GtransfoLin();
    }
  cerr << " GtransfoToLin : Unable to transform transfo to linear one " << endl;
  return NULL;
}
#endif /* __CINT__ */


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
   radians to dgrees convertion. They are explicitely 
   done in apply (for TanPix2RaDec and TanRaDec2Pix).
   This is a minor concern though....
*/

TanPix2RaDec::TanPix2RaDec(const GtransfoLin &Pix2Tan, const Point &TangentPoint, 
	       const GtransfoQuad* Corrections)
{
  /* the angles returned by linPix2Tan should be in 
     degrees. */
  linPix2Tan = Pix2Tan;
  ra0  = deg2rad(TangentPoint.x);
  dec0 = deg2rad(TangentPoint.y);
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = NULL;
  if (Corrections) corr = new GtransfoCub(*Corrections);
}

// ": Gtransfo" suppresses a warning
TanPix2RaDec::TanPix2RaDec(const TanPix2RaDec &Original) : Gtransfo()
{
  corr = NULL;
  *this = Original;
}

void TanPix2RaDec::operator = (const TanPix2RaDec &Original)
{
  if (corr) delete corr;
  corr = NULL;
  linPix2Tan = Original.linPix2Tan;
  ra0 = Original.ra0;
  dec0 = Original.dec0;
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = NULL;
  if (Original.corr) corr = new GtransfoCub(*Original.corr);
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

void TanPix2RaDec::SetCorrections(const GtransfoQuad *Corrections)
{
  if (corr) delete corr;
  corr = NULL;
  if (Corrections)
    {
      corr = (GtransfoQuad*) Corrections->Clone();
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


const GtransfoQuad* TanPix2RaDec::Corr() const
{
  return corr;
}


Point TanPix2RaDec::CrPix() const
{
  /* CRPIX's are defined by:
                    ( CD1_1  CD1_2 )   (X - crpix1)
     transformed =  (              ) * (          )
                    ( CD2_1  CD2_2 )   (Y - crpix2)

     so that CrPix is the point which transforms to (0,0) 
  */
  GtransfoLin inverse = linPix2Tan.invert();
  return Point(inverse.dX(),inverse.dY());
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
  cerr << "TanPix2RaDec::fit is NOT implemented (although it is doable) " << endl;
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
  if (&Precision || &Region) {}
  return new TanPix2RaDec(LinPart().invert(),TangentPoint());
}


Gtransfo *TanRaDec2Pix::Clone() const
{
  return new TanRaDec2Pix(*this);
}

double TanRaDec2Pix::fit(const StarMatchList &List)
{
  if (&List) {} // warning killer
  cerr << " no way yet to fit TanRaDec2Pix (because it seemed useless!)" << endl;
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
  cerr << " no UserTransfo::fit  function defined " << endl;
  return -1;
}

Gtransfo *UserTransfo::Clone() const
{
  return new UserTransfo(*this);
}


/*************************************************************/

#ifdef USE_ROOT

ClassImp(Gtransfo);
ClassImp(GtransfoIdentity);
ClassImp(GtransfoLin);
ClassImp(GtransfoQuad);
ClassImp(GtransfoQuad);
ClassImp(GtransfoCub);
ClassImp(GtransfoComposition);

/*
RUN_ROOTCINT


LINKDEF_CONTENT : #pragma link C++ class Gtransfo+;
LINKDEF_CONTENT : #pragma link C++ class GtransfoLin+;
LINKDEF_CONTENT : #pragma link C++ class GtransfoIdentity+;
LINKDEF_CONTENT : #pragma link C++ class GtransfoQuad+;
LINKDEF_CONTENT : #pragma link C++ class GtransfoCub+;
LINKDEF_CONTENT : #pragma link C++ class GtransfoComposition+;
LINKDEF_CONTENT : #pragma link C++ function GtransfoCompose(const Gtransfo* Left, const Gtransfo* Right);
LINKDEF_CONTENT : #pragma link C++ function operator << (class ostream&, const Gtransfo&);
LINKDEF_CONTENT : #pragma link C++ function GtransfoRead(const char*, const char* ObjectName);
*/


#include "root_dict/gtransfodict.cc"
#endif



