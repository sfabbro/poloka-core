// This may look like C code, but it is really -*- C++ -*-
#ifndef GENERALFIT_SEEN
#define GENERALFIT_SEEN

#include <string>
#include "exceptions.h"

//================================================================
// GeneralFunction
//================================================================

class GeneralFunction EXC_AWARE {
public:
  GeneralFunction(unsigned int nVar, unsigned int nPar) {}
  virtual ~GeneralFunction() {}

  virtual double Value(double const xp[], double const* parm)=0;
  virtual double Val_Der(double const xp[], double const* parm
    , double* DgDpar) {}

  void SetDeltaParm(int numPar, double delta=0.);
  void SetDeltaParm(double const* dparam);

  inline int     NVar() const {return mNVar;}
  inline int     NPar() const {return mNPar;}

protected:
  const int mNVar;  // nombre de variables f(x,y,z,...)
  const int mNPar;  // nombre de parametres

  double *deltaParm;
  double *tmpParm;
};

//================================================================
// GeneralFunc
//================================================================

class GeneralFunc : public GeneralFunction {
public:
  GeneralFunc(unsigned int nvar, unsigned int npar
        ,double (*fun)(double const*,double const*) 
        ,double (*funder)(double const*, double const*, double*)=NULL);
  virtual ~GeneralFunc();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* parm
                        , double* DgDpar);

protected:
  double (*tmpFun)   (double const*, double const*);
  double (*tmpFunDer)(double const*, double const*, double*);
};

#ifdef IS_IT_USEFUL
//#include "matrix.h"
//#include "cvector.h"
//#include "generaldata.h"


//================================================================
// GeneralXi2
//================================================================

class GeneralFitData;

class GeneralXi2 EXC_AWARE {
public:
  GeneralXi2(unsigned int nPar);
  virtual ~GeneralXi2();

  virtual double Value(GeneralFitData& data, double* parm, int& ndataused)=0;
  virtual double Derivee(GeneralFitData& data, int i, double* parm);
  virtual double Derivee2(GeneralFitData& data, int i,int j, double* parm);

  void SetDeltaParm(int numPar, double delta=0.);
  void SetDeltaParm(double const* dparam);

  inline int     NPar() const {return mNPar;}

protected:
  const int mNPar;  // nombre de parametres

  double *deltaParm;
};

//================================================================
// GENERALFIT
//================================================================

class GeneralFit EXC_AWARE {
public:
  GeneralFit(GeneralFunction* f);
  GeneralFit(GeneralXi2* f);
  ~GeneralFit();

  void            WriteStep(char *filename = NULL);
  void            SetDebug(int level = 0);
  void            SetMaxStep(int n = 100);
  void            SetLambda_Fac(double fac = 10.);
  void            SetStopChi2(double s = 0.01);
  void            SetEps(double ep = 1.e-8);
  void            SetEps(int n,double ep = 1.e-8);
  void            SetStopMx(int nstopmx = 3, double stopchi2 = -1.);
  void            SetStopLent(int nstoplent = 3);
  void            SetFunction(GeneralFunction*);
  void            SetFuncXi2(GeneralXi2*);
  void            SetData(GeneralFitData*);
  void            SetParam(int n,double value, double step
                          ,double min=1., double max=-1.);
  void            SetParam(int n,string const&
                          ,double value, double step
                          ,double min=1., double max=-1.);
  void            SetParam(int n,double value);
  void            SetStep(int n,double value);
  void            SetMinStepDeriv(int i,double val = 0.);
  void            SetMinStepDeriv(double val = 0.);
  void            SetBound(int n,double min,double max);
  void            SetBound(int n);
  void            SetUnBound(int n);
  void            SetUnBound();
  void            SetFix(int n,double v);
  void            SetFix(int n);
  void            SetFree(int n);
  void            SetFree();

  double          GetParm(int);
  Vector          GetParm();
  double          GetParmErr(int);
  double          GetCoVar(int,int);
  double          GetStep(int n);
  double          GetMax(int n);
  double          GetMin(int n);
  inline double   GetChi2()    const {return Chi2;};
  inline double   GetChi2Red() const {
                        if(mNddl<=0) return (double) mNddl;
                        return Chi2/(double) mNddl;
                                     };
  inline double   GetEps(int i) const {return Eps(i);};
  inline int      GetNddl()    const {return mNddl;};
  inline int      GetNStep()   const {return nStep;};
  inline int      GetNVar()    const {return mNVar;};
  inline int      GetNPar()    const {return mNPar;};
  inline int      GetNFree()   const {return mNParFree;};
  inline int      GetNBound()  const {return mNParBound;};
  inline int      GetNStop()   const {return nStop;};
  inline int      GetNStopLent()   const {return nStopLent;};
  inline GeneralFunction*  GetFunction() const {return mFunction;};
  inline GeneralFitData*   GetGData() const {return mData;};

  void            PrintStatus();
  void            PrintFit();
  void            PrintParm(int n);
  void            PrintParm();

  int             Fit();
  double          ReCalChi2(int& nddl, double* par = NULL);
  GeneralFitData* DataResidus(bool clean=true);
  GeneralFitData* DataFunction(bool clean=true);
  void            PrintFitErr(int rc);

protected:
  int             mNtry;       // numero d'appel de la routine de fit.
  int             mNVar;       // nombre de variables f(x,y,z,...)
  int             mNPar;       // nombre de parametres
  int             mNParFree;   // nombre de parametres libres
  int             mNParBound;  // nombre de parametres bornes
  GeneralFunction*  mFunction;
  GeneralXi2*       mFuncXi2;
  GeneralFitData*   mData;

  Vector          Param;
  Vector          errParam;
  Vector          stepParam;
  Vector          minParam;
  Vector          maxParam;
  Vector          minStepDeriv;
  Vector          Eps;
  unsigned short int* fixParam;
  unsigned short int* boundParam;
  string*             nameParam;
  
  double          Lambda_Fac;
  double          stopChi2;
  int             maxStep;
  int             nStopMx;
  double          stopChi2SMx;
  int             nStopLent;
  int             debugLevel;
  FILE            *FileStep;
  
  Matrix          ATGA;
  Vector          BETA;
  Matrix          ATGA_Try;
  Vector          BETA_Try;
  Vector          C;
  Vector          D;
  
  double          Chi2;
  int             mNddl;
  int             nStep;
  int             nStop, nStopL;
  double          Lambda;
  
  // Fonctions privees
  void       write_in_step(double ci2,Vector& par);
  void       General_Init(void);
  void       TryFunc(Vector& par,Vector& par_tr);
  void       TryXi2(Vector& par,Vector& par_tr);
  void       CheckSanity();
  void       Set_Bound_C_D(int i);
  void       Set_Bound_C_D();
  double     p_vers_tr(int i,double p);
  Vector     p_vers_tr(Vector const& p);
  void       p_vers_tr(Vector const& p,Vector& tr);
  double     tr_vers_p(int i,double tr);
  Vector     tr_vers_p(Vector const& tr);
  void       tr_vers_p(Vector const& tr,Vector& p);
  double     c_dp_vers_dtr(int i,double tr);
  Vector       dp_vers_dtr(Vector const& dp,Vector const& tr);
  void         dp_vers_dtr(Vector const& dp,Vector const& tr,Vector& dtr);
  double     c_dtr_vers_dp(int i,double tr);
  Vector       dtr_vers_dp(Vector const& dtr,Vector const& tr);
  int        put_in_limits_for_deriv(Vector const& p,Vector& dp,double dist=0.66);
  inline void dtr_vers_dp(Vector const& dtr,Vector const& tr,Vector& dp)
              {  for(int i=0;i<mNPar;i++)
                   { if( fixParam[i] ) continue;
                     if( ! boundParam[i] ) dp(i) = dtr(i);
                       else dp(i) = D(i)/(1.+tr(i)*tr(i)) * dtr(i); }
              };
};

#endif
#endif /* IS_IT_USEFUL */
