#ifdef IS_IT_USEFUL
#include "defs.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

#if defined(__KCC__)
using std::string ;
#endif

#include "generalfit.h"
//#include "perrors.h"
//#include "cvector.h"
#include "nbconst.h"
//#include "tabmath.h"

#define EPS_FIT_MIN 1.e-8

//================================================================
// GeneralFunction
//================================================================

//++
// Class	GeneralFunction
// Lib	Outils++ 
// include	generalfit.h
//
//	Classe de fonctions parametrees a plusieurs variables.
//|           F[x1,x2,x3,...:a1,a2,a3,...]
//--

//////////////////////////////////////////////////////////////////////
//++
GeneralFunction::GeneralFunction(unsigned int nVar, unsigned int nPar)
//
//	Creation d'une fonction de `nVar' variables et `nPar' parametres.
//|  F[x(1),x(2),x(3),...x(nVar) : a(1),a(2),a(3),...,a(nPar)]
//--
  : mNVar(nVar), mNPar(nPar)
{
  // DBASSERT( nVar > 0 && nPar > 0 );
 deltaParm = new double[nPar];
 tmpParm   = new double[nPar];
 // END_CONSTRUCTOR
}

//++
GeneralFunction::~GeneralFunction()
//
//--
{
 delete[] deltaParm;
 delete[] tmpParm;
}

//////////////////////////////////////////////////////////////////////
//++
double GeneralFunction::Val_Der(double const xp[], double const* parm
                               , double *DgDpar)
//
//	Valeur et Derivees de la fonction (fct virtuelle par defaut).
//--
{
 for(int i=0;i<mNPar;i++) tmpParm[i] = parm[i];
 {for(int i=0;i<mNPar;i++) {
   double d = deltaParm[i];
   if(d==0.) { DgDpar[i] = 0.; continue;}
   tmpParm[i] -= d/2.;
   double vg = Value(xp,tmpParm);
   tmpParm[i] += d;
   double vd = Value(xp,tmpParm);
   DgDpar[i] = (vd - vg)/d;
   tmpParm[i] = parm[i];
 }}
 return Value(xp, parm);
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFunction::SetDeltaParm(int numPar, double d)
//
//	Definition de la variation du parametre numPar
//	pour calculer la derivee automatiquement.
//--
{
  // DBASSERT(numPar >= 0 && numPar < mNPar);
 deltaParm[numPar] = d;
}

//++
void GeneralFunction::SetDeltaParm(double const* dparam)
//
//	Idem precedente fonction mais pour tous les parametres
//--
{
 for(int i=0;i<mNPar;i++) deltaParm[i] = dparam[i];
}

//////////////////////////////////////////////////////////////////////
// Rappel des inline functions pour commentaires
//++
// virtual double Value(double const xp[], double const* parm)=0;
//	Valeur de la fonction a definir par l'utilisateur (fct virtuelle pure)
//--
//++
// inline int     NVar() const
//	Retourne le nombre de variables Xi
//--
//++
// inline int     NPar() const
//	Retourne le nombre de parametres Ai
//--

//================================================================
// GeneralFunc
//================================================================

//++
// Class	GeneralFunc
// Lib	Outils++ 
// include	generalfit.h
//
//	Classe de fonctions parametrees a plusieurs variables
//	derivant de ``GeneralFunction''. Permet de definir
//	une fonction a fiter sans passer par une classe derivee
//	en utilisant l'ecriture courante du C. La fonction
//	retournant les derivees par rapport aux parametres du fit
//	peut etre egalement fournie (optionnel).
//--

/////////////////////////////////////////////////////////////////
//++
GeneralFunc::GeneralFunc(unsigned int nvar, unsigned int npar, double (*fun) (double const*, double const*)
                        , double (*funder) (double const*, double const*, double*) )
//
//	Createur, on passe le nom ``fun'' de la fonction a la mode C.
//	On peut optionellement egalement passer le nom de la fonction
//	``funder'' qui retourne les valeurs des derivees par rapport
//	aux parametres du fit.
//--
//++
//| ----------------------
//| Exemple d'utilisation:
//| ----------------------
//| include "generalfit.h"
//| ...
//| double   gaussc(double const* x,double const* p);
//| double d_gaussc(double const* x,double const* p,double* dp);
//| ...
//| main {
//|  ...
//|  // Fit SANS calcul automatique des derivees
//|  GeneralFunc      myfunc(2,7,gaussc);
//|  GeneralFit       myfit(&myfunc);
//|  ...
//|  myfit.Fit();
//|  ...
//|  // Fit AVEC calcul automatique des derivees
//|  GeneralFunc      myfunc(2,7,gaussc,d_gaussc);
//|  GeneralFit       myfit(&myfunc);
//|  ...
//|  myfit.Fit();
//| }
//--
//++
//| // Definition de la fonction a fitter a la mode C
//| double gaussc(double const* x,double const* p)
//| // Fonction: X=(x[0]-p[1])/p[3], Y=(x[1]-p[2])/p[4],
//| //  f = p[0]*exp{-0.5*[X^2+Y^2-2*p[5]*X*Y]} + p[6]
//| {
//|  double X = (x[0]-p[1])/p[3];
//|  double Y = (x[1]-p[2])/p[4];
//|  return p[0]*exp(-(X*X+Y*Y-2*p[5]*X*Y)/2)+p[6];
//| }
//| // Definition de la fonction des derivees / parametres
//| // Cette fonction retourne aussi la valeur de la fonction a fitter.
//| double d_gaussc(double const* x,double const* p,double* dp)
//| {
//|  dp[0] = derivee de gaussc par rapport au parametre p[0]
//|  ...
//|  dp[6] = derivee de gaussc par rapport au parametre p[6]
//|  return gaussc(x,p);
//| }
//--
: GeneralFunction(nvar,npar), tmpFun(fun), tmpFunDer(funder)
{
}

GeneralFunc::~GeneralFunc()
{
}

double GeneralFunc::Value(double const xp[], double const* Par)
{
return tmpFun(xp,Par);
}

double GeneralFunc::Val_Der(double const xp[],double const* parm, double* DgDpar)
{
if(tmpFunDer) return tmpFunDer(xp,parm,DgDpar);
  else        return GeneralFunction::Val_Der(xp,parm,DgDpar);
}

//================================================================
// GeneralXi2
//================================================================

//++
// Class	GeneralXi2
// Lib	Outils++ 
// include	generalfit.h
//
//	Classe de Xi2 a plusieurs parametres.
//|           Xi2[a1,a2,a3,...]
//--

//////////////////////////////////////////////////////////////////////
//++
GeneralXi2::GeneralXi2(unsigned int nPar)
//
//	Creation d'un Xi2 de `nPar' parametres.
//|  Xi2[a(1),a(2),a(3),...,a(nPar)]
//--
  : mNPar(nPar)
{ 
  // DBASSERT( nPar>0 );
 deltaParm = new double[nPar];
 // END_CONSTRUCTOR
}

//++
GeneralXi2::~GeneralXi2()
//
//--
{
 delete[] deltaParm;
}

//////////////////////////////////////////////////////////////////////
//++
double GeneralXi2::Derivee(GeneralFitData& data, int i, double* parm)
//
//	Derivee du Xi2 par rapport au parametre `i'
//	pour les valeurs `parm' des parametres.
//--
{
 int dum;
 double d = deltaParm[i];
 parm[i] -= d/2.;
 double vg = Value(data, parm,dum);
 parm[i] += d;
 double vd = Value(data, parm,dum);
 parm[i] -= d/2.;
 return (vd - vg)/d;
}

//++
double GeneralXi2::Derivee2(GeneralFitData& data, int i, int j, double* parm)
//
//	Derivee seconde du Xi2 par rapport aux parametres `i' et `j'
//	pour les valeurs `parm' des parametres. Attention, cette fonction
//	calcule d/di(dC2/dj), valeur qui est numeriquement differente
//	de d/dj(dC2/di).
//--
//++
//|
//| **** Remarque: Derivee2 = dXi2/dPi.dPj represente le Hessien.
//| Derivee2(k,l)= dXi2/dPk.dPl
//|              = 2*SUMi{1/Si^2*[df(xi;P)/dPk * df(xi;P)/dPl]
//|                       + [yi-f(xi;P)] * df(xi;P)/dPk.dPl }
//| ou (xi,yi) sont les points de mesure. "Si" l'erreur sur le point i
//|    SUMi represente la somme sur les points de mesure
//|    f(x;P) represente le modele parametrique a fitter
//|    "P" represente l'ensemble des parametres et "Pi" le ieme parametre
//| Les composantes du Hessien dependent des derivees 1ere et 2sd du modele
//| a fitter f(x;P) selon les parametres "Pi". La prise en compte des derivees
//| secondes est un facteur destabilisant. De plus le facteur [yi-f(xi;P)]
//| devant la derivee 2sd est seulement l'erreur de mesure aleatoire qui
//| n'est pas correlee avec le modele. Le terme avec la derivee 2sd
//| tend donc a s'annuler et peut donc etre omis.
//| (cf. Numerical Recipes in C, chap 15 Modeling of Data, Nonlinear Models,
//|  Calculation of the Gradient and Hessian p682,683)
//|
//| **** Conseil: Il est conseille a l'utilisateur de sur-ecrire
//| la fonction virtuelle Derivee2 et de la remplacer par:
//| Derivee2(k,l) = 2*SUMi{1/Si^2*[df(xi;P)/dPk * df(xi;P)/dPl]}
//--
{
 double d = deltaParm[i];
 parm[i] -= d/2.;
 double vg = Derivee(data,j,parm);
 parm[i] += d;
 double vd = Derivee(data,j,parm);
 parm[i] -= d/2.;
 d = (vd - vg)/d;
 return d;
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralXi2::SetDeltaParm(int numPar, double d)
//
//	Definition de la variation du parametre numPar
//	pour calculer la derivee automatiquement.
//--
{
  // DBASSERT(numPar >= 0 && numPar < mNPar);
  
 deltaParm[numPar] = d;
}

//++
void GeneralXi2::SetDeltaParm(double const* dparam)
//
//	Idem precedente fonction mais pour tous les parametres.
//--
{
 for(int i=0;i<mNPar;i++) deltaParm[i] = dparam[i];
}

//////////////////////////////////////////////////////////////////////
// Rappel des inline functions pour commentaires
//++
// virtual double Value(GeneralFitData& data, double const* parm, int& ndataused)=0;
//	Valeur du Xi2 a definir par l'utilisateur (fct virtuelle pure)
//	a partir des donnees de `data'. l'utilisateur doit egalement
//	retourner le nombre de points de mesure utilises dans le calcul
//	du Xi2 (`ndataused').
//--
//++
// inline int     NPar() const
//	Retourne le nombre de parametres Ai.
//--

//================================================================
// GeneralFit
//================================================================
//                                Christophe 8/11/93 La Silla
//                                re-codage C++ 16/01/96 Saclay

//++
// Class	GeneralFit
// Lib		Outils++
// include	generalfit.h
//
//	Classe de fit d'une GeneralFunction sur une GeneralFitData
//--

//////////////////////////////////////////////////////////////////////
//++
GeneralFit::GeneralFit(GeneralFunction* f)
//
//	Creation d'une classe de fit pour la `GeneralFunction f'.
//--
  : mNVar         (f->NVar()),
    mNPar         (f->NPar()),
    mFunction     (f),
    mFuncXi2      (NULL),

    Param         (f->NPar()),
    errParam      (f->NPar()),
    stepParam     (f->NPar()),
    minParam      (f->NPar()),
    maxParam      (f->NPar()),
    minStepDeriv  (f->NPar()),
    Eps           (f->NPar()),
  
    ATGA          (f->NPar(), f->NPar()),
    BETA          (f->NPar()),
    ATGA_Try      (f->NPar(), f->NPar()),
    BETA_Try      (f->NPar()),
    C             (f->NPar()),
    D             (f->NPar())
{
  // DBASSERT(mNVar>0 && mNPar>0);
  // DBASSERT(mNPar<1000000);

 TRY {
   General_Init();
 } CATCHALL {
   THROW_SAME;
 } ENDTRY

 END_CONSTRUCTOR
}

//++
GeneralFit::GeneralFit(GeneralXi2* f)
//
//	Creation d'une classe de fit pour le `GeneralXi2 f'.
//	L'emploi de cette methode n'est pas conseillee car elle
//	calcule automatiquement la derivee 2sd du Xi2 par rapport
//	aux parametres, ce qui entraine un manque de robustesse
//	et qui ne garanti pas que la matrice de covariance soit
//	definie positive (il est possible de surecrire
//	la methode virtuelle Derivee2 pour palier ce probleme).
//--
  : mNVar         (0),
    mNPar         (f->NPar()),
    mFunction     (NULL),
    mFuncXi2      (f),

    Param         (f->NPar()),
    errParam      (f->NPar()),
    stepParam     (f->NPar()),
    minParam      (f->NPar()),
    maxParam      (f->NPar()),
    minStepDeriv  (f->NPar()),
    Eps           (f->NPar()),
  
    ATGA          (f->NPar(), f->NPar()),
    BETA          (f->NPar()),
    ATGA_Try      (f->NPar(), f->NPar()),
    BETA_Try      (f->NPar()),
    C             (f->NPar()),
    D             (f->NPar())
{
  // DBASSERT( mNPar>0 );
  // DBASSERT( mNPar < 1000000 );

 TRY {
   General_Init();
 } CATCHALL {
   THROW_SAME;
 } ENDTRY

 END_CONSTRUCTOR
}

//
void GeneralFit::General_Init(void)
// Initialisation des diverses variables
{
 mNtry      = 0;
 mNParFree  = mNPar;
 mNParBound = 0;

 mData      = NULL;

 fixParam   = NULL;
 boundParam = NULL;
 nameParam  = NULL;

 Lambda_Fac = 10.;
 stopChi2   = 0.01;
 maxStep    = 100;
 nStopMx    = 3;
 stopChi2SMx = stopChi2;
 nStopLent  = 0;
 debugLevel = 0;
 FileStep = NULL;
  
 Chi2       = 0.;
 mNddl      = -1;
 nStep      = 0;
 nStop      = 0;
 nStopL     = 0;
 Lambda     = 0.001;

 GetIntEnv("PDEBUG_GENERALFIT",debugLevel);

 TRY {
   fixParam   = new unsigned short int[mNPar];
   boundParam = new unsigned short int[mNPar];
   nameParam  = new string[mNPar];
 } CATCHALL {
   cout<<"GeneralFit::GeneralFit Impossible d'allouer l'espace"<<endl;
   THROW_SAME;
 } ENDTRY

 Param        = (double)  0.;
 errParam     = (double)  0.;
 stepParam    = (double)  1.;
 minParam     = (double)  1.;
 maxParam     = (double) -1.;
 minStepDeriv = (double) 0.;
 Eps          = (double) EPS_FIT_MIN;
 char str[8];
 for(int i=0;i<mNPar;i++) {
   sprintf(str,"P%d",i);
   fixParam[i]   = 0;
   boundParam[i] = 0;
   nameParam[i]  = str;
 }
}

//++
GeneralFit::~GeneralFit()
//
//--
{
 delete[] fixParam;
 delete[] boundParam;
 delete[] nameParam;
 if(FileStep!=NULL) fclose(FileStep);
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::WriteStep(char *filename)
//
//	Pour ecrire les iterations dans le fichier filename
//--
{

//#if defined(__DECCXX) || defined(__KCC__) || defined(__aCC__) 
if(filename==NULL) filename = const_cast<char *>("generalfit.iter");
//#else
//if(filename==NULL) filename = "generalfit.iter";
//#endif
FileStep = fopen(filename,"w");
if(FileStep==NULL) THROW(nullPtrErr);
}

//++
void GeneralFit::SetDebug(int level)
//
//	Niveau de debug
//	(voir aussi la variable d'environnement PDEBUG_GENERALFIT).
//--
{
 debugLevel = ( level < 0 ) ? 0: level;
 if(debugLevel>0) cout<<"SetDebug_level "<<debugLevel<<endl;
}

//++
void GeneralFit::SetMaxStep(int n)
//
//	Nombre maximum d'iterations permis.
//--
{
 maxStep = ( n <= 1 ) ? 100: n;
 if(debugLevel>0) cout<<"SetMaxStep "<<maxStep<<endl;
}

//++
void GeneralFit::SetLambda_Fac(double fac)
//
//	Facteur de multiplication/division de Lambda selon
//	que le Chi2 a augmente ou diminue.
//--
{
  Lambda_Fac = (fac>1.) ? fac : 10.;
}

//++
void GeneralFit::SetStopChi2(double s)
//
//	Critere de convergence sur le Chi2.
//--
{
 stopChi2 = ( s <= 0. ) ? 0.01: s;
 if(debugLevel>0) cout<<"SetStopChi2 "<<stopChi2<<endl;
}

//++
void GeneralFit::SetEps(double ep)
//
//	Precision des calculs (cf descriptif general).
//--
{
 ep = (ep<=0.) ? EPS_FIT_MIN: ep;
 if(debugLevel>0) cout<<"SetEps "<<ep<<endl;
 for(int i=0;i<mNPar;i++) SetEps(i,ep);
}

//++
void GeneralFit::SetEps(int n,double ep)
//
//	Precision des calculs pour le parametre n.
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 Eps(n) = (ep<=0.) ? EPS_FIT_MIN: ep;
 if(debugLevel>0) cout<<"SetEps("<<n<<") = "<<Eps(n)<<endl;
}

//++
void GeneralFit::SetStopMx(int nstopmx,double stopchi2)
//
//	Critere de convergence sur le nombre de stop en chi2
//	dans le cas ou le chi2 augmente de moins de stopchi2
//	(cf descriptif general).
//	Si nstopmx<=0, alors ce critere de convergence n'est pas applique.
//	Si stopchi2<=0, alors la valeur generale mise par SetStopChi2()
//	est utilisee.
//--
{
 nStopMx = (nstopmx>0) ? nstopmx : 0;
 stopChi2SMx = (stopchi2>0.) ? stopchi2 : stopChi2;
 if(debugLevel>0) cout<<"SetStopMx: nStopMx="<<nStopMx
                      <<" stopChi2SMx="<<stopChi2SMx<<endl;
}

//++
void GeneralFit::SetStopLent(int nstoplent)
//
//	Critere de convergence sur le nombre de stop en chi2
//	dans le cas ou le chi2 diminue (cf descriptif general).
//	Si nstopl<=0, alors ce critere de convergence n'est pas applique.
//--
{
 nStopLent = (nstoplent>0) ? nstoplent : 0;
 if(debugLevel>0) cout<<"SetStopLent "<<nStopLent<<endl;
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetFunction(GeneralFunction* f)
//
//	Pour changer la fonction a fitter en cours de route
//	(On ne peut passer d'un fit sur une GeneralFunction
//	a un fit sur un GeneralXi2 sans recreer la classe).
//--
{
  // DBASSERT( mFuncXi2  == NULL );
  // DBASSERT( f != NULL );
  // DBASSERT( f->NVar() == mNVar );
  // DBASSERT( f->NPar() == mNPar );
 mFunction = f;
 if(debugLevel>0) cout<<"SetFunction "<<mFunction<<endl;
}

//++
void GeneralFit::SetFuncXi2(GeneralXi2* f)
//
//	Pour changer le Xi2 a fitter en cours de route
//	(On ne peut passer d'un fit sur une GeneralFunction
//	a un fit sur un GeneralXi2 sans recreer la classe).
//--
{
  // DBASSERT( mFunction == NULL );
  // DBASSERT( f != NULL );
  // DBASSERT( f->NPar() == mNPar );
 mFuncXi2  = f;
 if(debugLevel>0) cout<<"SetFuncXi2 "<<mFuncXi2<<endl;
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetData(GeneralFitData* data)
//
//	Pour connecter une structure de donnees.
//--
{
  // DBASSERT( data->NVar()==mNVar || mFunction==NULL );
 mData = data;
 mNddl = mData->NDataGood() - mNParFree;
 if(debugLevel>0)
   cout<<"SetData "<<mData<<" data pour "<<mNddl<<" ddl"<<endl;
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetParam(int n,double value,double step
                         ,double min,double max)
//
//	Definition du parametre "n" a fitter.
//--
{
  // DBASSERT(n>=0 && n<mNPar);

 Param(n)     = value;
 if(step>0.) {
   if( fixParam[n] ) { fixParam[n]=0; mNParFree++;}
 } else {
   if( ! fixParam[n] ) { fixParam[n]=1; mNParFree--;}
 }
 stepParam(n) = step;
 minParam(n)  = min;
 maxParam(n)  = max;
 if(max>min) {
   if( ! boundParam[n] ) {boundParam[n]=1; mNParBound++;}
 } else {
   if( boundParam[n] ) {boundParam[n]=0; mNParBound--;}
 }
  
 if(debugLevel) {cout<<"Set_"; PrintParm(n);}
}

//++
void GeneralFit::SetParam(int n, string const& name
                        ,double value,double step,double min,double max)
//
//	Definition du parametre "n" a fitter
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 SetParam(n,value,step,min,max);
 nameParam[n] = name;
 if(debugLevel) {cout<<"Set_Param "; PrintParm(n);}
}

//++
void GeneralFit::SetParam(int n,double value)
//
//	Definition du parametre "n" a fitter
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 Param(n) = value;
 if(debugLevel) {cout<<"Set_Param "; PrintParm(n);}
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetStep(int n,double step)
//
//	Definition du pas de depart du parametre "n"
//	Si negatif ou nul, parametre fixe.
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 if(step>0.) {
   if( fixParam[n] ) { fixParam[n]=0; mNParFree++;}
 } else {
   if( ! fixParam[n] ) { fixParam[n]=1; mNParFree--;}
 }
 stepParam(n) = step;
 if(debugLevel) {cout<<"Set_Step"; PrintParm(n);}
}

//++
void GeneralFit::SetMinStepDeriv(int i,double val)
//
//	Definition du pas minimum `val' pour le parametre `i'
//	pouvant etre utilise dans le calcul automatique des derivees
//	(soit de la fonction, soit du Xi2 selon les parametres du fit).
//	Si nul pas de limite, si negatif alors `EPS(i)' (cf SetEps).
//	Inutile dans le cas ou les derivees sont donnees
//	par l'utilisateur.
//--
{
  // DBASSERT(i>=0 && i<mNPar);
 if(val<0.) minStepDeriv(i) = Eps(i);
   else     minStepDeriv(i) = val;
 if(debugLevel>0) cout<<"SetMinStepDeriv("<<i<<") = "<<minStepDeriv(i)<<endl;
}

//++
void GeneralFit::SetMinStepDeriv(double val)
//
//	Definition du pas minimum `val' pour tout les parametres
//	(voir description SetMinStepDeriv ci-dessus).
//--
{
 if(debugLevel>0) cout<<"SetMinStepDeriv "<<val<<endl;
 for(int i=0;i<mNPar;i++) SetMinStepDeriv(i,val);
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetBound(int n, double min, double max)
//
//	Definition des bornes du parametre "n"
//	Si max<=min, parametre non-borne.
//--
{
  // DBASSERT(n>=0 && n<mNPar && max>min);

 minParam(n)  = min;
 maxParam(n)  = max;
 if( ! boundParam[n] ) {
   boundParam[n] = 1;
   mNParBound++;
   if(debugLevel>0)
     cout<<"SetBound "<<n<<" min="<<min<<" max="<<max
         <<" (Nbound="<<mNParBound<<")"<<endl;
 }
}

//++
void GeneralFit::SetBound(int n)
//
//	Pour re-borner le parametre "n" aux bornes par defaut
//--
{
  // DBASSERT(n>=0 && n<mNPar && maxParam(n)>minParam(n));
 SetBound(n,minParam(n),maxParam(n));
}

//++
void GeneralFit::SetUnBound(int n)
//
//	Pour ne plus borner le parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);

 if( boundParam[n] ) {
   boundParam[n] = 0;
   mNParBound++;
   if(debugLevel>0) cout<<" SetUnBound "<<n
                        <<" (Nbound="<<mNParBound<<")"<<endl;
 }
}

//++
void GeneralFit::SetUnBound()
//
//	Pour ne plus borner tous les parametres
//--
{
 for(int i=0;i<mNPar;i++) SetUnBound(i);
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::SetFix(int n,double v)
//
//	Pour fixer le parametre "n" a la valeur "v"
//--
{
  // DBASSERT(n>=0 && n<mNPar);

 Param(n) = v;
 if( ! fixParam[n] ) {
   fixParam[n] = 1;
   mNParFree--;
 }
 if(debugLevel>0) cout<<" SetFix "<<n
                      <<" v="<<v
                      <<" (Nfree="<<mNParFree
                      <<")"<<endl;
}

//++
void GeneralFit::SetFix(int n)
//
//	Pour fixer le parametre "n" a la valeur par defaut
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 SetFix(n,Param(n));
}

//++
void GeneralFit::SetFree(int n)
//
//	Pour liberer le parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);

 if( fixParam[n] ) {
   fixParam[n] = 0;
   mNParFree++;
   if(debugLevel>0) cout<<" SetFree "<<n
                        <<"   Step "<<stepParam(n)
                        <<" (Nfree="<<mNParFree<<")"<<endl;
   if(stepParam(n)<=0.)
     cout<<"ATTENTION SetFree["<<n<<"] avec step<=0 "
         <<stepParam(n)<<endl;
 }
}

//++
void GeneralFit::SetFree()
//
//	Pour liberer tous les parametres
//--
{
 for(int i=0;i<mNPar;i++) SetFree(i);
}

//////////////////////////////////////////////////////////////////////
//++
double GeneralFit::GetParm(int n)
//
//	Retourne la valeur du parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 return Param(n);
}

//++
Vector GeneralFit::GetParm()
//
//	Retourne les valeurs des parametres dans un vecteur.
//--
{
return Param;
}

//++
double GeneralFit::GetParmErr(int n)
//
//	Retourne la valeur de l'erreur du parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 return errParam(n);
}

//++
double GeneralFit::GetCoVar(int i,int j)
//
//	Retourne la covariance pour les parametre `i' et `j'
//--
{
  // DBASSERT(i>=0 && i<mNPar && j>=0 && j<mNPar);
 return ATGA(i,j);
}

//////////////////////////////////////////////////////////////////////
//++
double GeneralFit::GetStep(int n)
//
//	Retourne la valeur du pas du parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 return stepParam(n);
}

//++
double GeneralFit::GetMax(int n)
//
//	Retourne la valeur de la borne superieure du parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 return maxParam(n);
}

//++
double GeneralFit::GetMin(int n)
//
//	Retourne la valeur de la borne inferieure du parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);
 return minParam(n);
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::PrintStatus()
//
//	Impression du status du fit
//--
{
 cout<<"GeneralFit::PrintStatus"
     <<" mData="<<mData
     <<" mFunction="<<mFunction
     <<" mFuncXi2="<<mFuncXi2
     <<endl;
 cout<<" mNVar="<<mNVar
     <<" mNPar="<<mNPar
     <<" mNParFree="<<mNParFree
     <<" mNParBound="<<mNParBound
     <<endl;
 cout<<" Lambda_Fac="<<Lambda_Fac
     <<" stopChi2="<<stopChi2
     <<" maxStep="<<maxStep
     <<" nStopMx="<<nStopMx<<" stopChi2SMx="<<stopChi2SMx
     <<" nStopLent="<<nStopLent
     <<" debugLevel="<<debugLevel
     <<endl;
 PrintParm();
}

//++
void GeneralFit::PrintFit()
//
//	Impression des resultats du fit
//--
{
 cout<<"PrintFit: Chi2="<<Chi2
     <<" Lambda="<<Lambda
     <<" nStep="<<nStep
     <<" nStop="<<nStop
     <<" nStopL="<<nStopL
     <<" nDDL="<<mNddl
     <<endl;
 PrintParm();
}

//++
void GeneralFit::PrintParm(int n)
//
//	Impression des informations relatives au parametre "n"
//--
{
  // DBASSERT(n>=0 && n<mNPar);

 cout<<"Par["<<n<<"] "<<nameParam[n]
     <<" F"<<fixParam[n]
     <<" B"<<boundParam[n]
     <<" : "<<Param(n)
     <<" +/- "<<errParam(n)
     <<" : "<<stepParam(n)
     <<" "<<minParam(n)
     <<" "<<maxParam(n)
     <<" : "<<Eps(n)
     <<" "<<minStepDeriv(n)
     <<endl;
}

//++
void GeneralFit::PrintParm()
//
//	Impression des informations relatives a tous les parametres
//--
{
 cout<<"*** Parametres : fix bnd : par err : step min max : eps dmin\n";
 for (int i=0; i<mNPar; i++) PrintParm(i);
 cout<<endl;
}

//////////////////////////////////////////////////////////////////////
//++
int GeneralFit::Fit()
//
//--
//++
//|  Fonction de fit de la fonction f(x,y,z,...:p1,p2,...,pn)
//|          sur les donnees x[i],y[i],z[i],...,F[i],ErrF[i]
//|  - Methode:   fit des moindres carres dans le cas non lineaire
//|  - Reference: Statistical and Computational Methods in Data Analysis
//|              Siegmund Brandt, North-Holland 1970  p 204-206.
//|              Introduction des limites pour la variation des parametres (cmv).
//|              Increment des parametres selon la methode de Levenberg-Marquardt
//|  (Numerical Recipes in C, chap 15 Modeling of Data, Nonlinear Models,
//|   Levenberg-Marquardt Method p683)
//--
//++
//|  - Gestion des parametres bornes:
//|     si p est un parametre borne entre pmin et pmax, le parametre fitte est q
//|     tel que     q = tang((p-C)/D)    ....   p = C + D*atan(q)
//|     ou   C = (pmin+pmax)/2.  et   D = (pmax-pmin)/Pi
//|     On a  dq = (1+q**2)/D * dp    ....   dp = D/(1+q**2) * dq
//|     et    dF/dq = dF/dp * dp/dq = D/(1+q**2) * dF/dp
//|           dF/dp = dF/dq * dq/dp = (1+q**2)/D * dF/dp
//|                                   ^ q
//|                   |               |              *| "tang()"
//|                   |               |              *|
//|                   |               |              *|
//|                   |               |             * |
//|                   |               |            *  |
//|                   |               |          *    |
//|                   |               |       *       |
//|               Pmin|              C|   *           |Pmax
//|     --------------|---------------*---------------|--------------> p
//|              -Pi/2|           *   |0              |Pi/2
//|                   |       *       |               |
//|                   |    *          |               |
//|                   |  *            |               |
//|                   | *             |               |
//|                   |*              |               |
//|                   |*              |               |
//|                   |*              |               |
//|                   <------------------- D --------->
//--
//++
//|  - Criteres de convergence, arrets standards:
//|    - SOIT: le Chi2 est descendu de moins de stopChi2
//|            entre l'iteration n et n+1
//|            (stopChi2 est change par SetStopChi2)
//|    - SOIT: 1. le chi2 est remonte de moins de stopChi2SMx et
//|            2. les parametres libres ont varie de moins de Eps(i)
//|               pendant les nStopmx dernieres iterations
//|            Si nStopmx<=0, alors ce critere n'est pas applique (def=3).
//|            (nStopmx,stopChi2SMx sont changes par SetStopMx, Eps par SetEps)
//|
//|  - Criteres de convergence, arrets par non-convergence:
//|    - plus de "maxStep" iterations.
//|
//|  - Criteres de convergence, arrets speciaux:
//|    - Si l'utilisateur a demande explicitement la methode d'arret
//|      "SetStopLent()", arret si :
//|        1. le Chi2 est descendu et
//|        2. les parametres libres ont varies de moins de Eps
//|           pendant les nStopLent dernieres iterations.
//|           (nStopLent est change par SetStopLent, Eps par SetEps)
//|
//--
//++
//|  - Remarques diverses:
//|     Les points avec erreurs <=0 ne sont pas utilises dans le fit.
//|     Les bornes des parametres ne peuvent etre atteintes
//|  - entrees:
//|     la fonction est definie par une classe GeneralFunction
//|     les donnees sont passees par une classe GeneralFitData
//|     le nombre de parametres et le nombre de variables doivent etre
//|        coherents entre GeneralFunction GeneralFitData GeneralFit
//|  - Return: 
//|       la function elle meme retourne le nombre d'iterations  du fit si succes
//|       -1  : si le nombre de degre de liberte est <0
//|       -10 : si l'inversion de la matrice des erreurs n'a pas ete possible
//|       -11 : si un element diagonal de la matrice des covariances est <=0
//|       -20 : si le fit n'a pas converge (nstep>nNstepMX)
//|       -100-N : si le parametre "N" est initialise hors limites
//|       -200-N : si le parametre "N" atteint sa limite inferieure
//|       -300-N : si le parametre "N" atteint sa limite superieure
//--
{
 volatile double oldChi2;
 Matrix COVAR(mNPar,mNPar);
 Vector DA(mNPar);
 Vector dparam(mNPar);
 Vector paramTry(mNPar);
 Vector param_tr(mNPar);
 Vector paramTry_tr(mNPar);
 Vector step_tr(mNPar);
 nStop = nStopL = nStep = 0;
 Chi2 = oldChi2 = 0.;
 Lambda = 0.001;
 mNddl = mData->NDataGood() - mNParFree;
 if(mNddl<0) return -1;
 mNtry++;

 if(debugLevel>= 2)
   cout<<"\n********* DEBUT GENERALFIT.FIT() **************"<<endl;

 // set matrices C,D dans le cas de parametres bornes
 if(mNParBound>0) Set_Bound_C_D();

 if(debugLevel>= 2) PrintStatus();

 // check de la coherence des operations et assignations
 CheckSanity();

 // Pour les parametres bornes on verifie
 // qu'ils sont initialises dans leurs limites
 {for(int i=0;i<mNPar;i++) 
   {
     if( !boundParam[i] || fixParam[i] ) continue;
     if( minParam(i)<Param(i) && Param(i)<maxParam(i) ) continue;
     /* if(debugLevel>= 1) */
     cout<<"Parametre "<<i<<" initialise hors limites "
	 <<minParam(i)<<" < "<<Param(i)
	 <<" < "<<maxParam(i)<<endl;
     return(-100-i);
 }}

 // premier essai d'initialisation
 param_tr = p_vers_tr(Param);
 dparam = stepParam / 2.;
 put_in_limits_for_deriv(Param,dparam);
 if(mFunction!=NULL) mFunction->SetDeltaParm(dparam.Data());
   else if(mFuncXi2!=NULL) mFuncXi2->SetDeltaParm(dparam.Data());
 step_tr = dp_vers_dtr(stepParam,param_tr);

 if(debugLevel>= 2) {
   cout<<"ESSAI numero 1: Param:"<<endl;
   cout<<Param;
   cout<<"param_tr:"<<endl;
   cout<<param_tr;
   cout<<"step_tr:"<<endl;
   cout<<step_tr;
 }
 if(mFunction!=NULL) TryFunc(Param,param_tr);
   else if(mFuncXi2!=NULL) TryXi2(Param,param_tr);
 ATGA = ATGA_Try;
 BETA = BETA_Try;
 oldChi2 = Chi2;

 // Iterations
 while (1) {
   nStep++;

   // un nouvel essai (si Lambda!=0)
   {for(int i=0; i<mNPar; i++)
       if(! fixParam[i] ) ATGA(i,i) *= 1 + Lambda;
          else  ATGA(i,i) = 1.;}

   // Calcul de la matrice des covariances
#ifdef __mac__
   COVAR = ATGA.Inverse();
#else
   TRY {
     COVAR = ATGA.Inverse();
   } CATCHALL {
     if(debugLevel>0) {
       cout<<"Pb inversion matrice ATGA:"<<endl;
       cout<<ATGA;
     }
     return(-10);
   } ENDTRY
#endif

   if (debugLevel >= 3) {
     cout<<"Matrice (tA G A)^-1 = \n";
     cout<<COVAR;
   }

   // calculs des deplacements a effectuer
   DA = COVAR * BETA;
   if (debugLevel >=2) {
     cout<<"Correction parametres DA : \n";
     cout<<DA;
   }
   

   //////////////////////////////////////////////////
   ////////////////// Arret du Fit //////////////////
   // si Lambda = 0, le fit a converge on s'arrete
   //                ou bien on a trop d'iterations
   //////////////////////////////////////////////////
   if(Lambda == 0 || nStep > maxStep) {
     // trop d'iterations
     if(nStep>maxStep) 
       cout<<"GeneralFit : pas de convergence"<<endl;
     // Probleme de matrice de covariance non-definie positive?
     bool bad_covar = false;
     {for(int i=0; i<mNPar; i++) {
       if( fixParam[i] ) errParam(i) = 0.;
       else {
         stepParam(i) = DA(i);
         if( COVAR(i,i)<=0. ) {
            if( debugLevel>0 )
              cout<<"Erreur: Par["<<i<<"]="<<param_tr(i)
                  <<" ("<<Param(i)<<") COVAR()="<<COVAR(i,i)
                  <<" step="<<DA(i)<<endl;
             errParam(i) = 0.;
             bad_covar = true;
         } else {
           errParam(i) = sqrt( COVAR(i,i) );
	 }
       }
     }}
     // print de debug pour parametres bornes
     if(debugLevel>= 2) {
       cout<<"param_tr:"<<endl;
       cout<<param_tr;
       cout<<"stepParam_tr:"<<endl;
       cout<<stepParam;
       cout<<"errParam_tr:"<<endl;
       cout<<errParam;
     }
     // Calcul de la matrice des covariances
     {for(int i=0; i<mNPar; i++) {
        for(int j=0; j<mNPar; j++) {
          if( fixParam[i] || fixParam[j] ) {
            // Parametre fixe, on retourne l'identite
            if(i==j) ATGA(i,j) = 1.; else ATGA(i,j) = 0.;
	  } else if( errParam(i)<=0. || errParam(j)<=0.) {
            // parametres avec mauvaise variance, on retourne 0
            ATGA(i,j) = 0;
	  } else {
            // parametres OK
            ATGA(i,j) = COVAR(i,j)/(errParam(i)*errParam(j));
	  }
	}
     }}
     if (debugLevel >= 1) {
       cout<<">>> Matrice des Covariances = \n";
       cout<<ATGA;
     }
     // Calcul du step et de l'erreur finale en tenant
     // compte des parametres bornes
     stepParam = dtr_vers_dp(stepParam,param_tr);
     errParam  = dtr_vers_dp(errParam,param_tr);
     // Print si demande et code de retour.
     if (debugLevel>0 ) PrintFit();
     if(nStep>maxStep) return(-20);
       else if(bad_covar) return(-11);
         else return(nStep);
   }
   ////////////////////////////////////////////////////////
   ////////////////// Fin d'Arret du Fit //////////////////
   ////////////////////////////////////////////////////////

   // Gestion des deplacements
   {for (int i=0; i<mNPar; i++) {
     if( fixParam[i] ) { DA(i) = 0; continue;}
     // le premier deplacement ne peut etre plus grand que stepParam
     if( nStep == 1 && fabs(DA(i)) > step_tr(i) ) {
       DA(i) = DA(i) < 0. ? -step_tr(i) : step_tr(i);
       if(debugLevel>1 ) cout<<"Excursion parametre "<<i
                             <<" limitee a "<<DA(i)<<endl;
     }
   }}
   paramTry_tr = param_tr + DA;
   paramTry = tr_vers_p(paramTry_tr);
   dparam = dtr_vers_dp(DA,paramTry_tr);
   dparam /= 2.;
   put_in_limits_for_deriv(paramTry,dparam);
   {for(int i=0; i<mNPar; i++) {
     if( ! boundParam[i] ) continue;
     if(paramTry(i) <= minParam(i)) {
       if(debugLevel>0) cout<<"Parametre "<<i
                            <<" limite au minimum"<<endl;
       Param(i) = minParam(i);
       return(-200-i);
     } else if (paramTry(i) >= maxParam(i)) {
       if(debugLevel>0) cout<<"Parametre "<<i
                            <<" limite au maximum"<<endl;
       Param(i) = maxParam(i);
       return(-300-i);
     }
   }}

   // Nouvel essai
   if(mFunction!=NULL) mFunction->SetDeltaParm(dparam.Data());
     else if(mFuncXi2!=NULL) mFuncXi2->SetDeltaParm(dparam.Data());
   if(debugLevel >= 2) {
     cout<<">>>>>>>>>>> ESSAI avec nouveaux parametres\n";
     cout<<"paramTry:\n";
     cout<<paramTry;
     cout<<"paramTry_tr:\n";
     cout<<paramTry_tr;
     cout<<"dparam:\n";
     cout<<dparam;
   }
   if(mFunction!=NULL) TryFunc(paramTry,paramTry_tr);
     else if(mFuncXi2!=NULL) TryXi2(paramTry,paramTry_tr);

   if (debugLevel >= 2)
     cout<<"step "<<nStep<<" Chi2 : old="<<oldChi2
         <<" new="<<Chi2<<" d="<<Chi2-oldChi2<<endl;
   if(FileStep) write_in_step(Chi2,paramTry);

   // *************************************************************
   // ****************** quelle strategie sur Lambda ???? *********
   // *************************************************************
   if (Chi2 < oldChi2) {
     // ****************** le Chi2 est descendu ******************
     nStop = 0;
     if(nStopLent>0) {
       // Arret special demande, comment se comporte les parametres?
       int k=0;
       for (int i=0; i<mNPar; i++) if( (!fixParam[i]) &&
              (fabs(param_tr(i)-paramTry_tr(i))<Eps(i))) k++;
       if (k==mNParFree) nStopL++; // Tous les parametres ont peu varies
          else nStopL=0;
       if (debugLevel>=2) cout<<k<<" parametres sur "<<mNParFree
                              <<" ont peu varies, nStopL="<<nStopL<<endl;
     } else nStopL = 0;
     // Preparation des parametres pour iteration suivante
     ATGA = ATGA_Try;
     BETA = BETA_Try;
     param_tr = paramTry_tr;
     Param = paramTry;
     Lambda *= 1./Lambda_Fac;
     // Arret ?
     if (oldChi2-Chi2<stopChi2) {
       // arret normal, convergence
        Lambda = 0;
        if (debugLevel >= 2)
          cout<<"Arret>> demande car Chi2 decroit et oldChi2-Chi2= "
	      <<oldChi2-Chi2<<"<"<<stopChi2<<endl;
     } else if (nStopLent>0 && nStopL >= nStopLent) {
       // arret demande par SetStopLent, variation lente des parametres
       Lambda = 0.;
       if (debugLevel >= 2) 
         cout<<"Arret>> demande car Chi2 decroit et nStop(lent)= "
             <<nStopL<<">="<<nStopLent<<endl;
     }
     oldChi2 = Chi2;
     if (debugLevel >= 2) cout<<"Succes essai: Lambda divided by "
                              <<Lambda_Fac<<" -> "<<Lambda<<endl;
   } else {
     // ****************** le Chi2 est remonte ******************
     nStopL = 0;
     if(nStopMx>0 && Chi2-oldChi2<stopChi2SMx) {
       // Il est remonte tres peu, comment se comporte les parametres?
       int k=0;
       for (int i=0; i<mNPar; i++) if( (!fixParam[i]) &&
              (fabs(param_tr(i)-paramTry_tr(i))<Eps(i))) k++;
       if (k==mNParFree) nStop++; // Tous les parametres ont peu varies
          else nStop=0;
       if (debugLevel>=2) cout<<k<<" parametres sur "<<mNParFree
                              <<" ont peu varies, nStop="<<nStop<<endl;
     } else nStop = 0;
     // Preparation des parametres pour iteration suivante
     Lambda *= Lambda_Fac;
     // Arret ?
     if (nStopMx>0 && nStop>=nStopMx) {
       // arret normal, convergence car ci2 varie peu et parametres aussi
       Lambda = 0.;
       if (debugLevel >= 2) 
         cout<<"Arret>> demande car Chi2 croit et nstop= "
             <<nStop<<">="<<nStopMx<<endl;
     }
     Chi2 = oldChi2;
     if (debugLevel >= 2)
       cout<<"Echec essai: Lambda multiplied by "<<Lambda_Fac
           <<" -> "<<Lambda<<" nStop="<<nStop<<endl;
   }

 } // fin des iterations
}

//////////////////////////////////////////////////////////////////////
//++
double GeneralFit::ReCalChi2(int& nddl, double *par)
//
//	Recalcul du Chi2 a partir des parametres courants (`par==NULL')
//	ou a partir du tableau de parametres `par'.
//	Retourne le chi2 et le nombre de degres de liberte.
//	Si nddl<0 probleme.
//--
{
double c2 = -.1;
if(par==NULL) par = Param.Data();
if( mData->NData() <= 0 ) {nddl = -100; return 0.;}

if( mFunction != NULL ) {

  double e,result;

  nddl = 0; c2  = 0.;
  for(int k=0; k<mData->NData(); k++) {
    if (! mData->mOK[k]) continue;
    e = mData->mErr[k];
    result = mFunction->Value(&mData->mXP[mNVar*k],par);
    c2 += (mData->mF[k]-result)*(mData->mF[k]-result)/(e*e);
    nddl++;
  }
  nddl -= mNParFree;

  return c2;

} else if( mFuncXi2 != NULL ) {

  c2 = mFuncXi2->Value(*mData,par,nddl);
  nddl -= mNParFree;
  return c2;

} else {

  cout<<"GeneralFit::ReCalChi2_Erreur: mFunction && mFuncXi2 == NULL"<<endl;
  nddl = -1;
  return c2;
}

}

//////////////////////////////////////////////////////////////////////
//++
GeneralFitData* GeneralFit::DataResidus(bool clean)
//
//	Retourne une structure ``GeneralFitData'' contenant
//	les residus du fit (val-func) pour les points du fit.
//	Si ``clean'' est ``true''
//	seules les donnees valides de ``data'' sont copiees.
//	Si ``clean'' est ``false'' (defaut) toutes les donnees
//	sont copiees et la taille totale de ``data'' est allouee
//	meme si elle est plus grande que la taille des donnees stoquees.
//--
{
if(!mData) return NULL;
if(!mFunction) return NULL;
GeneralFitData* datres = new GeneralFitData(*mData,clean);
for(int k=0; k<mData->NData(); k++)
  datres->mF[k] -= mFunction->Value(&datres->mXP[mNVar*k],Param.Data());
return datres;
}

//////////////////////////////////////////////////////////////////////
//++
GeneralFitData* GeneralFit::DataFunction(bool clean)
//
//	Retourne une structure ``GeneralFitData'' contenant
//	les valeurs de la fonction fittee pour les points du fit.
//	(voir commentaires pour ``clean'' dans ``DataResidus'')
//--
{
if(!mData) return NULL;
if(!mFunction) return NULL;
GeneralFitData* datres = new GeneralFitData(*mData,clean);
for(int k=0; k<mData->NData(); k++)
  datres->mF[k] = mFunction->Value(&datres->mXP[mNVar*k],Param.Data());
return datres;
}

//////////////////////////////////////////////////////////////////////
//++
void GeneralFit::PrintFitErr(int rc)
//
//	Imprime le commentaire lie a l'erreur rc retournee par Fit()
//	(voir le commentaire de la methode `Fit()')
//--
{
int n;
if(rc>0) return;

if(rc==-1)
  cout<<"rc = "<<rc<<"  : le nombre de degre de liberte est <0"<<endl;

else if(rc==-10)
  cout<<"rc = "<<rc<<" : l'inversion de la matrice des erreurs n'a pas ete possible"<<endl;

else if(rc==-11)
  cout<<"rc = "<<rc<<" : un element diagonal de la matrice des covariances est <=0"<<endl;

else if(rc==-20)
  cout<<"rc = "<<rc<<" : le fit n'a pas converge (nstep>nNstepMX="<<maxStep<<")"<<endl;

else if(rc>-200 && rc<=-100) {
  n = -100-rc;
  cout<<"rc = "<<rc<<" : le parametre "<<n<<" ("<<nameParam[n]
      <<") est initialise hors limites"<<endl;
}

else if(rc>-300 && rc<=-200) {
  n = -200-rc;
  cout<<"rc = "<<rc<<" : le parametre "<<n<<" ("<<nameParam[n]
      <<") atteint sa limite inferieure"<<endl;
}

else if(rc>-400 && rc<=-300) {
  n = -300-rc;
  cout<<"rc = "<<rc<<" : le parametre "<<n<<" ("<<nameParam[n]
      <<") atteint sa limite superieure"<<endl;
}

else cout<<"rc = "<<rc<<" : type d'erreur inconnue"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Fonctions privees
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
void GeneralFit::write_in_step(double ci2,Vector& par)
{
if(FileStep==NULL) return;
fprintf(FileStep,"%d %d %f",mNtry,nStep,ci2);
for(int i=0; i<mNPar; i++) fprintf(FileStep," %f",par(i));
fprintf(FileStep,"\n");
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::TryFunc(Vector& par,Vector& par_tr)
{
  BETA_Try = 0;
  ATGA_Try = 0;
  Chi2  = 0;
  Vector deriv(mNPar);
  Vector derivtr(mNPar);
  double result;

  for(int k=0; k<mData->NData(); k++) {
    if (! mData->mOK[k]) continue;
    double e = mData->mErr[k];
    if(mNParBound==0)
       result = mFunction->Val_Der(&mData->mXP[mNVar*k]
                                  ,par.Data(),derivtr.Data());
    else {
       result = mFunction->Val_Der(&mData->mXP[mNVar*k]
                                  ,par.Data(),deriv.Data());
       dtr_vers_dp(deriv,par_tr,derivtr);
    }
    double Gkk = 1/(e*e);
    double Ck  = mData->mF[k] - result;
    Chi2 += Ck*Ck*Gkk;
    for(int j=0; j<mNPar; j++) {
      if( fixParam[j] ) continue;
      for(int i=0; i<mNPar; i++)
        if(!fixParam[i]) ATGA_Try(i,j) += derivtr(i)*Gkk*derivtr(j);
      BETA_Try(j) += derivtr(j) * Gkk * Ck;
    }
  }

  if (debugLevel >= 3) {
    cout<<"Try: matrice ( At * G * A )_Try\n";
    cout<<ATGA_Try;
    cout<<"Try: beta_Try:\n";
    cout<<BETA_Try;
  }
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::TryXi2(Vector& par,Vector& par_tr)
{
  double c, *parloc;
  BETA_Try = 0;
  ATGA_Try = 0;
  Chi2  = 0;

  parloc = par.Data();  // He oui, encore ces ... de const*
  Chi2 = mFuncXi2->Value(*mData,parloc,mNddl);
  mNddl -= mNParFree;

  // Calcul des derivees du Xi2 (vecteur du gradient)
  {for(int i=0;i<mNPar; i++) {
    if( fixParam[i] ) continue;
    c = c_dtr_vers_dp(i,par_tr(i));
    BETA_Try(i) = -0.5 * mFuncXi2->Derivee(*mData,i,parloc) * c;
  }}

  // Calcul des derivees 2sd du Xi2 (matrice de courbure ou 0.5*Hessien)
  double c1,c2;
  {for(int i=0;i<mNPar; i++) {
    if( fixParam[i] ) continue;
    c1 = c_dtr_vers_dp(i,par_tr(i));
    for(int j=0;j<mNPar; j++) {
      if( fixParam[j] ) continue;
      c2 = c_dtr_vers_dp(j,par_tr(j));
      ATGA_Try(i,j) = 0.5 * mFuncXi2->Derivee2(*mData,i,j,parloc) *c1*c2;
    }
  }}
  // et on symetrise car d/di(dC2/dj) =  d/dj(dC2/di) mathematiquement
  // mais malheureusement pas numeriquement.
  if( mNPar>1) {
    for(int i=0;i<mNPar-1; i++) {
      if( fixParam[i] ) continue;
      for(int j=i+1;j<mNPar; j++) {
        if( fixParam[j] ) continue;
        c1 = 0.5*(ATGA_Try(i,j) + ATGA_Try(j,i));
        ATGA_Try(i,j) = c1;
        ATGA_Try(j,i) = c1;
      }
    }
  }
  
  if (debugLevel >= 3) {
    cout<<"Try: matrice ( At * G * A )_Try\n";
    cout<<ATGA_Try;
    cout<<"Try: beta_Try:\n";
    cout<<BETA_Try;
  }
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::CheckSanity()
{
  //  DBASSERT( mData != NULL );
  //  DBASSERT( mFunction != NULL || mFuncXi2 != NULL );
  if( mFunction != NULL ) {
    //    DBASSERT( mFunction->NVar() == mNVar );
    //    DBASSERT( mData->NVar() == mNVar );
  }
  //  DBASSERT( mNParFree > 0 && mNParFree <= mNPar );
  //  DBASSERT( mNParBound >= 0 && mNParBound <= mNPar );
  //  DBASSERT( mNParFree <= mData->NDataGood() );
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::Set_Bound_C_D(int i)
// C = (min+max)/2
// D = (max-min)/Pi
{
  // DBASSERT(i>=0 && i<mNPar);
 C(i) = D(i) = 0.;
 if( !boundParam[i] || fixParam[i] ) return;
 C(i) = (maxParam(i)+minParam(i))/2.;
 D(i) = (maxParam(i)-minParam(i))/M_PI;
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::Set_Bound_C_D()
{
 for(int i=0;i<mNPar;i++) Set_Bound_C_D(i);
 if(debugLevel>= 2) {
   cout<<"Set_Bound_C_D: C=\n";
   cout<<C;
   cout<<"Set_Bound_C_D: D=\n";
   cout<<D;
 }
}

//////////////////////////////////////////////////////////////////////
double GeneralFit::p_vers_tr(int i,double p)
// tr = tan( (p-C)/D )
{
 // DBASSERT(i>=0 && i<mNPar);
 double tr = p;
 if(boundParam[i]) tr = tan((p-C(i))/D(i));
 return(tr);
}

//////////////////////////////////////////////////////////////////////
Vector GeneralFit::p_vers_tr(Vector const& p)
{
 Vector tr(p);
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] || ! boundParam[i] ) continue;
   tr(i) = p_vers_tr(i,p(i));
 }
 return(tr);
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::p_vers_tr(Vector const& p,Vector& tr)
{
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] ) continue;
   if( ! boundParam[i] ) tr(i) = p(i);
     else tr(i) = tan((p(i)-C(i))/D(i));
 }
}

//////////////////////////////////////////////////////////////////////
double GeneralFit::tr_vers_p(int i,double tr)
// p = C+D*atan(tr)
{
 // DBASSERT(i>=0 && i<mNPar);
 double p = tr;
 if(boundParam[i]) p = C(i)+D(i)*atan(tr);
 return(p);
}

//////////////////////////////////////////////////////////////////////
Vector GeneralFit::tr_vers_p(Vector const& tr)
{
 Vector p(tr);
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] || ! boundParam[i] ) continue;
   p(i) = tr_vers_p(i,tr(i));
 }
 return(p);
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::tr_vers_p(Vector const& tr,Vector& p)
{
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] ) continue;
   if( ! boundParam[i] ) p(i) = tr(i);
     else p(i) = C(i)+D(i)*atan(tr(i));
 }
}

//////////////////////////////////////////////////////////////////////
double GeneralFit::c_dp_vers_dtr(int i,double tr)
// dtr = (1+tr**2)/D * dp = (1+tan( (p-C)/D )**2)/D * dp = coeff * dp
// attention: df/dp = (1+tr**2)/D * dF/dtr = coeff * dF/dtr
{
 // DBASSERT(i>=0 && i<mNPar);
 double coeff = 1.;
 if(boundParam[i]) coeff = (1.+tr*tr)/D(i);
 return(coeff);
}

//////////////////////////////////////////////////////////////////////
Vector GeneralFit::dp_vers_dtr(Vector const& dp,Vector const& tr)
{
 Vector dtr(dp);
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] || ! boundParam[i] ) continue;
   dtr(i) *= c_dp_vers_dtr(i,tr(i));
 }
 return(dtr);
}

//////////////////////////////////////////////////////////////////////
void GeneralFit::dp_vers_dtr(Vector const& dp,Vector const& tr,Vector& dtr)
{
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] ) continue;
   if( ! boundParam[i] ) dtr(i) = dp(i);
     else dtr(i) = (1.+tr(i)*tr(i))/D(i) * dp(i);
 }
}

//////////////////////////////////////////////////////////////////////
double GeneralFit::c_dtr_vers_dp(int i,double tr)
// dp = D/(1+tr**2) * dtr = coeff * dtr
// attention: df/dtr = D/(1+tr**2) * dF/dp = coeff * dF/dp
{
  // DBASSERT(i>=0 && i<mNPar);
 double coeff = 1.;
 if(boundParam[i]) coeff = D(i)/(1.+tr*tr);
 return(coeff);
}

//////////////////////////////////////////////////////////////////////
Vector GeneralFit::dtr_vers_dp(Vector const& dtr,Vector const& tr)
{
 Vector dp(dtr);
 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] || ! boundParam[i] ) continue;
   dp(i) *= c_dtr_vers_dp(i,tr(i));
 }
 return(dp);
}

//////////////////////////////////////////////////////////////////////
// inline fonction pour aller + vite dans le try()
//void GeneralFit::dtr_vers_dp(Vector const& dtr,Vector const& tr,Vector& dp)

//////////////////////////////////////////////////////////////////////
int GeneralFit::put_in_limits_for_deriv(Vector const& p,Vector& dp,double dist)
// 1-/ Redefinit dp pour qu'il soit superieur a minStepDeriv
// 2-/ Redefinit dp pour que p+/-dp reste dans les limites (parametre borne)
// Si hors limites alors:
//     p-dp <= min_p : dp = (p-min_p)*dist
//     p+dp >= max_p : dp = (max_p-p)*dist
{
 int nchanged = 0;
 bool changed;
 double dp_old;

 for(int i=0;i<mNPar;i++) {
   if( fixParam[i] ) {dp(i)=0.; continue;} // Pas calcul derivee pour param fixe

   if( fabs(dp(i))<minStepDeriv(i) ) {
     // On ne redefinit dp que si minStepDeriv>0.
     dp_old = dp(i);
     if(dp(i)>=0.) dp(i) = minStepDeriv(i); else dp(i) = -minStepDeriv(i);
     if(debugLevel>=2)
       cout<<"put_in_limits_for_deriv(range) dp["<<i<<"]=abs("<<dp_old
           <<") <"<<minStepDeriv(i)<<" changed to "<<dp(i)<<endl;
   }

   if( !boundParam[i] ) continue;

   changed = false;
   if( p(i)-dp(i)<=minParam(i) ) {
     dp_old = dp(i);
     dp(i) = dist*(p(i)-minParam(i));
     changed = true;
     if(debugLevel>=2)
       cout<<"put_in_limits_for_deriv(min) p["<<i<<"}="<<p(i)<<" >="
           <<minParam(i)<<" .. dp="<<dp_old<<" -> dp="<<dp(i)<<endl;
   }

   if( p(i)+dp(i)>=maxParam(i) ) {
     dp_old = dp(i);
     dp(i) = dist*(maxParam(i)-p(i));
     changed = true;
     if(debugLevel>=2)
       cout<<"put_in_limits_for_deriv(max) p["<<i<<"]="<<p(i)<<" <="
           <<maxParam(i)<<" .. dp="<<dp_old<<" -> dp="<<dp(i)<<endl;
   }

   if(changed) nchanged++;
 }

 return nchanged;
}


//////////////////////////////////////////////////////////////////////
// Rappel des inline functions pour commentaires
//++
// inline double   GetChi2()
//	Retourne le Chi2
//--
//++
// inline double   GetChi2Red() const
//	Retourne le Chi2 reduit
//--
//++
// inline int      GetNddl()    const
//	Retourne le nombre de degres de liberte
//--
//++
// inline int      GetNStep()   const
//	Retourne le nombre d'iterations
//--
//++
// inline int      GetNVar()    const
//	Retourne le nombre de variables
//--
//++
// inline int      GetNPar()    const
//	Retourne le nombre de parametres
//--
//++
// inline int      GetNFree()   const
//	Retourne le nombre de parametres libres
//--
//++
// inline int      GetNBound()  const
//	Retourne le nombre de parametres bornes
//--
//++
// inline int      GetNStop()   const
//	Retourne le nstop de convergence
//--
//++
// inline int      GetNStopLent()   const
//	Retourne le nstop de convergence lente.
//--
//++
// inline double   GetEps(int i)
//	Retourne la precision de convergence pour le parametre i.
//--
//++
// inline GeneralFunction*  GetFunction()
//	Retourne le pointeur sur la GeneralFunction utilisee.
//--
//++
// inline GeneralFitData*   GetGData()
//	Retourne le pointeur sur la GeneralFitData utilisee.
//--

#endif /* IS_IT_USEFUL */
