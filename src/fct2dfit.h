#ifndef FCT2DFIT_SEEN
#define FCT2DFIT_SEEN

#include "generalfit.h"

//================================================================
// GeneralFunction 2D pour PSF pixel taille 1x1
//================================================================

class GeneralPSF2D : public GeneralFunction {
public:
  GeneralPSF2D(unsigned int nPar);
  virtual ~GeneralPSF2D();

  virtual double ValueH(double const xp[], double const* parm);
  virtual double VolPSF(double const* parm);
  virtual void   DefaultParam(double *parm);

  void SetVolEps(double const prec) ;
protected:
  double VolEps;
 private:
   GeneralPSF2D & operator=( const GeneralPSF2D &);
};

//================================================================
// GeneralFunction 2D pour MULTI-PSF pixel taille 1x1
//================================================================

class GenMultiPSF2D : public GeneralPSF2D {
public:
  GenMultiPSF2D(GeneralPSF2D* psf2d,unsigned int nstar);
  virtual ~GenMultiPSF2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);

protected:
  GeneralPSF2D* mPsf2D; // Type de PSF generique a fiter.
  int mNStar;           // Nombre d etoiles a fiter.
  int mNParmTot;        // Nombre total de parametres pour le fit.
  int mNParm;           // Nombre de parametre pour la PSF generique.
  int mNForme;          // Nombre de parametres de forme autres que Sx,Sy,Rho
  double* mParm;        // pour la PSF generique. Buffer pour PSF generique
  double* mDer;         // et pour le stoquage temporaire des derivees.

private:
  GenMultiPSF2D(const GenMultiPSF2D &);
  GenMultiPSF2D &  operator=( const GenMultiPSF2D &);

};

//==============================================================================
// CLASSES DE FONCTIONS 2D type PSF AVEC PARAMETRES POUR LE FIT pixel taille 1x1
//==============================================================================

//////////////////////////////////////////////////////////////////
class GauRho2D : public GeneralPSF2D {
public:
  GauRho2D();
  virtual ~GauRho2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class GauRhInt2D : public GeneralPSF2D {
public:
  GauRhInt2D();
  virtual ~GauRhInt2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class GdlRho2D : public GeneralPSF2D {
public:
  GdlRho2D();
  virtual ~GdlRho2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class GdlRhInt2D : public GeneralPSF2D {
public:
  GdlRhInt2D();
  virtual ~GdlRhInt2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class Gdl1Rho2D : public GeneralPSF2D {
public:
  Gdl1Rho2D();
  virtual ~Gdl1Rho2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class Gdl1RhInt2D : public GeneralPSF2D {
public:
  Gdl1RhInt2D();
  virtual ~Gdl1RhInt2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class Gdl2Rho2D : public GeneralPSF2D {
public:
  Gdl2Rho2D();
  virtual ~Gdl2Rho2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class Gdl2RhInt2D : public GeneralPSF2D {
public:
  Gdl2RhInt2D();
  virtual ~Gdl2RhInt2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class MofRho2D : public GeneralPSF2D {
public:
  MofRho2D();
  virtual ~MofRho2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};

//////////////////////////////////////////////////////////////////
class MofRhInt2D : public GeneralPSF2D {
public:
  MofRhInt2D();
  virtual ~MofRhInt2D();

  virtual double Value(double const xp[], double const* Par);
  virtual double Val_Der(double const xp[],double const* Par,double* DgDpar);
  virtual double ValueH(double const xp[], double const* Par);
  virtual double VolPSF(double const* Par);
  virtual void   DefaultParam(double *Par);
};


#ifdef IS_IT_USEFUL
//==============================================================================
// CLASSES DE FONCTIONS 2D type Xi2 AVEC PARAMETRES POUR LE FIT pixel taille 1x1
//==============================================================================

//////////////////////////////////////////////////////////////////
class X2_GauRho2D : public GeneralXi2 {
public:
  X2_GauRho2D();
  virtual ~X2_GauRho2D();

  virtual double Value(GeneralFitData& data, double* parm, int& ndataused);
  virtual double Derivee2(GeneralFitData& data, int i,int j, double* parm);

protected:
  GauRho2D* gaurho2d;
};

#endif
#endif
