//#ifdef IS_IT_USEFUL
#include "defs.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "fct2dfit.h"
//#include "perrors.h"
#include "nbconst.h"
//#include "tabmath.h"

// define SIMPSON4  c'etait la prod 91-95 rcecile
#define SIMPSON9
#include "simps2d.h"

# define EXPO exp
//#define EXPO tabFExp
#define MINEXPM (100.)

using namespace std;

//================================================================
// GeneralFunction 2D pour PSF pixel taille 1x1
//================================================================

/////////////////////////////////////////////////////////////////
//++
// Class	GeneralPSF2D
// Lib	Outils++ 
// include	fct2dfit.h
//
//	Classe de definition d'une PSF 2D a nPar parametres
//	Pour definir une PSF, il faut creer une classe qui herite
//	de ``GeneralPSF2D'' (cf par exemple GauRho2D...).
//	La disposition des parametres definissant la PSF est indifferente,
//	toutefois il est conseille de suivre l'ordre:
//--
//++
//  - PSF 2D a NPar parametres:
//  p[0] = Volume (ou hauteur)
//  p[1] = centre X0, p[2] = centre Y0
//  p[3] = SigmaX   , p[4] = SigmaY,    p[5] = Rho
//  p[6],p[7],... = autres parametres (eventuels) definissant la PSF.
//                  (ex: pour la Moffat p[6] = exposant Beta et NPar=8).
//  p[NPar-1] = Fond
//--
//++
//	L'emploi de certaines classes comme par exemple ``GenMultiPSF2D''
//	necessite de suivre rigoureusement l'ordre indique ci-dessus
//	pour les parametres.
//--

//++
GeneralPSF2D::GeneralPSF2D(unsigned int nPar)
//
//--
: GeneralFunction(2,nPar), VolEps(1.e-4)
{
  //DBASSERT( nPar>0 );
}

GeneralPSF2D::~GeneralPSF2D()
{
}

//++
double GeneralPSF2D::ValueH(double const xp[], double const* parm)
//
//| ValueH = hauteur*forme(x,y)+fond tq forme(0,0)=1.
//|     alors que Value = volume*forme(x,y)+fond tq volume(forme)=1.
//| Dans notre convention le dernier parametre est le fond,
//| le premier le volume et les 2 suivants le centrage x0,y0
//| ---> Ici parm[0] = hauteur
//--
{
double x0[2];
int mm1 = mNPar - 1;

// point central en [x0,y0]
x0[0] = parm[1];  x0[1] = parm[2];

// retour avec hauteur = 1
return   (Value(xp,parm) - parm[mm1]) / (Value(x0,parm) - parm[mm1])
         * parm[0] + parm[mm1];
}

//++
double GeneralPSF2D::VolPSF(double const* parm)
//
//| Cette fonction calcule le volume d'une PSF de hauteur=1
//| avec une precision de "VolEps"
//| dans le but de connaitre le coefficient permettant
//| de convertir le volume d'une PSF en son amplitude
//| ou vice-versa: " volume = VolPSF * hauteur "
//|  L'integration se fait 1/4 de pixel par 1/4 de pixel
//|  ATTENTION: Il s'agit de PSF donc x,y,x0,y0,Sigma.. sont en pixels
//--
{
 double x[2],step;
 double vol,volprec;
 int ecart,i,j,k;
 int mm1 = mNPar-1;

 step = 1. / 4.;
 vol = volprec = 0.;
 ecart = 1;

 /* pixel central */
 for(k=0;k<nd2d;k++) {
   x[0] = parm[1] + dx2d[k]*step;
   x[1] = parm[2] + dy2d[k]*step;
   vol += (ValueH(x,parm)-parm[mm1]) * w2d[k];
 }

/* increment en couronnes carrees de 2*ecart+1 de cote */
 while ( ecart < 2 || fabs((vol-volprec)/vol) > VolEps ) {
   volprec = vol;
   for (i= -ecart;i<=ecart;i++) for(k=0;k<nd2d;k++) {
     x[0] = parm[1] + (i+dx2d[k])*step;
     x[1] = parm[2] + (-ecart+dy2d[k])*step;
     vol += (ValueH(x,parm)-parm[mm1]) * w2d[k];
     x[1] = parm[2] + ( ecart+dy2d[k])*step;
     vol += (ValueH(x,parm)-parm[mm1]) * w2d[k];
   }
   for (j= -ecart+1;j<=ecart-1;j++) for(k=0;k<nd2d;k++) {
     x[1] = parm[2] + (j+dy2d[k])*step;
     x[0] = parm[1] + (-ecart+dx2d[k])*step;
     vol += (ValueH(x,parm)-parm[mm1]) * w2d[k];
     x[0] = parm[1] + ( ecart+dx2d[k])*step;
     vol += (ValueH(x,parm)-parm[mm1]) * w2d[k];
   }
   ecart++;
   // printf("ec=%d v=%f prec=%f %f\n",ecart,vol,fabs((vol-volprec)/vol),VolEps);
 }

 vol *= step * step / parm[0];
 return vol;
}

//++
void GeneralPSF2D::DefaultParam(double *parm)
//
//	Definition des defauts des parametres
//--
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
}

//++
void GeneralPSF2D::SetVolEps(double const prec)
//
//	Definition de la precision sur le calcul du volume
//--
{
  VolEps = prec;
}

//================================================================
// GeneralFunction 2D pour MULTI-PSF pixel taille 1x1
//================================================================

/////////////////////////////////////////////////////////////////
//++
// Class	GenMultiPSF2D
// Lib	Outils++ 
// include	fct2dfit.h
//
//	Classe de definition d'un ensemble de PSF2D
//	pour fiter simultanement plusieurs etoiles et un fond constant.
//	Les parametres de forme de la PSF (Sx, Sy, Rho etc... et Fond)
//	sont les memes pour toutes les etoiles, seuls le centre
//	(X0,Y0) et le volume (ou la hauteur) V varient pour chaque etoile.
//	La disposition des parametres definissant la PSF generique
//	est obligatoirement la suivante:
//--
//++
//|  - PSF 2D a NPar parametres:
//|  p[0] = Volume (ou hauteur)
//|  p[1] = centre X0, p[2] = centre Y0
//|  p[3] = SigmaX   , p[4] = SigmaY,    p[5] = Rho
//|  p[6],p[7],... = autres parametres (eventuels) definissant la PSF.
//|                  (ex: pour la Moffat p[6] = exposant Beta et NPar=8).
//|  p[NPar-1] = Fond
//|
//--
//++
//|  - La Multi-PSF a ses parametres arranges dans l'ordre suivant:
//|    Soit NStar le nombre d'etoiles a fiter simultanement
//|         NP = le nombre de parametres de la PSF 2D generique
//|    On a NF = NP-7 parametres de forme supplementaires
//|         (ex: nf=0 pour GauRho2D, nf=1 pour MofRho2D)
//|  p[0],p[1],p[2] = V0,X0,Y0 pour la premiere etoile
//|  p[3],p[4],p[5] = V1,X1,Y1 pour la deuxieme etoile
//|  ...
//|  p[3*i],p[3*i+1],p[3*i+2] = Vi,Xi,Yi pour la (i+1) ieme etoile
//|  ...
//|  p[m*i],p[m*i+1],p[m*i+2] = Vm,Xm,Ym   ;   m = NStar-1
//|                       pour la NStar ieme et derniere etoile
//|  p[3*NStar],p[3*NStar+1],p[3*NStar+2] = SigmaX, SigmaY et Rho
//|  p[3*NStar+3],...,p[3*NStar+2+NF] = parametres de forme
//|                        supplementaires pour definir la PSF 2D
//|  p[3*NStar+2+NF+1] = Fond
//--

//++
GenMultiPSF2D::GenMultiPSF2D(GeneralPSF2D* psf2d,unsigned int nstar)
//
//	Createur. ``psf2d'' est le nom de la PSF generique a utiliser,
//	et ``nstar'' est le nombre d'etoiles a fiter simultanement.
//--
  : GeneralPSF2D((psf2d!=NULL) ? 3*nstar+4+psf2d->NPar()-7: 0)
  , mPsf2D(psf2d), mNStar(nstar)
{
  //DBASSERT( nstar>0 && psf2d!=NULL );
mNForme = mPsf2D->NPar() - 7;
//DBASSERT( mNForme>=0 );
mNParm = mPsf2D->NPar();
mParm = new double[mNParm];
mDer = new double[mNParm];
mNParmTot = GeneralPSF2D::NPar();
cout<<"mNStar="<<mNStar<<" mNParmTot="<<mNParmTot
    <<" mNParm="<<mNParm<<" mNForme="<<mNForme<<endl;
}

GenMultiPSF2D::~GenMultiPSF2D()
{
delete [] mParm; mParm = NULL;
delete [] mDer;  mDer = NULL;
}

double GenMultiPSF2D::Value(double const xp[], double const* Par)
{
// Fond commun
double val = Par[mNParmTot-1];

// Remplissage le tableau des parametres pour la PSF generique
// ... Communs a toutes les PSF individuelles: Sx,Sy,Rho,[Forme],Fond
const double *pt = &Par[3*mNStar];
double *p = &mParm[3];
{for(int i=0;i<3+mNForme;i++) *(p++) = *(pt++);}  // Sx,Sy,Rho,[Forme...]
*(p++) = 0.;   // Fond

// ... Propres a chaque etoiles: Vi,Xi,Yi
pt = Par;
{for(int i=0;i<mNStar;i++) {
  mParm[0] = *(pt++);  // Vi (ou Hi)
  mParm[1] = *(pt++);  // Xi
  mParm[2] = *(pt++);  // Yi
  val += mPsf2D->Value(xp,mParm);
}}

return val;
}

double GenMultiPSF2D::Val_Der(double const xp[], double const* Par
                             ,double *DgDpar)
{
{for(int i=3*mNStar;i<mNParmTot-1;i++) DgDpar[i] = 0.;}

// Fond commun
double val = Par[mNParmTot-1];
DgDpar[mNParmTot-1] = 1.;  // D./DFond

// Remplissage le tableau des parametres pour la PSF generique
// ... Communs a toutes les PSF individuelles: Sx,Sy,Rho,[Forme],Fond
const double *pt = &Par[3*mNStar];
double *p = &mParm[3];
{for(int i=0;i<3+mNForme;i++) *(p++) = *(pt++);}  // Sx,Sy,Rho,[Forme...]
*(p++) = 0.;   // Fond

// ... Propres a chaque etoiles: Vi,Xi,Yi
double *dpt = DgDpar, *dpt2 = &DgDpar[3*mNStar];
pt = Par;
{for(int i=0;i<mNStar;i++) {
  mParm[0] = *(pt++);  // Vi (ou Hi)
  mParm[1] = *(pt++);  // Xi
  mParm[2] = *(pt++);  // Yi
  val += mPsf2D->Val_Der(xp,mParm,mDer);
  {for(int j=0;j<3;j++) *(dpt++) = mDer[j];}  // D./DVi,D./DXi,D./DYi
  {for(int j=0;j<3+mNForme;j++) *(dpt2+j) += mDer[3+j];} // D./DSx,D./DSy,D./DRho,[D./DForme]
}}

return val;
}

//==============================================================================
// CLASSES DE FONCTIONS 2D type PSF AVEC PARAMETRES POUR LE FIT pixel taille 1x1
// la taille du pixel est importante quand on utilise les PSF integrees
//    (x,y x0,y0 sigmaX.... sont en unites de pixels !!!)
//==============================================================================

#define _x0_   Par[1]
#define _y0_   Par[2]
#define _sigx_ Par[3]
#define _sigy_ Par[4]
#define _rho_  Par[5]
#define _Gm_   Par[6]
#define _B4_   Par[6]
#define _B6_   Par[7]
#define _B2_   Par[8]

/////////////////////////////////////////////////////////////////
//++ 
// Module	Classes de PSF 2D
// Lib	Outils++ 
// include	fct2dfit.h
//--
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
//++
// Titre	GauRho2D
// \index{GauRho2D}
//
//| gaussienne+fond 2D
//| Par [0]=vol [1]=x0 [2]=y0 [3]=sigx [4]=sigy [5]=rho [6]=fond
//|   sigx,sigy,rho = sigma et rho de la gaussienne
//|   x0,y0 = centre de la gaussienne
//|  PSF(x,y) = N * exp[ - 1/2 (X**2 + Y**2 -2*rho*X*Y) ]
//|           avec X = (x-x0)/sigx et Y = (y-y0)/sigy
//|                N = sqrt(1-rho**2)/(2*Pi*sigx*sigy)
//|  le volume de cette gaussienne est V=1.
//|  F(x,y) = Par[0]*PSF(x,y)+Par[6] (volume=Par[0],fond=Par[6])
//--
//++
//| -*- Remarque: De la facon dont est ecrite la PSF gaussienne
//| sigx,sigy representent les sigmas des gaussiennes 1D
//| qui sont les coupes de la gaussienne 2D pour y=0 et x=0.
//| Les moments centres d'ordre 2 sont
//|   sx = sigx/sqrt(1-ro^2) et sy = sigy/sqrt(1-ro^2)
//--
/////////////////////////////////////////////////////////////////

//++
GauRho2D::GauRho2D()
//
//	Createur
//--
: GeneralPSF2D(7)
{
}

GauRho2D::~GauRho2D()
{
}

double GauRho2D::Value(double const xp[], double const* Par)
{
 double N = sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 if( z2<MINEXPM ) return Par[0] * N*EXPO(-z2) + Par[6];
   else return Par[6];
}

double GauRho2D::ValueH(double const xp[], double const* Par)
{
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 //double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 double z2 = 0.5*(X-Y)*(X-Y) - (_rho_ - 1)*X*Y;
 if( z2<MINEXPM ) return Par[0] * EXPO(-z2) + Par[6];
   else return Par[6];
}

double GauRho2D::VolPSF(double const* Par)
{
 return DeuxPi * _sigx_ * _sigy_ / sqrt(1.-_rho_*_rho_);
}

double GauRho2D::Val_Der(double const xp[], double const* Par
                        , double *DgDpar)
{
 double unmr2 = 1.-_rho_*_rho_;
 double N = sqrt(unmr2)/(DeuxPi*_sigx_*_sigy_);
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;

 double XmrY = X-_rho_*Y;
 double YmrX = Y-_rho_*X;
 double z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;

 /* g(x,y) */
 double PSF = 0.;
 if( z2<MINEXPM ) PSF = N * EXPO(-z2);
 /* dg(x,y)/d(x0) */
 DgDpar[1] =  Par[0]* PSF* XmrY/_sigx_;
 /* dg(x,y)/d(y0) */
 DgDpar[2] =  Par[0]* PSF* YmrX/_sigy_;
 /* dg(x,y)/d(sx)*/
 DgDpar[3] =  Par[0]* PSF* (X*XmrY-1.)/_sigx_;
 /* dg(x,y)/d(sy) */
 DgDpar[4] =  Par[0]* PSF* (Y*YmrX-1.)/_sigy_;
 /* dg(x,y)/d(rho) */
 DgDpar[5] =  Par[0]* PSF* (X*Y-2.*_rho_/unmr2);
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = PSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[6] = 1.;

 return Par[0] * PSF + Par[6];
}

void GauRho2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
}


/////////////////////////////////////////////////////////////////
//++
// Titre	GauRhInt2D
// \index{GauRhInt2D}
//
//| Cette fonction calcule une approximation a l'integrale d'une
//| gaussienne 2D sur un carre de longueur unite (x,y-05 -> x,y+0.5)
//--
/////////////////////////////////////////////////////////////////

//++
GauRhInt2D::GauRhInt2D()
//
//	Createur
//--
: GeneralPSF2D(7)
{
}

GauRhInt2D::~GauRhInt2D()
{
}

double GauRhInt2D::Value(double const xp[], double const* Par)
{
  double N = sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
  double xx = xp[0] - _x0_;
 double yy = xp[1] - _y0_;
 double SPSF = 0.;
 double X,Y,z2;
 for(int i=0; i<nd2d; i++) {
   X = (xx+dx2d[i])/_sigx_;
   Y = (yy+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;  
   if( z2<MINEXPM ) 
     {
       SPSF += EXPO(-z2) * w2d[i];
     }
 }
 return Par[0]* N*SPSF + Par[6];
}

double GauRhInt2D::ValueH(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF = 0.;
 double X,Y,z2;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   if( z2<MINEXPM ) SPSF += EXPO(-z2) * w2d[i];
 }
 return Par[0] *SPSF + Par[6];
}

double GauRhInt2D::VolPSF(double const* Par)
{
 return DeuxPi * _sigx_ * _sigy_ / sqrt(1.-_rho_*_rho_);
}

double GauRhInt2D::Val_Der(double const xp[], double const* Par
                          ,double *DgDpar)
{
 for(int i=0; i<=6; i++) DgDpar[i] = 0.;
 double unmr2 = 1.-_rho_*_rho_;
 double N = sqrt(unmr2)/(DeuxPi*_sigx_*_sigy_);
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;

 double z2,PSF,X,Y,XmrY,YmrX;
 double SPSF = 0.;
 {
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   XmrY = X-_rho_*Y;
   YmrX = Y-_rho_*X;
   z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
   /* g(x,y) */
   if(z2<MINEXPM) PSF = N * EXPO(-z2) * w2d[i]; else PSF = 0.;
   SPSF += PSF;
   /* dg(x,y)/d(x0) */
   DgDpar[1] += Par[0] * PSF* XmrY/_sigx_;
   /* dg(x,y)/d(y0) */
   DgDpar[2] += Par[0] * PSF* YmrX/_sigy_;
   /* dg(x,y)/d(sx)*/
   DgDpar[3] += Par[0] * PSF* (X*XmrY-1.)/_sigx_;
   /* dg(x,y)/d(sy) */
   DgDpar[4] += Par[0] * PSF* (Y*YmrX-1.)/_sigy_;
   /* dg(x,y)/d(rho) */
   DgDpar[5] += Par[0] * PSF* (X*Y-2.*_rho_/unmr2);
 }
 }
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = SPSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[6] = 1.;

 return Par[0] *SPSF + Par[6];
}

void GauRhInt2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
}

#define B4 1.
#define B6 1.
#define KB4B6 0.136887

/////////////////////////////////////////////////////////////////
//++
// Titre	GdlRho2D
// \index{GdlRho2D}
//
//| Cette fonction calcule une gaussienne 2D de volume 1 approchee
//|  par son developpement limite au 3ieme ordre (see dophot)
//|  Meme commentaire que GauRho2D, cf plus haut sauf que:
//|  Par [0]=vol [1]=x0 [2]=y0 [3]=sigx [4]=sigy [5]=rho [6]=fond
//|    z**2 = 1/2 (X**2 + Y**2 -2*rho*X*Y)
//|    PSF(x,y) = N / [ 1 + z**2 + B4/2 *z**4 + B6/6 *z**6 ]
//|               N = KB4B6
//|    le coefficient KB4B6 etant trop dur a calculer analytiquement
//|    Il doit etre calcule numeriquement et entre dans ce programme
//|  ATTENTION: dans cette routine B4 et B6 sont imposes et pas fites!
//|  - DL de la gaussienne:  B4=1., B6=1., KB4B6=0.13688679
//|  le volume de cette gaussienne est V=1.
//|  F(x,y) = Par[0]*PSF(x,y)+Par[6] (volume=Par[0],fond=Par[6])
//--
/////////////////////////////////////////////////////////////////

//++
GdlRho2D::GdlRho2D()
//
//	Createur
//--
: GeneralPSF2D(7)
{
}

GdlRho2D::~GdlRho2D()
{
}

double GdlRho2D::Value(double const xp[], double const* Par)
{
 double N = KB4B6*sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_);
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 double D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
 return Par[0] *N/D + Par[6];
}

double GdlRho2D::ValueH(double const xp[], double const* Par)
{
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 double D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
 return Par[0] /D + Par[6];
}

double GdlRho2D::VolPSF(double const* Par)
{
 return _sigx_*_sigy_/(KB4B6*sqrt(1.-_rho_*_rho_));
}

double GdlRho2D::Val_Der(double const xp[], double const* Par
                        , double *DgDpar)
{
 double unmr2 = 1.-_rho_*_rho_;
 double N = KB4B6*sqrt(unmr2)/(_sigx_*_sigy_);
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double XmrY = X-_rho_*Y;
 double YmrX = Y-_rho_*X;
 double z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
 double D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
 double dDsD = (1.+z2*(B4+B6/2.*z2))/D;

 /* g(x,y) */
 double PSF = N / D;
 /* dg(x,y)/d(x0) */
 DgDpar[1] = Par[0]* PSF* dDsD*XmrY/_sigx_;
 /* dg(x,y)/d(y0) */
 DgDpar[2] = Par[0]* PSF* dDsD*YmrX/_sigy_;
 /* dg(x,y)/d(sx)*/
 DgDpar[3] = Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
 /* dg(x,y)/d(sy) */
 DgDpar[4] = Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
 /* dg(x,y)/d(rho) */
 DgDpar[5] = Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = PSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[6] = 1.;

 return Par[0] *PSF + Par[6];
}

void GdlRho2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
}


/////////////////////////////////////////////////////////////////
//++
// Titre	GdlRhInt2D
// \index{GdlRhInt2D}
//
//	fonction integree de  GdlRho2d
//--
/////////////////////////////////////////////////////////////////

//++
GdlRhInt2D::GdlRhInt2D()
//
//	Createur
//--
: GeneralPSF2D(7)
{
}

GdlRhInt2D::~GdlRhInt2D()
{
}

double GdlRhInt2D::Value(double const xp[], double const* Par)
{
 double N = KB4B6*sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_);
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;

 double z2,X,Y,D;
 double SPSF = 0.;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
   SPSF += w2d[i] / D;
  }
 return Par[0] *N*SPSF + Par[6];
}

double GdlRhInt2D::ValueH(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;

 double z2,X,Y,D;
 double SPSF = 0.;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
   SPSF += w2d[i] / D;
  }
 return Par[0] *SPSF + Par[6];
}

double GdlRhInt2D::VolPSF(double const* Par)
{
 return _sigx_*_sigy_/(KB4B6*sqrt(1.-_rho_*_rho_));
}

double GdlRhInt2D::Val_Der(double const xp[], double const* Par
                          , double *DgDpar)
{
 for(int i=0; i<=6; i++) DgDpar[i] = 0.;
 double unmr2 = 1.-_rho_*_rho_;
 double N = KB4B6*sqrt(unmr2)/(_sigx_*_sigy_);
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;

 double z2,PSF,X,Y,XmrY,YmrX,D,dDsD;
 double SPSF = 0.;
 {
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   XmrY = X-_rho_*Y;
   YmrX = Y-_rho_*X;
   z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
   D = 1.+z2*(1.+z2*(B4/2.+B6/6.*z2));
   dDsD = (1.+z2*(B4+B6/2.*z2))/D;
   /* g(x,y) */
   PSF = N / D  * w2d[i];
   SPSF += PSF;
   /* dg(x,y)/d(x0) */
   DgDpar[1] += Par[0]* PSF* dDsD*XmrY/_sigx_;
   /* dg(x,y)/d(y0) */
   DgDpar[2] += Par[0]* PSF* dDsD*YmrX/_sigy_;
   /* dg(x,y)/d(sx)*/
   DgDpar[3] += Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
   /* dg(x,y)/d(sy) */
   DgDpar[4] += Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
   /* dg(x,y)/d(rho) */
   DgDpar[5] += Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
 }
 }
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = SPSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[6] = 1.;

 return Par[0] *SPSF + Par[6];
}

void GdlRhInt2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
}

#undef B4
#undef B6
#undef KB4B6

/////////////////////////////////////////////////////////////////
//++
// Titre	Gdl1Rho2D
// \index{Gdl1Rho2D}
//
//| Cette fonction calcule une gaussienne 2D approchee
//| par son developpement limite au 2sd ordre (see dophot)
//| Meme commentaire que GauRho2D, cf plus haut sauf que:
//|   z**2 = 1/2 (X**2 + Y**2 -2*rho*X*Y)
//|   PSF(x,y) = N / [ 1 + z**2 + B4/2 *z**4 ]
//|   Le coefficient B4 est fitte (6ieme parametres)
//| ATTENTION: les normalisations N dependent de B4
//|  1-/ B4 est suppose etre toujours positif pour que la PSF tendent
//|      vers 0+ quand z2 tend vers l'infini
//|  2-/ Il y a 3 cas de calcul de K(B4) = int(PSF(x,y)) de 0 a l'infini
//|             0<B4<1/2, 1/2<B4, et B4=1/2
//|      mais pour des raisons d'analyse
//|      numerique j'ai pris 3 intervalles:
//|             0.<B4<0.499, 0.501<B4, 0.499<=B4<=0.501
//|      dans le 3ieme intervalle, comme K est continue est derivable
//|      en B4=1/2, j'ai represente K par la droite tangeante
//|      ce qui, apres verification dans paw est une tres bonne approx.
//|      (je tiens les calculs a disposition.. me demander)
//|  Par [0]=vol [1]=x0 [2]=y0 [3]=sigx [4]=sigy [5]=rho [6]=B4 [7]=fond
//|  F(x,y) = Par[0]*PSF(x,y)+Par[7]
//--
/////////////////////////////////////////////////////////////////

//++
Gdl1Rho2D::Gdl1Rho2D()
//
//	Createur
//--
: GeneralPSF2D(8)
{
}

Gdl1Rho2D::~Gdl1Rho2D()
{
}

double Gdl1Rho2D::Value(double const xp[], double const* Par)
{
 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   return(0.);
 }
 double N = sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_*DeuxPi*K);

 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 double D = 1.+z2*(1.+z2*_B4_/2.);
 return Par[0] *N/D + Par[7];
}

double Gdl1Rho2D::ValueH(double const xp[], double const* Par)
{
 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
 double D = 1.+z2*(1.+z2*_B4_/2.);
 return Par[0] /D + Par[7];
}

double Gdl1Rho2D::VolPSF(double const* Par)
{
 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   return(0.);
 }
 double N = sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_*DeuxPi*K);
 return 1./N;
}

double Gdl1Rho2D::Val_Der(double const xp[], double const* Par
                         , double *DgDpar)
{
 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   for(int i=0;i<=7;i++) DgDpar[i] = 0.;
   return(0.);
 }
 double unmr2 = 1.-_rho_*_rho_;
 double N = sqrt(unmr2)/(_sigx_*_sigy_*DeuxPi*K);

 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double XmrY = X-_rho_*Y;
 double YmrX = Y-_rho_*X;
 double z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
 double D = 1.+z2*(1.+z2*_B4_/2.);
 double dDsD = (1.+z2*_B4_)/D;

 /* g(x,y) */
 double PSF = N / D;
 /* dg(x,y)/d(x0) */
 DgDpar[1] =  Par[0]* PSF* dDsD*XmrY/_sigx_;
 /* dg(x,y)/d(y0) */
 DgDpar[2] =  Par[0]* PSF* dDsD*YmrX/_sigy_;
 /* dg(x,y)/d(sx)*/
 DgDpar[3] =  Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
 /* dg(x,y)/d(sy) */
 DgDpar[4] =  Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
 /* dg(x,y)/d(rho) */
 DgDpar[5] =  Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
 /* dg(x,y)/d(B4) */
 DgDpar[6] =  Par[0]* PSF* (-dKdB4/K-z2*z2/2./D);
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = PSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[7] = 1.;

 return Par[0] *PSF + Par[7];
}

void Gdl1Rho2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = 1.;
}

/////////////////////////////////////////////////////////////////
//++
// Titre	Gdl1RhInt2D
// \index{Gdl1RhInt2D}
//
//	fonction integree de Gdl1Rho2D
//--
/////////////////////////////////////////////////////////////////

//++
Gdl1RhInt2D::Gdl1RhInt2D()
//
//	Createur
//--
: GeneralPSF2D(8)
{
}

Gdl1RhInt2D::~Gdl1RhInt2D()
{
}

double Gdl1RhInt2D::Value(double const xp[], double const* Par)
{
 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   return(0.);
 }
 double N = sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_*DeuxPi*K);

 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF=0.;
 double z2,X,Y,D;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y -2.*_rho_*X*Y)/2.;
   D = 1.+z2*(1.+z2*_B4_/2.);
   SPSF += w2d[i] / D;
 }
 return Par[0] *N*SPSF + Par[7];
}

double Gdl1RhInt2D::ValueH(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF=0.;
 double z2,X,Y,D;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y -2.*_rho_*X*Y)/2.;
   D = 1.+z2*(1.+z2*_B4_/2.);
   SPSF += w2d[i] / D;
 }
 return Par[0] *SPSF + Par[7];
}

double Gdl1RhInt2D::VolPSF(double const* Par)
{
 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   return(0.);
 }
 double N = sqrt(1.-_rho_*_rho_)/(_sigx_*_sigy_*DeuxPi*K);
 return 1./N;
}

double Gdl1RhInt2D::Val_Der(double const xp[], double const* Par
                           , double *DgDpar)
{
 for(int i=0; i<7; i++) DgDpar[i] = 0.;

 double K,W,V,dKdB4;
 if ( 0. < _B4_ && _B4_ < 0.499 ) {
   V = 1.-2.*_B4_;
   W = sqrt(V);
   K = -log( (1.-W)/(1.+W) )/W;
   dKdB4 = ( K-2./(1.-V) )/V;
 } else if ( 0.501 < _B4_ ) {
   V = 1./(2.*_B4_-1.);
   W = sqrt(V);
   K = 2.*W*( Pi/2.-atan(W) );
   dKdB4 = V*( 2.*V/(1.+V) - K );
 } else if ( 0.499 <= _B4_ && _B4_ <= 0.501 ) {
   dKdB4 = -4./3.;
   K = dKdB4 * ( _B4_ - 0.5 ) + 2.;
 } else {
   return(0.);
 }
 double unmr2 = 1.-_rho_*_rho_;
 double N = sqrt(unmr2)/(_sigx_*_sigy_*DeuxPi*K);

 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double z2,PSF,X,Y,XmrY,YmrX,D,dDsD;
 double SPSF=0.;
 {for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   XmrY = X-_rho_*Y;
   YmrX = Y-_rho_*X;
   z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
   D = 1.+z2*(1.+z2*_B4_/2.);
   dDsD = (1.+z2*_B4_)/D;
   /* dg(x,y) */
   PSF = N / D  * w2d[i];
   SPSF += PSF;
   /* dg(x,y)/d(x0) */
   DgDpar[1] += Par[0]* PSF* dDsD*XmrY/_sigx_;
   /* dg(x,y)/d(y0) */
   DgDpar[2] += Par[0]* PSF* dDsD*YmrX/_sigy_;
   /* dg(x,y)/d(sx)*/
   DgDpar[3] += Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
   /* dg(x,y)/d(sy) */
   DgDpar[4] += Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
   /* dg(x,y)/d(rho) */
   DgDpar[5] += Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
   /* dg(x,y)/d(B4) */
   DgDpar[6] += Par[0]* PSF* (-dKdB4/K-z2*z2/2./D);
 }}
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = SPSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[7] = 1.;

 return Par[0] *SPSF + Par[7];
}

void Gdl1RhInt2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = 1.;
}

/////////////////////////////////////////////////////////////////
//++
// Titre	Gdl2Rho2D
// \index{Gdl2Rho2D}
//
//| Cette fonction calcule une gaussienne 2D de hauteur 1 approchee
//| par son developpement limite ordre 3 (see dophot)
//| Meme commentaire que GauRho2D, cf plus haut sauf que:
//|   z**2 = 1/2 (X**2 + Y**2 -2*rho*X*Y)
//|   Z**2 = B2*z2
//|   PSF(x,y) = h / [ 1 + Z**2 + B4**2/2 *Z**4 + B6**2/6 *Z**6 ]
//|   B2,B4,B6 peuvent etre fittes
//| - DL de la gaussienne:  B2=B4=B6=1.
//|  Par [0]=hauteur [1]=x0 [2]=y0 [3]=sigx [4]=sigy [5]=rho
//|      [6]=B4 [7]=B6 [8]=B2 [9]= fond
//|  F(x,y) = Par[0]*PSF(x,y)+Par[9]
//--
/////////////////////////////////////////////////////////////////


//++
Gdl2Rho2D::Gdl2Rho2D()
//
//	Createur
//--
: GeneralPSF2D(10)
{
}

Gdl2Rho2D::~Gdl2Rho2D()
{
}

double Gdl2Rho2D::Value(double const xp[], double const* Par)
{
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2*_rho_*X*Y)/2.;
 double Z2 = _B2_ * _B2_ * z2;
 double Z4 = Z2*Z2;
 double Z6 = Z4*Z2;
 double D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
 return Par[0] /D + Par[9];
}

double Gdl2Rho2D::ValueH(double const xp[], double const* Par)
{
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double z2 = (X*X + Y*Y - 2*_rho_*X*Y)/2.;
 double Z2 = _B2_ * _B2_ * z2;
 double Z4 = Z2*Z2;
 double Z6 = Z4*Z2;
 double D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
 return Par[0] /D + Par[9];
}

double Gdl2Rho2D::Val_Der(double const xp[], double const* Par
                         , double *DgDpar)
{
 double unmr2 = 1.-_rho_*_rho_;
 double X = (xp[0]-_x0_)/_sigx_;
 double Y = (xp[1]-_y0_)/_sigy_;
 double XmrY = X-_rho_*Y;
 double YmrX = Y-_rho_*X;
 double z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
 double Z2 = _B2_ * _B2_ * z2;
 double Z4 = Z2*Z2;
 double Z6 = Z4*Z2;
 double D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
 double dDsDB2 = (1. + _B4_*_B4_*Z2    + _B6_*_B6_/2.*Z4 )/D;
 double dDsD = _B2_*_B2_ * dDsDB2;
 /* g(x,y) */
 double PSF = 1. / D;
 /* dg(x,y)/d(x0) */
 DgDpar[1] = Par[0]* PSF* dDsD*XmrY/_sigx_;
 /* dg(x,y)/d(y0) */
 DgDpar[2] = Par[0]* PSF* dDsD*YmrX/_sigy_;
 /* dg(x,y)/d(sx)*/
 DgDpar[3] = Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
 /* dg(x,y)/d(sy) */
 DgDpar[4] = Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
 /* dg(x,y)/d(rho) */
 DgDpar[5] = Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
 /* dg(x,y)/d(B4) */
 DgDpar[6] = Par[0]* PSF* (-_B4_*Z4/D);
 /* dg(x,y)/d(B6) */
 DgDpar[7] = Par[0]* PSF* (-_B6_*Z6/3./D);
 /* dg(x,y)/d(B2)  */
 DgDpar[8] = Par[0]* PSF* (-2.*_B2_*z2*dDsDB2);
 /* dg(x,y)/d(hauteur) */
 DgDpar[0] = PSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[9] = 1.;

 return Par[0] *PSF + Par[9];
}

void Gdl2Rho2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = parm[7] = parm[8] = 1.;
}

/////////////////////////////////////////////////////////////////
//++
// Titre	Gdl2RhInt2D
// \index{Gdl2RhInt2D}
//
//	fonction integree de Gdl2Rho2d
//--
/////////////////////////////////////////////////////////////////

//++
Gdl2RhInt2D::Gdl2RhInt2D()
//
//	Createur
//--
: GeneralPSF2D(10)
{
}

Gdl2RhInt2D::~Gdl2RhInt2D()
{
}

double Gdl2RhInt2D::Value(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF=0.;
 double X,Y,z2,Z2,Z4,Z6,D;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   Z2 = _B2_ * _B2_ * z2;
   Z4 = Z2*Z2;
   Z6 = Z4*Z2;
   D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
   /* g(x,y) */
   SPSF += w2d[i] / D;
 }
 return Par[0] *SPSF + Par[9];
}

double Gdl2RhInt2D::ValueH(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF=0.;
 double X,Y,z2,Z2,Z4,Z6,D;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   Z2 = _B2_ * _B2_ * z2;
   Z4 = Z2*Z2;
   Z6 = Z4*Z2;
   D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
   /* g(x,y) */
   SPSF += w2d[i] / D;
 }
 return Par[0] *SPSF + Par[9];
}

double Gdl2RhInt2D::Val_Der(double const xp[], double const* Par
                           , double *DgDpar)
{
 for(int i=0; i<=9; i++) DgDpar[i] = 0.;

 double unmr2 = 1.-_rho_*_rho_;
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF=0.;
 double X,Y,XmrY,YmrX,z2,Z2,Z4,Z6,D,dDsD,dDsDB2,PSF;
 {for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   XmrY = X-_rho_*Y;
   YmrX = Y-_rho_*X;
   z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
   Z2 = _B2_ * _B2_ * z2;
   Z4 = Z2*Z2;
   Z6 = Z4*Z2;
   D = 1. + Z2 + _B4_*_B4_/2.*Z4 + _B6_*_B6_/6.*Z6;
   dDsDB2 = (1. + _B4_*_B4_*Z2    + _B6_*_B6_/2.*Z4 )/D;
   dDsD = _B2_*_B2_ * dDsDB2;
   /* dg(x,y) */
   PSF = w2d[i] / D;
   SPSF += PSF;
   /* dg(x,y)/d(x0) */
   DgDpar[1] += Par[0]* PSF* dDsD*XmrY/_sigx_;
   /* dg(x,y)/d(y0) */
   DgDpar[2] += Par[0]* PSF* dDsD*YmrX/_sigy_;
   /* dg(x,y)/d(sx)*/
   DgDpar[3] += Par[0]* PSF* (dDsD*X*XmrY-1.)/_sigx_;
   /* dg(x,y)/d(sy) */
   DgDpar[4] += Par[0]* PSF* (dDsD*Y*YmrX-1.)/_sigy_;
   /* dg(x,y)/d(rho) */
   DgDpar[5] += Par[0]* PSF* (dDsD*X*Y-2.*_rho_/unmr2);
   /* dg(x,y)/d(B4) */
   DgDpar[6] += Par[0]* PSF* (-_B4_*Z4/D);
   /* dg(x,y)/d(B6) */
   DgDpar[7] += Par[0]* PSF* (-_B6_*Z6/3./D);
   /* dg(x,y)/d(B2)  */
   DgDpar[8] += Par[0]* PSF* (-2.*_B2_*z2*dDsDB2);
 }}
 /* dg(x,y)/d(hauteur) */
 DgDpar[0] = SPSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[9] = 1.;

 return Par[0] *SPSF + Par[9];
}

void Gdl2RhInt2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = parm[7] = parm[8] = 1.;
}

/////////////////////////////////////////////////////////////////
//++
// Titre	MofRho2D
// \index{MofRho2D}
//
//| Cette fonction calcule une Moffat 2D
//|   Par [0]=hauteur [1]=x0 [2]=y0 [3]=sigx [4]=sigy [5]=rho
//|       [6]=Gm [7]= fond
//|  PSF(x,y)  = valeur de la Moffat normalisee a un volume = 1
//|  PSF(x,y) = N / [ 1. +  0.5*(X**2 + Y**2 -2*rho*X*Y) ]**Gm
//|           avec X = (x-x0)/sigx et Y = (y-y0)/sigy et Gm>1
//|                N = (1-Gm)*sqrt(1-rho**2)/(2*Pi*sigx*sigy)
//|  le volume de cette Moffat est V=1.
//|  F(x,y) = Par[0]*PSF(x,y)+Par[7]
//--
/////////////////////////////////////////////////////////////////

//++
MofRho2D::MofRho2D()
//
//	Createur 
//--
: GeneralPSF2D(8)
{
}

MofRho2D::~MofRho2D()
{
}

double MofRho2D::Value(double const xp[], double const* Par)
{
 double N = (_Gm_-1.)*sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double z2 = (X*X + Y*Y -2.*_rho_*X*Y)/2.;
 z2 = _Gm_*log(1. + z2);
 if( z2<MINEXPM ) return Par[0] *N*EXPO(-z2) + Par[7];
    else return Par[7];
}

double MofRho2D::ValueH(double const xp[], double const* Par)
{
 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double z2 = (X*X + Y*Y -2.*_rho_*X*Y)/2.;
 z2 = _Gm_*log(1. + z2);
 if( z2<MINEXPM ) return Par[0] *EXPO(-z2) + Par[7];
    else return Par[7];
}

double MofRho2D::VolPSF(double const* Par)
{
 double N = (_Gm_-1.)*sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
 return 1./N;
}

double MofRho2D::Val_Der(double const xp[], double const* Par
                        , double *DgDpar)
{
 double unmr2 = 1.-_rho_*_rho_;
 double N = (_Gm_-1.)*sqrt(unmr2)/DeuxPi/_sigx_/_sigy_;
 double X = (xp[0] - _x0_)/_sigx_;
 double Y = (xp[1] - _y0_)/_sigy_;
 double XmrY = X-_rho_*Y;
 double YmrX = Y-_rho_*X;
 double z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
 double D = 1. + z2;
 double lD = log(D);
 /* g(x,y) */
 double PSF = _Gm_*lD;
 if( PSF<MINEXPM ) PSF = N * EXPO(-PSF); else PSF = 0.;
 /* dg(x,y)/d(x0) */
 DgDpar[1] = Par[0]* PSF* XmrY/_sigx_ * _Gm_/D;
 /* dg(x,y)/d(y0) */
 DgDpar[2] = Par[0]* PSF* YmrX/_sigy_ * _Gm_/D;
 /* dg(x,y)/d(sx)*/
 DgDpar[3] = Par[0]* PSF* (X*XmrY*_Gm_/D - 1.)/_sigx_;
 /* dg(x,y)/d(sy) */
 DgDpar[4] = Par[0]* PSF* (Y*YmrX*_Gm_/D - 1.)/_sigy_;
 /* dg(x,y)/d(rho) */
 DgDpar[5] = Par[0]* PSF* (X*Y*_Gm_/D - 2.*_rho_/unmr2);
 /* dg(x,y)/d(Gm) */
 DgDpar[6] = Par[0]* PSF* (1./(_Gm_-1.) - lD);
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = PSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[7] = 1.;

 return Par[0] *PSF + Par[7];
}

void MofRho2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = 3.;
}

/////////////////////////////////////////////////////////////////
//++
// Titre	MofRhInt2D
// \index{MofRhInt2D}
//
//	fonction integree de MofRho2d
//--
/////////////////////////////////////////////////////////////////

//++
MofRhInt2D::MofRhInt2D()
//
//	Createur 
//--
: GeneralPSF2D(8)
{
}

MofRhInt2D::~MofRhInt2D()
{
}

double MofRhInt2D::Value(double const xp[], double const* Par)
{
 double N = (_Gm_-1.)*sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF = 0.;
 double z2,X,Y;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   /* g(x,y) */
   z2  = _Gm_*log(1. + z2);
   if( z2<MINEXPM ) SPSF += EXPO(-z2) * w2d[i];
 }
 return Par[0] * N*SPSF + Par[7];
}

double MofRhInt2D::ValueH(double const xp[], double const* Par)
{
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF = 0.;
 double z2,X,Y;
 for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   z2 = (X*X + Y*Y - 2.*_rho_*X*Y)/2.;
   /* g(x,y) */
   z2  = _Gm_*log(1. + z2);
   if( z2<MINEXPM ) SPSF += EXPO(-z2) * w2d[i];
 }
 return Par[0] *SPSF + Par[7];
}

double MofRhInt2D::VolPSF(double const* Par)
{
 double N = (_Gm_-1.)*sqrt(1.-_rho_*_rho_)/(DeuxPi*_sigx_*_sigy_);
 return 1./N;
}

double MofRhInt2D::Val_Der(double const xp[], double const* Par
                          , double *DgDpar)
{
 for(int i=0; i<=7; i++) DgDpar[i] = 0.;

 double unmr2 = 1.-_rho_*_rho_;
 double N = (_Gm_-1.)*sqrt(unmr2)/(DeuxPi*_sigx_*_sigy_);
 double x = xp[0] - _x0_;
 double y = xp[1] - _y0_;
 double SPSF = 0.;
 double X,Y,XmrY,YmrX,z2,D,lD,PSF;
 {for(int i=0; i<nd2d; i++) {
   X = (x+dx2d[i])/_sigx_;
   Y = (y+dy2d[i])/_sigy_;
   XmrY = X-_rho_*Y;
   YmrX = Y-_rho_*X;
   z2 = (X*(XmrY-_rho_*Y)+Y*Y)/2.;
   D = 1. + z2;
   lD = log(D);
   /* g(x,y) */
   PSF = _Gm_*lD;
   if( PSF<MINEXPM ) PSF = N * EXPO(-PSF) * w2d[i]; else PSF = 0.;
   SPSF += PSF;
   /* dg(x,y)/d(x0) */
   DgDpar[1] += Par[0]* PSF* XmrY/_sigx_ * _Gm_/D;
   /* dg(x,y)/d(y0) */
   DgDpar[2] += Par[0]* PSF* YmrX/_sigy_ * _Gm_/D;
   /* dg(x,y)/d(sx)*/
   DgDpar[3] += Par[0]* PSF* (X*XmrY*_Gm_/D - 1.)/_sigx_;
   /* dg(x,y)/d(sy) */
   DgDpar[4] += Par[0]* PSF* (Y*YmrX*_Gm_/D - 1.)/_sigy_;
   /* dg(x,y)/d(rho) */
   DgDpar[5] += Par[0]* PSF* (X*Y*_Gm_/D - 2.*_rho_/unmr2);
   /* dg(x,y)/d(Gm) */
   DgDpar[6] += Par[0]* PSF* (1./(_Gm_-1.) - lD);
 }}
 /* dg(x,y)/d(Vol) */
 DgDpar[0] = SPSF;
 /* dg(x,y)/d(Fond) */
 DgDpar[7] = 1.;

 return Par[0] *SPSF + Par[7];
}

void MofRhInt2D::DefaultParam(double *parm)
{
  for (int i=0; i<mNPar; i++) parm[i] = 0.;
  parm[3] = parm[4] = 1.;  // Sigx Sigy
  parm[6] = 3.;
}


#ifdef IS_IT_USEFUL
#undef _sigx_
#undef _sigy_
#undef _rho_
#undef _x0_
#undef _y0_
#undef _Gm_
#undef _B4_
#undef _B6_
#undef _B2_

//==============================================================================
// CLASSES DE FONCTIONS 2D type Xi2 AVEC PARAMETRES POUR LE FIT pixel taille 1x1
// la taille du pixel est importante quand on utilise les PSF integrees
//    (x,y x0,y0 sigmaX.... sont en unites de pixels !!!)
//==============================================================================

/////////////////////////////////////////////////////////////////
//++
// Titre	X2-GauRho2D
// \index{X2-GauRho2D}
//
//	Chi2 pour une Gaussienne+fond 2D (voir detail dans GauRho2D).
//--
/////////////////////////////////////////////////////////////////

//++
// X2-GauRho2D::X2-GauRho2D()
//	Createur 
//--
X2_GauRho2D::X2_GauRho2D()
: GeneralXi2(7)
{
  gaurho2d = new GauRho2D;
}

X2_GauRho2D::~X2_GauRho2D()
{
  delete gaurho2d;
}

double X2_GauRho2D::Value(GeneralFitData& data, double* parm, int& ndataused)
{
  // DBASSERT( data.NVar()==2 );
 double x[2],z;

 double c2 = 0.;
 ndataused = 0;
 for(int k=0;k<data.NData();k++) {
   if( ! data.IsValid(k) ) continue;
   x[0] = data.X(k); x[1] = data.Y(k);
   z = (data.Val(k)-gaurho2d->Value(x,parm))/data.EVal(k);
   c2 += z*z;
   ndataused++;
 }
 return c2;
}

double X2_GauRho2D::Derivee2(GeneralFitData& data, int i,int j, double* parm)
{
  // DBASSERT( data.NVar()==2 && i<7 && j<7);
 double x[2],dparm[7];

 double d2c2 = 0.;
 for(int k=0;k<data.NData();k++) {
   if( ! data.IsValid(k) ) continue;
   x[0] = data.X(k); x[1] = data.Y(k);
   gaurho2d->Val_Der(x,parm,dparm);
   d2c2 += 2.*dparm[i]*dparm[j]/(data.EVal(k)*data.EVal(k));
 }
 return d2c2;
}
#endif /* IS_IT_USEFUL */
