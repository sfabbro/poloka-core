#include <fstream>
#include <iomanip>
#include "fct2dfit.h"
//#include "fct1dfit.h"
//#include "perrors.h"
#include "psf.h"
#define SIMPSON4

// Converter :
BaseStarList* PSF2Base(PSFStarList * This)
{ return (BaseStarList*) This;}
 
const BaseStarList* PSF2Base(const PSFStarList * This)
{ return (BaseStarList*) This;}  


//#include "nbinteg.h"

//#include "fmath.h" // for floorf


//++
// Class	FuncPtrGFF
//
//	Une GeneralPSF2D fabriquée à partir d'un pointeur sur fonction.
//	Le type du pointeur est "FuncPtrGFF::Func" qui correspond à
//|	typedef double (*Func)(double dx, double dy, double const* par)
//
//	Le nombre de paramètres indiqués ne tient pas compte de x et y.
//	En général, le premier paramètre sera le volume, et le suivant le
//	fond.
//
//	Le vecteur de paramètres passé à la FuncPtrGFF sera {x0, y0, autres...},
//	et elle le retransmet au pointeur de fonction en séparant x0 et y0 des
//	autres paramètres.
//--

//++
// Links	Parents
// GeneralPSF2D
// GeneralFunction
//--

//++
// Titre	Constructeurs
//--

//	Crée un FuncPtrGFF à partir d'un pointeur de fonction.
//	nPar est le nombre de paramètres en plus de x et y.
FuncPtrGFF::FuncPtrGFF(Func f, int nPar)
  :GeneralPSF2D(nPar+2), mFunc(f), mNPar(nPar), 
 mTmpParm(new double[nPar])
{}

FuncPtrGFF::~FuncPtrGFF()
{
  delete[] mTmpParm;
}

//++
// Titre	Méthodes
//--

//++
double
FuncPtrGFF::Value(double const xp[], double const* parm)
//
//	Appelle f(xp[0] - parm[1], xp[1] - parm[2], parm[0], parm[3]...)
//--
{
  memcpy(mTmpParm+1, parm+3, (mNPar-1) * sizeof(double));
  mTmpParm[0] = parm[0];
  return mFunc(xp[0] - parm[1], xp[1] - parm[2], mTmpParm);
}

//********************   DEFINITION PSF   *********************


//	Une Point Spread Function. Elle est toujours associée à une
//	GeneralPSF2D, qui peut être spécifiée soit via une FuncPtrGFF,
//	soit dans une liste de PSF standards.
//

//Links	Voir aussi
// FuncPtrGFF
// GeneralPSF2D
// GeneralFunction
//--
// GauRho2D
// GauRhInt2D
// GdlRho2D
// GdlRhInt2D
// Gdl1Rho2D
// Gdl1RhInt2D
// Gdl2Rho2D
// Gdl2RhInt2D
// MofRho2D
// MofRhInt2D
////--


//Les différents types sont:
//	* psfGaussienne
//	* psfGaussInteg : gaussienne intégrée (sur la surface du pixel)
//	* psfDelGaussienne : gaussienne, DL 1er ordre
//	* psfDelGaussinteg : gaussienne intégrée, DL 1er ordre
//	* psfDel1Gaussienne : gaussienne, DL 1er ordre
//	* psfDel1Gaussinteg : gaussienne intégrée, DL 1er ordre
//	* psfMoffat : Moffat
//	* psfMoffInteg : Moffat intégrée
//	* psfDel2Gaussienne : gaussienne, DL 2e ordre
//	* psfDel2Gaussinteg : gaussienne intégrée, DL 2e ordre
//
//	Voir les classes correspondantes pour la description des
//	paramètres.
//
//	ns est utilisé pour la fabrication d'images synthétiques : 
// nombre de sigmas à partir duquel on néglige la contribution.




//utilise  mTypePSF pour creer une fonction Gaussienne , integree etc.
// MakeStdGFF() gere l'allocation unique des differents type
PSF::PSF(int psfKind, double ns)
: mTypePSF(psfKind), mNParm(0), mNSig0(ns)
{
  MakeStdGFF();
}

//  Constructeur, a partir d'un pointeur sur une classe GeneralPSF2D.
PSF::PSF(GeneralPSF2D *psf2d, double ns)
: mTypePSF(psfCustom), mNParm(0), mNSig0(ns)
{
  mGFF = psf2d;
  mNParm = mGFF->NPar();
}

//  Constructeur, à partir d'un pointeur sur fonction. Voir FuncPtrGFF.
PSF::PSF(FuncPtrGFF::Func f, int nPar, double ns)
: mTypePSF(psfCustom), mNParm(nPar+2), mNSig0(ns)
{
   mGFF = new FuncPtrGFF(f, nPar);
}

// n'alloue pas la PSF, pour pouvoir lire l'etoile
// il faut allouer ensuite avec SetPSF(PSF *)
PSF::PSF()
: mTypePSF(psfCustom), mNParm(0), mNSig0(3)
{}

PSF::~PSF()
{
  delete mGFF;
}

// 	Méthodes




//	Retourne la valeur de la PSF en (x,y).
double
PSF::Value(double xx, double yy, double const* par) const
{
  double xp[2] = {xx,yy};
  return mGFF->Value(xp, par);
}




// allocation de la PSF avec le type connu


GeneralPSF2D* PSF::_instance[10]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,
			      NULL,NULL,NULL};


void
PSF::MakeStdGFF()
{
  switch (mTypePSF) {
    case psfGaussienne :
      if (_instance[psfGaussienne] == NULL)
	_instance[psfGaussienne]  = new GauRho2D;
      break;
    case psfGaussInteg :
      if (_instance[psfGaussInteg] == NULL)
	_instance[psfGaussInteg]  = new GauRhInt2D;
      break;
    case psfDelGaussienne :
      if (_instance[psfDelGaussienne] == NULL)
	_instance[psfDelGaussienne]  = new GdlRho2D;
      break;
    case psfDelGaussinteg :
      if (_instance[psfDelGaussinteg] == NULL)
	_instance[psfDelGaussinteg]  = new GdlRhInt2D;
      break;
    case psfDel1Gaussienne :
      if (_instance[psfDel1Gaussienne] == NULL)
	_instance[psfDel1Gaussienne]  = new Gdl1Rho2D;
      break;
    case psfDel1GaussInteg :
      if (_instance[psfDel1GaussInteg] == NULL)
	_instance[psfDel1GaussInteg]  = new Gdl1RhInt2D;
      break;
    case psfDel2Gaussienne :
      if (_instance[psfDel2Gaussienne] == NULL)
	_instance[psfDel2Gaussienne]  = new Gdl2Rho2D;
      break;
    case psfDel2GaussInteg :
      if (_instance[psfDel2GaussInteg] == NULL)
	_instance[psfDel2GaussInteg]  = new Gdl2RhInt2D;
      break;
    case psfMoffat :
      if (_instance[psfMoffat] == NULL)
	_instance[psfMoffat]  = new MofRho2D;
      break;
    case psfMoffInteg :
      if (_instance[psfMoffInteg] == NULL)
	_instance[psfMoffInteg]  = new MofRhInt2D;
      break;
    default:
	{
	    cerr << "parmErr" << endl ;
	    exit(0);
	    //THROW(parmErr);
	}
  }
  mGFF = _instance[mTypePSF];
  mNParm = mGFF->NPar();
}









//********************   FIN DEFINITION PSF   *********************


//********************   DEFINITION PSFStar   *********************



//++
// Class	PSFStar
// include	psf.h
//
//	Une étoile avec une PSF, essentiellement pour faire de la simulation
//	d'images. Pour fitter une PSF, utiliser PSFFitStar.
//
// les createurs mettent les parametres aux valeurs 
// par defaut (cf fitfct2d.cc)  sauf pour la creation a partir d'une 
// PSFStar (copie)

// attention la variable flux de SEStar est sur-ecrite 

// CREATEURS


PSFStar::PSFStar(PSF* psf)
: SEStar()
{
  if ( psf == NULL )
    psf = new PSF(PSF::psfGaussInteg );
  mPSF=psf;
  mParm = new double[mPSF->NParm()];
  fond_fit = 0. ;
  mPSF->GetGFF()->DefaultParam(mParm);
  for(int j = 0; j< mPSF->NParm() ; j++)
       mParm[j] = 0 ;
}

PSFStar::PSFStar(SEStar const & sestar, PSF* psf)
: SEStar(sestar)
{
  if ( psf == NULL )
    psf = new PSF(PSF::psfGaussInteg );
  mPSF=psf;
  mParm = new double[mPSF->NParm()];
  fond_fit = 0. ;
  // mise des parametres aux valeurs par defaut
  mPSF->GetGFF()->DefaultParam(mParm);
}

PSFStar::PSFStar(PSFStar const & star)
  : SEStar(star)
{
  mPSF= star.GetPSF();
  mParm = new double[mPSF->NParm()];
  fond_fit = star.Fond_fit();
  for(int i = 0 ; i < mPSF->NParm() ; i++)
    mParm[i] =  star.Parm(i);
}


PSFStar::PSFStar(double xx, double yy, double ff,  PSF* aPsf)
: SEStar(xx,yy,ff)
{ 
  if ( aPsf == NULL )
    aPsf = new PSF(PSF::psfGaussInteg );     
  mPSF=aPsf;
  mParm = new double[mPSF->NParm()];
  fond_fit = 0. ;
  // mise des parametres aux valeurs par defaut
  mPSF->GetGFF()->DefaultParam(mParm);
}

PSFStar::~PSFStar()
{
  if(mParm) delete[] mParm;
}









void
PSFStar::SetPSF(PSF* psf)
{
  // copie pointeur de PSF
  mPSF = psf;
  // reallocation des variables de PSFStar mParm
  delete[] mParm;
  mParm = new double[mPSF->NParm()];
  // mise des parametres aux valeurs par defaut
  mPSF->GetGFF()->DefaultParam(mParm);
}

void
PSFStar::UpdParmVect() const
{
  mParm[0]   = flux;
  mParm[1]   = x;
  mParm[2]   = y;
  mParm[mPSF->NParm()-1] = fond_fit;
}

void
PSFStar::UpdFromParm() const
{
  const_cast<double&>(flux) = mParm[0];
  const_cast<double&>(x) = mParm[1];
  const_cast<double&>(y) = mParm[2];
  const_cast<double&>(fond_fit) = mParm[mPSF->NParm()-1];
}

// read & write with or without end-of-line
void
PSFStar::dumpn(ostream& s) const
{
  SEStar::dumpn(s);
  s << " type PSF: " << (GetPSF())->Kind();
  s << " fonf fitte: " << Fond_fit() ;
  s << " sigma X: " << SigX() ;
  s << " sigma Y: " <<  SigY() ;
  s << "  rho: " <<  Rho() ;
}

void
PSFStar::dump(ostream& s)  const
{
  dumpn(s); s << endl ;
}

void
PSFStar::writen(ostream& s) const
{
  SEStar::writen(s); 
  s  << resetiosflags(ios::scientific) ;
  s  << setiosflags(ios::fixed) ;
  s  << setprecision(4) ;
  s << (GetPSF())->Kind()<< " " ;
  s << Fond_fit()<< " " ;
  s << SigX() << " " ;
  s <<  SigY() << " " ;
  s <<  Rho() << " " ;
}

void
PSFStar::write(ostream& s) const
{
  writen(s); s << endl ;
}




// pour garder la possibilite d'emboitement
// des Read
void
PSFStar::Read(istream& r, const char *Format, int & type_psf )
{
  SEStar::read_it(r, Format);
  r >> type_psf ;
  r >> Fond_fit() ;
  r >> SigX()  ;
  r >>  SigY()  ;
  r >>  Rho()  ;
}


// creation avec PSF donnee, quelque soit le type lu
PSFStar*  PSFStar::read(istream& r, const char *Format, PSF* mypsf)
{
  PSFStar *pstar = new PSFStar(mypsf);
  int type_lu = 0 ;
  pstar->Read(r, Format, type_lu);  
  int type_psf = mypsf->Kind();
  if ( type_lu  != type_psf)
    {
      cout << "type psf impose a la lecture incompatible avec le type lu " 
<< endl ;
      //THROW non fatal
    }
  return(pstar);
}



PSFStar*  PSFStar::read(istream& r, const char *Format)
{ 
  int type_psf = 0 ;
  SEStar *pse = SEStar::read(r, Format);
  r >> type_psf ;
  PSF * mypsf = new PSF(type_psf);
  PSFStar * pstar = new PSFStar(*pse, mypsf);
  delete pse ;
  
  r >> pstar->Fond_fit() ;
  r >> pstar->SigX()  ;
  r >>  pstar->SigY()  ;
  r >>  pstar->Rho()  ;

  return(pstar);

}


std::string
PSFStar::WriteHeader_(ostream & pr, const char *i ) const
{
  if (i==NULL) i = "";
  std::string SEStarFormat = SEStar::WriteHeader_(pr,i);
  pr    << "# tpsf"<< i <<" : type PSF" << endl 
        << "# fd_fit"<< i <<" : background" << endl 
        << "# sigx"<< i <<" : sigma x in pixels" << endl 
        << "# sigy"<< i <<" : sigma y in pixels" << endl 
        << "# rho"<< i <<"  : " << endl ;
  static char format[256];
  sprintf(format,"%s PSFStar %d",SEStarFormat.c_str(), 0);
  return format;return(format);
}











//++
double
PSFStar::operator() (double xx, double yy) const
//
//	Valeur de la PSF en ce point.
//--
{
  // Synchronisation du vecteur de parametres.
  // La valeur retournee (ne) tient (pas) compte du fond de l'etoile.
  // Est-ce le bon choix ?
  UpdParmVect();
  return mPSF->Value(xx,yy,mParm);
}

//++
// void PSFStar::ImgAdd(Image & img) const
//	Ajoute/soustrait l'étoile sur l'image.
//--


// utilise a la fois x et y pour la position, et operateur (xx,yy)
// qui fait UpdParmVect.
// la star est en coordonnees XY 

void  PSFStar::ImgAdd(Image& img) const
{
  double r = mPSF->mNSig0 * Sig();
  int imin = (int) floor(x - r);
  int jmin = (int) floor(y - r);
  int imax = (int)  ceil(x + r);
  int jmax = (int)  ceil(y + r);

  int XOrg = 0 ;
  int YOrg = 0 ;

  if (imin - XOrg < 0) imin = XOrg;
  if (jmin - YOrg < 0) jmin = YOrg;
  if (imax - XOrg > img.Nx()) imax = img.Nx() + XOrg;
  if (jmax - YOrg > img.Ny()) jmax = img.Ny() + YOrg;

  for (int i=imin; i<imax; i++)
    for (int j=jmin; j<jmax; j++)
      {
	img(i - XOrg,j - YOrg) += (*this)(i+DECALAGE_IJ_XY, j+DECALAGE_IJ_XY) ;
      }
}

void  PSFStar::ImgAdd_2(Image& img) const
{
  double r = mPSF->mNSig0 * Sig();
  int imin = (int) floor(x - r);
  int jmin = (int) floor(y - r);
  int imax = (int)  ceil(x + r);
  int jmax = (int)  ceil(y + r);

  int XOrg = 0 ;
  int YOrg = 0 ;

  if (imin - XOrg < 0) imin = XOrg;
  if (jmin - YOrg < 0) jmin = YOrg;
  if (imax - XOrg > img.Nx()) imax = img.Nx() + XOrg;
  if (jmax - YOrg > img.Ny()) jmax = img.Ny() + YOrg;

  for (int i=imin; i<imax; i++)
    for (int j=jmin; j<jmax; j++)
      {
	img(i - XOrg,j - YOrg) += (*this)(i, j) ;
      }
}

Image& operator += (Image& img, PSFStar const& psf)
{
    psf.ImgAdd(img);
    return img;
}





void  PSFStar::ImgSub(Image& img) const
{
  double r = mPSF->mNSig0 * Sig();
  int imin = (int) floor(x - r);
  int jmin = (int) floor(y - r);
  int imax = (int)  ceil(x + r);
  int jmax = (int)  ceil(y + r);

  int XOrg = 0 ;
  int YOrg = 0 ;

  if (imin - XOrg < 0) imin = XOrg;
  if (jmin - YOrg < 0) jmin = YOrg;
  if (imax - XOrg > img.Nx()) imax = img.Nx() + XOrg;
  if (jmax - YOrg > img.Ny()) jmax = img.Ny() + YOrg;

  for (int i=imin; i<imax; i++)
    for (int j=jmin; j<jmax; j++)
      img(i - XOrg,j - YOrg) -= (*this)(i+DECALAGE_IJ_XY, j+DECALAGE_IJ_XY);
}

Image& operator -= (Image& img, PSFStar const& psf)
{
    psf.ImgSub(img);
    return img;
}

//********************   FINDEFINITION PSFStar   *********************


#ifdef IS_IT_USEFUL

//********************   DEFINITION PSFFitStar   *********************



//++
// Class	PSFFitStar
// Lib		
// include	psf.h
//
//	Une étoile avec une PSF, pour faire des fits.
//--

//++
// Links	Parents
// PSFStar
// SEStar
//--

//++
// Links	Voir aussi
// PSF
//--

//++
// double& Err(int i)        
//	L'erreur sur le paramètre i. Version const.
// double& ErrXXX()
//	avec XXX = (Flux, X, Y, Fond, SigX, SigY, Rho) : erreur sur le
//	paramètre.
// double& Xi2()
//	Le Xi2 du fit
//--



PSFFitStar::PSFFitStar(PSF* psf)
: PSFStar(psf)
{
  mErr = new double[mPSF->NParm()];
  for(int j = 0; j< mPSF->NParm() ; j++)
       mErr[j] = 0 ;
  mXi2 = -1.;
}

PSFFitStar::PSFFitStar(PSFFitStar const & star)
: PSFStar(star)
{
  mErr = new double[mPSF->NParm()];
  for(int j = 0; j< mPSF->NParm() ; j++)
       mErr[j] = star.Err(j) ;
  mXi2 = star.Xi2();
}

PSFFitStar::PSFFitStar(PSFStar const & star)
: PSFStar(star)
{
  mErr = new double[mPSF->NParm()];
  for(int j = 0; j< mPSF->NParm() ; j++)
       mErr[j] = 0 ;
  mXi2 = -1.;
}

PSFFitStar::PSFFitStar(SEStar const & star, PSF *psf)
: PSFStar(star,psf)
{
  mErr = new double[mPSF->NParm()];
  for(int j = 0; j< mPSF->NParm() ; j++)
       mErr[j] = 0 ;
  mXi2 = -1.;
}



PSFFitStar::~PSFFitStar()
{
   if(mErr) delete[] mErr;
}




// read & write with or without end-of-line
void
PSFFitStar::dumpn(ostream& s)  const
{
  PSFStar::dumpn(s);
  s << " Erreur fit position X: " << ErrX() ;
  s << " Erreur fit position Y: " << ErrY() ;
  s << " Erreur fit : " << ErrFlux() ;
  s << " Erreur fit fond : " <<ErrFond() ;
  s << " Erreur fit sigma X : "<< ErrSigX()  ;
  s << " Erreur fit sigma Y : "<< ErrSigY()  ;
  s << " Erreur fit rho: " << ErrRho()  ;
  s << " Xi2 fit: " << Xi2()  ;
}

void
PSFFitStar::dump(ostream& s)  const
{
  dumpn(s); s << endl ;
}
void
PSFFitStar::writen(ostream& s) const
{
  PSFStar::writen(s);
  s  << resetiosflags(ios::scientific) ;
  s  << setiosflags(ios::fixed) ;
  s  << setprecision(4) ;
  s << ErrX() << " " ;
  s << ErrY() << " " ;
  s << ErrFlux() << " " ;
  s << ErrFond() << " " ;
  s << ErrSigX()  << " " ;
  s << ErrSigY()  << " " ;
  s << ErrRho()  << " " ;
  s << Xi2()  << " " ;
}

void
PSFFitStar::write(ostream& s) const
{
  writen(s); s << endl ;
}



void
PSFFitStar::Read(istream& r, const char *Format, int & type_psf)
{
  PSFStar::Read(r, Format, type_psf);
  r >> ErrX() ;
  r >> ErrY() ;
  r >> ErrFlux() ;
  r >> ErrFond() ;
  r >> ErrSigX()  ;
  r >> ErrSigY()  ;
  r >> ErrRho()  ;
  r >> Xi2()  ;
}

PSFFitStar*  PSFFitStar::read(istream& r , const char *Format, PSF* mypsf )
{
  PSFFitStar *pstar = new PSFFitStar(mypsf);
  int type_lu = 0 ;
  pstar->Read(r, Format, type_lu); 
  int type_psf = mypsf->Kind();
  if ( type_lu  != type_psf)
    {
      cout << "type psf impose a la lecture incompatible avec le type lu " 
	   << endl ;
      //THROW non fatal
    }
  return(pstar);
}


PSFFitStar*  PSFFitStar::read(istream& r, const char *Format)
{
 PSFStar * pst =   PSFStar::read(r, Format);
 PSFFitStar * pstar = new PSFFitStar(*pst);
 delete pst ;
 r >> pstar->ErrX() ;
 r >> pstar->ErrY() ;
 r >> pstar->ErrFlux() ;
 r >> pstar->ErrFond() ;
 r >> pstar->ErrSigX()  ;
 r >> pstar->ErrSigY()  ;
 r >> pstar->ErrRho()  ;
 r >> pstar->Xi2()  ;
 return(pstar);
}

const char *  
PSFFitStar::WriteHeader_(ostream & pr, const char *i) const
{
  if (i==NULL) i = "";
  const char * format = PSFStar::WriteHeader_(pr,i);
  pr    << "# errx"<< i <<" : error on x fitted position " << endl
	<< "# erry"<< i <<" : error on y fitted position " << endl
	<< "# errflux"<< i <<" : error on fitted flux " << endl
	<< "# errfd"<< i <<" :  error on fitted background " << endl
	<< "# errsx"<< i <<" : error on fitted sigma x" << endl
	<< "# errs"<< i <<"y :  error on fitted sigma y" << endl
	<< "# errrho"<< i <<" : error on fitted rho " << endl
	<< "# xi2"<< i <<" : fit xi2" << endl ;
  return(format);
}







//++
int PSFFitStar::Fit(Image const& img, Image const& noiseimg, int lp)
//
//	Ajustement de l'étoile sur l'image.
//--
{
 int dx = img.Nx();
 int dy = img.Ny();
 int npar = mPSF->NParm();

 GeneralFit       psffit( mPSF->GetGFF() );

 GeneralFitData   psfdata(2,dx*dy+10);
 float xx[2];
int XOrg = 0, YOrg = 0 ;

// DH : la star est en XY, on mets les data (ie l'image) en XY

 for(int i=0;i<dx;i++) for(int j=0;j<dy;j++) {
   // le coin inferieur gauche du pixel [0,0] a la coordonnee [XOrg,YOrg]
   // DH : ???? [XOrg+0.5,YOrg+0.5] plutot ?
   xx[0] = (double) (XOrg + i) + DECALAGE_IJ_XY;
   xx[1] = (double) (YOrg + j) + DECALAGE_IJ_XY;
   psfdata.AddData(xx,img(i,j),noiseimg(i,j));
   // le contenu du pixel (0,0) est associe a (x=0.5, y=0.5)
 }

 if(lp>1) psffit.SetDebug(1);
 psffit.SetData(&psfdata);
 // initialisation des bornes des valeurs autorisees 
 // pour les parametres fittes avec le contenu de mParm
 psffit.SetParam(0,"Flux",mParm[0],fabs(mParm[0])/100.);
 psffit.SetParam(1,"X0",mParm[1],0.1
                ,mParm[1]-(double) dx/2.,mParm[1]+(double) dx/2.);
 psffit.SetParam(2,"Y0",mParm[2],0.1
                ,mParm[2]-(double) dy/2.,mParm[2]+(double) dy/2.);
 psffit.SetParam(3,"Sx",mParm[3],0.05,0.1,(double)dx);
 psffit.SetParam(4,"Sy",mParm[4],0.05,0.1,(double)dy);
 psffit.SetParam(5,"Rho",mParm[5],0.01,-1.,1.);
 psffit.SetParam(npar-1,"Fond",mParm[npar-1],1.);
 // Si pas fixe la strategie de fit doit etre plus subtile
 if(npar>7) for(int i=6;i<npar-1;i++) psffit.SetFix(i,mParm[i]);

 int rc = psffit.Fit();

 if(rc>=0) {
   rc = 0;
   Xi2() = psffit.GetChi2();
   for(int i=0;i<npar;i++) { mParm[i] = psffit.GetParm(i);
                             mErr[i]  = psffit.GetParmErr(i); }
   // remplace flux, x, y fd par les valeurs fittees
   UpdFromParm();
   double h = mParm[0] / mPSF->GetGFF()->VolPSF(mParm) + mParm[npar-1];
   if(lp>1) cout << "hauteur= " << h << endl;
 } else for(int i=0;i<npar;i++) mParm[i] = mErr[i]  = 0.;
 if(lp>1)
   {
     printf("Resultat Fit: Parametres ");
     for(int i=0;i<npar;i++) printf(" %d= %g,",i,mParm[i]);
     printf("\n");
   }
 return rc;
}


int PSFFitStar::IniFitSE(bool fond_force, double fond, int lp)
{
  int npar = mPSF->NParm();
  mParm[0] = flux; // flux mesure par SExtractor
  mParm[1] = x +  DECALAGE_SE_XY ;
  mParm[2] = y +  DECALAGE_SE_XY ;
  double ux= Mxx();
  if ( ux > 1.e-10)
    SigX() = sqrt(ux) ;
  else
    SigX() = 0. ;
  double uy = Myy();
  if ( uy > 1.e-10)
    SigY() = sqrt(uy) ;
  else
    SigY() = 0. ;
  if ( (ux > 1.e-10)&& ( uy > 1.e-10))
    {
      double ur = Mxy()/(ux*uy) ;
      Rho() = ur;
    }
  else
    Rho() = 0. ;
  // si rho en dehors de [-1;1] ie abs(rho)-1 > 0
  if ( (fabs(Rho())-1.) > 1.e-10 )
    Rho() = 0. ;
  double fond_star =  Fond(); // fond mesure par SExtractor (!= fond_fit)
  if (fond_force)
    fond_star = fond ;
  mParm[npar-1] = fond_star ;

  if(lp>0) {
    cout << "IniFitSE (SE parameters): " 
	 << x << " "   
	 << y << " "  
	 << Mxx() << " " 
	 << Myy() << " " 
	 << Mxy() << " "  
	 << A() << " "   
	 << B() << " " 
	 << Gyr_Angle() << " " << endl ;
    cout << "IniFitSE : " << SigX() << " " << SigY() << " " 
	 << Rho() << endl ;
    printf("IniFitSE: Parametres ");
    for(int i=0;i<npar;i++) printf(" %d= %g,",i,mParm[i]);
    printf("\n");
  }
 
  return 0;
}
 

//++
int PSFFitStar::IniFit(Image const& img, Image const& noiseimg, int lp)
//
//	Initialisation à des valeurs raisonnables pour appeler Fit.
//--
{
  int XOrg = 0, YOrg = 0 ;

  // on ne prend pas les defauts Peida
  //mPSF->GetGFF()->DefaultParam(mParm);


 int npar = mPSF->NParm();
 int x0[2] = { XOrg , YOrg  };
 int dx[2] = { img.Nx(), img.Ny() };
 if(lp>1) printf("PSFFitStar::IniFit  x0,y0=%d,%d  dx,dy=%d,%d  npar=%d\n"
                ,x0[0],x0[1],dx[0],dx[1],npar);
 // test sur taille image
 if(dx[0]<3 || dx[1]<3) return 1;

 // recherche du volume et du fond
 int nfond = 0, nvol = 0;
 double vol = 0., ffond = 0.;
 double* fondval = new double[2*(dx[0]+dx[1]-2)];
 
 
 for(int i=0;i<dx[0];i++) for(int j=0;j<dx[1];j++)
   if( noiseimg(i,j) > 0. ) 
     {
       nvol++; vol += img(i,j);
       // fond sur le bord du pave
       if( i==0 || i==dx[0]-1 || j==0 || j==dx[1]-1 ) 
	 {
	   fondval[nfond] = img(i,j);
	   nfond++;
	   ffond += img(i,j);
	 }
     }
 int d = (dx[0]+dx[1])/2;
 if( nfond > 2*d ) {
   HeapSort(nfond,fondval);
   ffond = 0.;
   nfond -= d;
   for(int i=0;i<nfond;i++) ffond += fondval[i];
 }
 delete[] fondval;
 if(nfond==0) 
   {
     cerr << "Echec fond " << endl ;
     return 2;
   }
 
 ffond /= (double) nfond;
 vol -= ffond * (double) nvol;


 if(lp>1)
   printf("      vol=%g nvol=%d   fond=%g nfond=%d\n",vol,nvol,ffond,nfond);

 
 // initialisation
 
 mParm[0] = vol;
 mParm[1] = x;
 mParm[2] = y;
 mParm[npar-1] = ffond;

 if(lp>0) {
   printf("IniFit: Parametre");
   for(int i=0;i<npar;i++) printf(" %d= %g,",i,mParm[i]);
   printf("\n");
 }
 
 return 0;
}

//********************   FIN DEFINITION PSFFFitStar   *********************



#include "starlist.cc" // since starlist is a template class 

template class StarList<PSFStar>;

template class StarList<PSFFitStar>;

int
Copy_BasetoPSFStarList(BaseStarList & stlse, PSFStarList & stlfit, 
		       PSF *psf)
{
  BaseStarCIterator it ;
  int n = 0 ;

  for(it = stlse.begin(); it != stlse.end(); it++)
    {

      BaseStar *star = *it;
      SEStar *starse = new SEStar(star->x, star->y, star->flux);
      PSFStar *starfit = new PSFStar(*starse, psf);
      stlfit.push_back(starfit);

      n++;
    }
  return(n);

}

int
Copy_PSFtoBaseStarList(BaseStarList & stlb, PSFStarList & stlfit)
{
  PSFStarCIterator it ;
  int n = 0 ;

  for(it = stlfit.begin(); it != stlfit.end(); it++)
    {

      PSFStar *star = *it;
      BaseStar *starb = new BaseStar(star->x, star->y, star->flux);
      stlb.push_back(starb);

      n++;
    }
  return(n);

}


int
Copy_PSFFittoSEStarList(PSFFitStarList & stlfit, SEStarList & stlse)
{
  PSFFitStarCIterator it ;
  int n = 0 ;

  for(it = stlfit.begin(); it != stlfit.end(); it++)
    {

      SEStar *star = (SEStar * ) *it;
      SEStar *starse = new SEStar( *star );
      stlse.push_back(starse);
      n++;
    }
  return(n);

}




// copie de SEStarlist ds PSFFitStarlist


int
Copy_SEtoPSFFit_StarList(SEStarList & stlse, PSFFitStarList & stlfit, 
			    PSF *psf)
{

  SEStarCIterator it ;
  int n = 0 ;

  for(it = stlse.begin(); it != stlse.end(); it++)
    {

      SEStar *starse = *it;
      PSFFitStar *starfit = new PSFFitStar(*starse, psf);
      stlfit.push_back(starfit);

      n++;
    }
  return(n);

}

#endif /*IS_IT_USEFUL */
