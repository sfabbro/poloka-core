// This may look like C code, but it is really -*- C++ -*-


/* DOC Definition of a PSF: it inherits from GeneralPSF2D, defined in
fct2dfit.h. GeneralPSF2D is a GeneralFunction daughter, at 2-Dimension`, and for wich the parameters Param are always classified in the order x, y, flux, other paramaters.

The different function that can be 
used as a PSF are defined in fct2dfit.h.

Definition of PSFStar: (SEStar daughter, + PSF) = a star
with a PSF 

Definition of PSFFitStar:
 (PSFStar Daughter, + Error Vector for a fit) = 
a star to be fitted

Associated StarList are also defined.

*/


#ifndef PSF_SEEN
#define PSF_SEEN


#include "image.h"
#include "sestar.h"

//#include "fct1dfit.h"
#include "fct2dfit.h"    
//#include "machine.h" /* for Int_4 */

// star qui herite de SEStar, donc a priori en coordonnees SE
// !!!! les routines de Add considere que la star est en XY
// les routines de fit retourne des coordonnees en  XY








//********************   DEFINITION PSF   *********************
/* The sole reasonable use of all this stuff is for simulation
of images, or to add fakes to real images. But anyway, do not
consider any serious simulation work without redesigning all that.
  In my view, one should really consider throwing all that stuff away:
It is far too complicated (a few thousand lines in total) to essentialy 
add Gaussians to images. I also suspect that all that code is not efficient.

P. Astier.
*/





class FuncPtrGFF : public GeneralPSF2D {
public:
  typedef double (*Func)(double dx, double dy, double const* par);
  FuncPtrGFF(Func f, int nPar);
  virtual        ~FuncPtrGFF();
  virtual double  Value(double const xp[], double const* parm);
  
protected:
  Func            mFunc;
  int             mNPar;
  double*		  mTmpParm;
private:
  FuncPtrGFF(const FuncPtrGFF &);
  FuncPtrGFF & operator=( const FuncPtrGFF  &);

};


//
// Une PSF est toujours associee a une GeneralFunction.
//

// la creation d'une PSF avec les types ci-dessous est 
// geree comme celle d'un singleton
// cela permet de l'utiliser ds les PSFStar, 
// ou la PSF n'est pas deletee par le destructeur


class PSF {
public:

  /*DOC {\bf Different kind of PSF are pre-defined}:\\
{\tt  enum PSFKind} {\\
   {\tt  psfCustom = 0}, (i.e. non pre-defined-type)\\ 
    {\tt   psfGaussienne = 1},\\
    {\tt   psfGaussInteg},
    {\tt   psfDelGaussienne},\\
     {\tt  psfDelGaussinteg},
    {\tt   psfDel1Gaussienne},\\
     {\tt  psfDel1GaussInteg},
     {\tt  psfMoffat},\\
     {\tt  psfMoffInteg},
     {\tt  psfDel2Gaussienne},\\
     {\tt  psfDel2GaussInteg}
  }; */


  enum PSFKind{
    psfCustom = 0, // PSF de type non determinee, perso par ex.
      psfGaussienne = 1,
      psfGaussInteg,
      psfDelGaussienne,
      psfDelGaussinteg,
      psfDel1Gaussienne,
      psfDel1GaussInteg,
      psfMoffat,
      psfMoffInteg,
      psfDel2Gaussienne,
      psfDel2GaussInteg
  };
  
  /*DOCF  allocation with pre-defined type.
   The PSF is then a singleton object (one allocation per pre-defined type one main program). */
  PSF(int psfKind, double ns=3); 


  /*DOCF  allocation with an already existing PSF.
No singleton management.*/
  PSF(GeneralPSF2D *psf2d, double ns=3);
  PSF(FuncPtrGFF::Func f, int nPar, double ns=3);
  virtual        ~PSF();
  
  /*DOCF returns PSF kind. */
  int           Kind() const {return mTypePSF;}

  /*DOC {\tt Value(double x, double y, double const* parV)}: 
returns PSF value. */
  virtual double  Value(double x, double y, double const* parV) const;
    
   /*DOCF returns PSF number of parameter. */
 int             NParm() const       {return mNParm;}
  /* returns the associated General Function */
 GeneralPSF2D* GetGFF() {return mGFF;}
 
friend class PSFStar;

protected:
  PSF(void);              // Pour relire une PSFStar

  int  mTypePSF;
  int  mNParm;          // Total, y compris x,y,fond_fit,flux.
  float    mNSig0;          // nb de sigmas pour negliger dans les images
                          // simulees
 
  GeneralPSF2D* mGFF;
  void  MakeStdGFF();

private:
  // pour gerer l'allocation d'une seule PSF par type pre-defini
  static GeneralPSF2D* _instance[10] ; // autant que de enum != 0 
 
  PSF(const PSF &);
  PSF& operator=(const PSF &); 
  
 };


//********************   FIN DEFINITION PSF   *********************


//********************   DEFINITION PSFStar   *********************






// Une etoile avec une PSF.
// Plusieurs etoiles peuvent partager la meme PSF

class PSFStar : public SEStar {

  /*DOC A SExtractor star with a PSF. Mainly for Monte-Carlo simulation. */

public:
  /*DOCF creation of a "blank" star, the PSF is given */
  PSFStar(PSF* psf=NULL);
  /*DOC {\tt   PSFStar(SEStar const & sestar, 
PSF* psf=NULL)}: copy of the SEStar, 
the given PSF is added */
  PSFStar(SEStar const & sestar, PSF* psf=NULL);
  PSFStar(double xx, double yy, double ampl, PSF* psf=NULL);
  PSFStar(PSFStar const & star);
  virtual        ~PSFStar();

  /* DOCF : background  */
  double&         Fond_fit()        {return fond_fit;}
  double const&   Fond_fit() const  {return fond_fit;}
 
  /* DOCF : PSF function value */
virtual double  operator() (double xx, double yy) const;
  /* DOCF : PSF vector parameters */
  double&         Parm(int i)        {return mParm[i];}
  double const&   Parm(int i) const  {return mParm[i];}
  /* DOCF : PSF Sigma X value */
  double&         SigX()       {return mParm[3];}
  double const&   SigX() const {return mParm[3];}
  /* DOCF : PSF Sigma Y value */
 double&         SigY()       {return mParm[4];}
  double const&   SigY() const {return mParm[4];}
  /* DOCF : PSF Rho value */
 double&         Rho()        {return mParm[5];}
  double const&   Rho()  const {return mParm[5];}
  double          Sig() const  {return sqrt(mParm[3]*mParm[3] + mParm[4]*mParm[4]);}
  /*DOC 
For the sake of convenience, some of the {\tt Parm(i)}
  have been renamed, as {\tt SigX, SigY, Rho, Fond_fit}. They are physically equivalent to {\tt Parm(3), (4), (5), (n-1)}.
Changing {\tt SigX} value is equivalent to change 
{\tt Parm(3)} value.
However, the star parameters as {\tt x, y, flux} 
{\em are not} are not re-named {\tt Parm(i)}.  
They are different variables, although they correspond 
to the same quantities as {\tt Parm(0), (1), (2)}. 
Thus, setting the star 
position x and y won't change the associated 
PSF parameters {\tt Parm(0)} and {\tt Parm(1)}.
This will be done only with calling the 
{\tt UpdParmVect()} method (see below). */


  /* DOCF PSF central (and not maximum) pixel value */
 virtual float     PixMax() const    {return (*this)(x,y);}

  /* DOCF add the star to an image (for simulation) */
  virtual void    ImgAdd(Image& img) const; 
  /* DOCF add the star to an image (for simulation) */
  virtual void    ImgAdd_2(Image& img) const; 
   /* DOCF substract the star from an image */
  virtual void    ImgSub(Image& img) const;
 
  /* DOCF recover the star's PSF */
  PSF*            GetPSF() const {return mPSF;}
   /* DOCF to set  star's PSF */
 void            SetPSF(PSF* psf);
  
   /* DOCF to update the parameters {\tt Parm} 
used by the PSF from
the parameters {\tt x, y, flux} of the star. */
  void            UpdParmVect()  const;
   /* DOCF update the parameters {\tt x, y, flux} of 
the star from
the PSF parameters {\tt Parm}. */
  void            UpdFromParm() const;

/*DOC read and write are similar to that of SEStar */
  virtual void    dump(ostream& s=cout) const;
  virtual void    dumpn(ostream& s=cout) const;

  virtual void    write(ostream& s=cout) const;
  virtual void    writen(ostream& s=cout)const ;

 
   virtual void    Read(istream& r, const char *Format, int &type_psf);
  // creation of the PSFStar with the given PSF
  static PSFStar* read(istream& r, const char *Format, PSF *mypsf);
  // creation of the PSFStar with the type given in the list
  static PSFStar* read(istream& r, const char *Format);

  //!
  //  const char * WriteHeader_(ostream & pr = cout, const char*i = NULL) const;
  std::string WriteHeader_(ostream & pr = cout, const char*i = NULL) const; // NRL 04/2004

  
protected:
  PSF*            mPSF;
  double*         mParm;
  double fond_fit ;


private:

  PSFStar& operator=(const PSFStar &); 

};

Image& operator += (Image& img, PSFStar const& psf);
Image& operator -= (Image& img, PSFStar const& psf);



// PSFStarList definition(s) : 
#include "starlist.h"
typedef StarList<PSFStar> PSFStarList;
//typedef list<PSFStar*>::const_iterator PSFStarCIterator;
//typedef list<PSFStar*>::iterator PSFStarIterator;
typedef PSFStarList::const_iterator PSFStarCIterator;
typedef PSFStarList::iterator PSFStarIterator;
typedef CountedRef<PSFStar> PSFStarRef;

BaseStarList* PSF2Base(PSFStarList * This);
const BaseStarList* PSF2Base(const PSFStarList * This); 



//********************   FINDEFINITION PSFStar   *********************


//********************   DEFINITION PSFFitStar   *********************



#ifdef IS_IT_USEFUL

// Une etoile pour faire des fits
class PSFFitStar : public PSFStar {
/* DOC A SEStar + PSF + method for fitting = a star for seeing computation ! */

public:
/* DOC The creation can be made giving a PSF, a
PSFStar, or a SEStar and a PSF.  */
  PSFFitStar(PSF* psf=NULL);
  PSFFitStar(PSFFitStar const & star);
  PSFFitStar(PSFStar const & star);
  PSFFitStar(SEStar const & star, PSF *psf);

  virtual        ~PSFFitStar();
  
/* DOCF: error vector associated to the {\tt Parm} parameters
vector. */
  double&         Err(int i)        {return mErr[i];}
  double const&   Err(int i) const  {return mErr[i];}
 /* DOCF: error on flux  */
 double&         ErrFlux()         {return mErr[0];}
  double const&   ErrFlux()  const  {return mErr[0];}
 /* DOCF: error on x  */
  double&         ErrX()            {return mErr[1];}
  double const&   ErrX()  const     {return mErr[1];}
  /* DOCF: error on y */
 double&         ErrY()            {return mErr[2];}
  double const&   ErrY()  const     {return mErr[2];}
  /* DOCF: error on Fond_fit  */
 double&         ErrFond()         {return mErr[mPSF->NParm()-1];}
  double const&   ErrFond()  const  {return mErr[mPSF->NParm()-1];}
  /* DOCF: error on SigX  */
 double&         ErrSigX()         {return mErr[3];}
  double const&   ErrSigX()  const  {return mErr[3];}
  /* DOCF: error on SigY  */
 double&         ErrSigY()         {return mErr[4];}
  double const&   ErrSigY()  const  {return mErr[4];}
 /* DOCF: error on Rho  */
  double&         ErrRho()          {return mErr[5];}
  double const&   ErrRho()   const  {return mErr[5];}
 /* DOCF: $\chi^2$ of fit  */
  double&         Xi2()             {return mXi2;}
  double const&   Xi2()      const  {return mXi2;}
  
/* DOC All these parameters are re-naming of  {\tt Err(i)}. */

/* DOC read and write are identical to that of SEStar */
  virtual void    dump(ostream& s=cout) const ;
  virtual void    dumpn(ostream& s=cout) const;

  virtual void    write(ostream& s=cout) const;
  virtual void    writen(ostream& s=cout) const;


  virtual void   Read(istream& r, const char *Format, int & type_psf);
  static PSFFitStar* read(istream& r, const char *Format, PSF *mypsf); 
  static PSFFitStar* read(istream& r, const char *Format);

  const char * WriteHeader_(ostream & pr = cout, const char *i= NULL) const;
 

 /* DOC //{\bf fitting of the PSF on an image}: //
  virtual int     Fit(Image const& img, 
Image const& noiseimg, int lp=0);
 */
  virtual int     Fit(Image const& img, Image const& noiseimg, int lp=0);
  /* DOCF initialisation of the 
fit with SExtractor computed x, y, flux, fond, Mxx, Myy, 
Mxy of the star.  */
 virtual int     IniFitSE(bool fond_force=false, double fond=0., int lp=0);
/* DOC // initialisation of the fit with x, y and a rough computation
of the flux and background on the image. */
  virtual int     IniFit(Image const& img, Image const& noiseimg, int lp=0);
protected:
  double          mXi2;
  double*         mErr;

private:

  PSFFitStar& operator=(const PSFFitStar &); 

};






//********************   FIN DEFINITION PSFFFitStar   *********************


/* what concerns the PSFFitStarList's : */





typedef StarList<PSFFitStar> PSFFitStarList;
typedef list<PSFFitStar*>::const_iterator PSFFitStarCIterator;
typedef list<PSFFitStar*>::iterator PSFFitStarIterator;



/*DOCF: to copy a SEStarList in a PSFFitStarList*/
int
Copy_SEtoPSFFit_StarList(SEStarList & stlse, PSFFitStarList & stlfit, PSF *psf);

int
Copy_BasetoPSFStarList(BaseStarList & stlse, PSFStarList & stlfit, 
		       PSF *psf);
int
Copy_PSFFittoSEStarList(PSFFitStarList & stlfit, SEStarList & stlse);

int
Copy_PSFtoBaseStarList(BaseStarList & stlb, PSFStarList & stlfit);

#endif /* IS_IT_USEFUL */


#endif
