// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMSNSTAR__H
#define SIMSNSTAR__H


#include <string>
#include "imagepsf.h"
#define NO_PERS
#ifndef NO_PERS
#include "persistence.h"
#endif
#include "basestar.h"
#include "sestar.h"

class Image;
class Gtransfo;

/*! \file */

/*! 

1) A model star to be added for simulating a supernova (ModelStar, inheriting from BaseStar).

The principle is to cut a stamp  around the 
model star, shift it to the the supernova position, and add 
it to the   image. However, as we are manipulating pixels, 
the shift in x and y has to be an INTEGER. 

The main problem is that a same supernova is going to be 
added on different images, e.g. all images of the same field
taken in a night (often refered to as "new" images).
If we wants that on EVERY new image, the position of the fake sn 
= the position of the model + an INTEGER pixel shift, 
then it means we are assuming that the geometrical transformation 
between 2 of these images is a translation. 
If the integer shift was computed
in the reference image frame of coordinate, then this assumption must holds
for the new images AND the ref image. That's why we provide the possibility to
compute this shift in the coordinates system of 1 of the new 
 images instead of the
ref system.
Of course, a minor rotation can be tolerated.



The model star x,y,flux of the model star refer
to a reference system, geometric and photometric.
The position of the fake supernova in the ref system is thus:
xsn = x + xShift, ysn = y + yShift   and its flux is :
flux * photFactor.

When adding the fake on a another image, one needs
the geometrical transformation tf(ref->image).
The model star is then in (X,Y) = tf(x,y),
the ture supernova position is (Xsn,Ysn)=  tf(xsn,ysn).
If tf is close enough to a translation, the supernova position 
on the image (Xsn,Ysn) is well approximated by:
X +  xShift, Y + yShift.
The pixels of the model star
around the position (X,Y) are cut, shift of an integer number of pixels
along x and y axis,  multiplied by photFactor an d paste.

In a nutshell,  to be added on an image, a ModelStar 
only needs a transfo between the 
coordinates system of its position and the coordinates system of the image.
This is done in AddToImage.


2) A Simulated supernova (SimSNStar, inheriting from BaseStar).
This regroups all the generated parameters: position, flux, 
magnitude, host caracteristics.

How to randomly generate a SimSNStar is put in "simulation.h".

3) A Simulated supernova + a Model Star (SimSNWModelStar, inheriting from SimSNStar)
Indeed, if one chooses to add
a supernova using a ModelStar, then the position of the 
supernova must be slightly
modified so that the shift between the SN position and the
Model position be an integer number of pixels.
Methods to marry a Simulated supernova to a Model Star are 
provided, along with the method to add it on an image.

How to randomly generate a SimSNWModelStar 
and marry it to a ModelStar is described in "simulation.h".


4)Simulated Supernova W Model Star List
A method to add it on an image, using of course the Model method.
It also writes a debug file, so that one can check 
the "integer shift" assumption : e.g., it gives 
the fractional part of the real shift between the model 
position and and the expected SN position according to the transfo.


*/

class ModelStar : public BaseStar
{
  private:
  int xShift, yShift; 
  double stampsize ;
  double photFactor; 



  // ******* for debug only, not written in output  
  // Flux of model star on current image
  double flux_i; // flux_i * coeff phot (ref->image_i) doit etre = flux_ref
  double eflux_i;
  double fwhm_i; 
  // pour check car tfref->image_i(xsn_ref,ysn_ref) 
  // = (xfake_i,yfake_i) doit etre = tf(xmod_ref,ymod_ref) + xShift .
  double xfake_i;
  double yfake_i;

  
 public:



  double PhotFactor() const {return photFactor;}
  double XShift() const {return xShift ;}
  double YShift() const {return yShift ;}

  // for debug only, not written in output
  double Flux_i() const {return flux_i ;}
  double& Flux_i() {return flux_i ;}
  double EFlux_i() const {return eflux_i ;}
  double& EFlux_i() {return eflux_i ;}
  double Fwhm_i() const {return fwhm_i ;}
  double& Fwhm_i() {return fwhm_i ;}
  double XFake_i() const {return xfake_i ;}
  double& XFake_i() {return xfake_i ;}
  double YFake_i() const {return yfake_i ;}
  double& YFake_i() {return yfake_i ;}
  
  ModelStar();
  ModelStar(const BaseStar &AStar, const double size, const int Dx = 0, 
	    const int Dy = 0, const double PhotFactor=0);

   //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;
  //! for dump
  virtual void    dump(ostream& s = cout) const ;

  virtual void writen(ostream & pr = cout) const;
  std::string WriteHeader_(ostream &pr = cout, const char* i = NULL) const;
  //  virtual void    read_it(istream& r, const char *Format); 
  //  static ModelStar* read(istream& r, const char* Format);

  //! Addition to an image. Needs a transfo as explained above. Transfo = tf(ref->image). Update the saturation map according to new pixel values.
  void AddToImage(const Image &image, Image & dest, const Gtransfo *Transfo, 
		  Image * psat=NULL, double satlevel=-1 ,double gain=1) const ;

};


class SimSNStar : public BaseStar 
{


  private:
  /*! concerning supernova: x, y, flux are in ref system */
  double mag_sn ;
  /*! concerning host galaxy */
  double mag_gal;
  double x_gal, y_gal ;
  double a_gal ;
  double  fluxmax_gal; 
 
 public:

  double Mag_SN() const {return mag_sn;}
  double& Mag_SN() {return mag_sn;}
  double Mag_Gal() const {return mag_gal;}
  double& Mag_Gal() {return mag_gal;}
  double X_Gal() const {return x_gal;}
  double& X_Gal() {return x_gal;}
  double& Y_Gal() {return y_gal;}
  double Y_Gal() const {return y_gal;}
  double A_Gal() const {return a_gal;}
  double& A_Gal() {return a_gal;}
  double Fluxmax_Gal() const {return fluxmax_gal;}
  double& Fluxmax_Gal() {return fluxmax_gal;}
  
  SimSNStar();
  SimSNStar(double xx, double yy, double ff);
  virtual void writen(ostream & pr = cout) const;
  std::string WriteHeader_(ostream &pr = cout, const char* i = NULL) const;
  virtual void    dumpn(ostream& s = cout) const;
  virtual void    dump(ostream& s = cout) const ;
  //  virtual void    read_it(istream& r, const char *Format); 
  //  static SimSNStar* read(istream& r, const char* Format);

  void NewFlux(double NewZeroPoint){flux = pow(10,(NewZeroPoint - mag_sn)*0.4);}

};










/* a Simulated Supernova, and a model star to add it on an image. */

class SimSNWModelStar : public SimSNStar 
{  

 public:

  ModelStar model_on_ref ;

  SimSNWModelStar() : SimSNStar(){};
  SimSNWModelStar(double xx, double yy, double ff) : SimSNStar(xx,yy,ff){};
  SimSNWModelStar(const SimSNStar &AStar) : SimSNStar(AStar){};
  virtual void writen(ostream & pr = cout) const;
  virtual void    dumpn(ostream& s = cout) const;
  virtual void    dump(ostream& s = cout) const ;
  std::string WriteHeader_(ostream &pr = cout, const char* i = NULL) const;
  //  virtual void    read_it(istream& r, const char *Format); 
  //  static SimSNWModelStar* read(istream& r, const char* Format);

  //! As indicated in Name : select the closest star in SEStarList
  // of beautiful stars, and compute the integer shift, the photometric ratio
  // !! Modifies the SimSNStar coordinates so that shift be integer.
  // the shift is computed not in the ref frame (coordinate of the 
  //  SimSNStar and of the StarList) but in the coordinates obtained
  // with the provided transfo.
  void MariageToAModelStar(SEStarList const & BellesEtoiles, 
			   const Gtransfo *Transfo,
			   const Gtransfo *TransfoInv);
  // same but transfos are set to identity
  void MariageToAModelStar(SEStarList const & BellesEtoile);

  //! Add fake SN to an image, with the model method. Has to be married before.
  //! Transfo is from the ref frame to the image.
  void AddWModelToImage(Image &image, const Gtransfo *Transfo,   
			Image * psat = NULL, 
			double satlevel=-1 ) const {
    model_on_ref.AddToImage(image,image,Transfo, psat, satlevel);}

  void AddWModelToImage(const Image &image, Image & dest, 
			const Gtransfo *Transfo,   
			Image * psat = NULL, 
			double satlevel=-1 ) const {
    model_on_ref.AddToImage(image,dest,Transfo, psat, satlevel);}

};






// should be static, but is rightnow used elsewhere for debug.
int integer_delta(double xsn, double xmodel);

/*** definitions  StarList **********/

#include "starlist.h"


class SimSNStarList : public StarList<SimSNStar> 
{
  public :
  void NewFlux(double NewZeroPoint);
  void AddWGaussianToImage(double sigmax, double sigmay, double rho,
			   Image & dest, 
			   const Gtransfo *Transfo,   
			   Image * psat = NULL, double satlevel=-1) const ;
  void AddWGaussianToImage(double sigmax, double sigmay, double rho,
			   Image & dest, 
			   Image * psat = NULL, double satlevel=-1) const ;

  
};


typedef SimSNStarList::const_iterator SimSNStarCIterator;
typedef SimSNStarList::iterator SimSNStarIterator;
typedef CountedRef<SimSNStar> SimSNStarRef;

#ifndef SWIG
//! type casting
BaseStarList* SimSN2Base(SimSNStarList * This);
const BaseStarList* SimSN2Base(const SimSNStarList * This);
#endif




class SimSNWModelStarList : public StarList<SimSNWModelStar> 
{
  public :
  //! Adds a fake supernova list to an image
  //! Transfo is from the frame where the model was chosen (i.e. ref) to the image
  //! saturation level is that of the image
  //! writes a debug file, to check if the model method is sound :
  //! gives the fractional part of the position between  the expected 
  //! sn position on image, computed with the geom. tf., 
  //! and the place where it is stuck. see code.
  void AddWModelToImage(const Image &image, Image & dest, 
			const Gtransfo *Transfo,   
			Image * psat = NULL, double satlevel=-1,
			bool print_debug=true) const ;
  void AddWModelToImage(Image &image, 
			const Gtransfo *Transfo,   
			Image * psat = NULL, double satlevel=-1,
			bool print_debug=true) const {
    AddWModelToImage(image,image,Transfo,psat,satlevel,print_debug);}

  void NewFlux(double NewZeroPoint);

};



typedef SimSNWModelStarList::const_iterator SimSNWModelStarCIterator;
typedef SimSNWModelStarList::iterator SimSNWModelStarIterator;
typedef CountedRef<SimSNWModelStar> SimSNWModelStarRef;




#ifndef SWIG
BaseStarList* SimSNWModel2Base(SimSNWModelStarList * This);
const BaseStarList* SimSNWModel2Base(const SimSNWModelStarList * This);
#endif
 


// DaoPhot utilitaries
// The place of this is to be discussed with Seb


#include "image.h"
#include "imagepsf.h"
class DaoPsf;


void AddWDaoPsfToImage(DaoPsf const & daopsf, double xc, double yc, 
		       double flux, Image & image, Image * psat, double saturation);

void AddListWDaoPsfToImage(DaoPsf const & daopsf, BaseStarList *List, 
		       Image & img, 
		       Image * psat=NULL, double saturation=-1);
		       
void AddWPsfToImage(ImagePSF &psf ,double xc, double yc, 
		       double flux,Image & image,
		       Image * psat,
		       double saturation);
void AddListWPsfToImage(ImagePSF &psf, BaseStarList *List,
		       Image & img,  Image * psat,
		       double saturation);		    
#endif /* SIMSNSTAR__H */
