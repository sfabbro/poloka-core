#include <iostream>

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif

//#define DEBUG
//#define TEST

#include <iostream>
#include "basestar.h"
#include "sestar.h"
#include "image.h"
#include "fitsimage.h"
#include "gtransfo.h"




// =================================================================
// =================================================================
//
// definition de la classe AperturePhotomSEStar
//
// =================================================================
// =================================================================

class AperturePhotomSEStar : public SEStar
{
 public:

#define nb_rayons 20

  float rayon[nb_rayons]; // serie de nb_rayons rayons d'ouverture differents
  float aperture_flux[nb_rayons];  // les flux dans ces rayons

  double ra;
  double dec; 

  float seeing;
  float exptime;
  float fluxmax;
  float fluxmax2;

  float delta_x; // decalage en x (pixels) du barycentre de l'objet par rapport a la valeur sextractor
  float delta_y; // decalage en y (pixels) du barycentre de l'objet par rapport a la valeur sextractor
  
 public:
  AperturePhotomSEStar(){return;};
  AperturePhotomSEStar(const SEStar &sestar) : SEStar(sestar)
    {delta_x=0.; delta_y=0.; };   //! constructor
  ~AperturePhotomSEStar(){};
  
  // initialisation des rayons d'ouverture, calcules en fonction du seeing
  void SetSeeing(float bob)
    {
      seeing = bob;
      for (int i=0; i<nb_rayons; ++i) 
	{
	  rayon[i] = (float) (i+1) *seeing; 
	}
    }
  bool FillFlux(FitsImage& fitsimage, FitsImage& fitsweight, int, int);
  bool FillFlux_noweight(FitsImage& fitsimage, double, int s=0, int i=0);
  bool ComputeBarycenter(FitsImage& fitsimage, FitsImage& fitsweight);
  bool ComputeBarycenter_noweight(FitsImage& fitsimage, int, int);

  float compute_platitude();
  void write(ostream& s = cout, const bool write_header=false);
  void write_short(ostream& s = cout, const bool write_header=false);
  void write_profile(ostream& s = cout, const bool write_header=false);
};


typedef StarList<AperturePhotomSEStar> AperturePhotomSEStarList;
typedef AperturePhotomSEStarList::const_iterator AperturePhotomSEStarCIterator;
typedef AperturePhotomSEStarList::iterator AperturePhotomSEStarIterator;
typedef CountedRef<AperturePhotomSEStar> AperturePhotomSEStarRef;

// =================================================================
// =================================================================
//
// definition de la classe AperturePhotomBaseStar
//
// =================================================================
// =================================================================

class AperturePhotomBaseStar : public BaseStar
{
 public:

#define nb_rayons 20

  float rayon[nb_rayons]; // serie de nb_rayons rayons d'ouverture differents
  float aperture_flux[nb_rayons];  // les flux dans ces rayons

  double ra;
  double dec; 

  float seeing;
  float exptime;
  float fluxmax;
  float fluxmax2;

  float delta_x; // decalage en x (pixels) du barycentre de l'objet par rapport a la valeur sextractor
  float delta_y; // decalage en y (pixels) du barycentre de l'objet par rapport a la valeur sextractor
  
 public:
  AperturePhotomBaseStar(){return;};
  AperturePhotomBaseStar(const BaseStar &sestar) : BaseStar(sestar)
    {delta_x=0.; delta_y=0.; };   //! constructor
  ~AperturePhotomBaseStar(){};
  
  // initialisation des rayons d'ouverture, calcules en fonction du seeing
  void SetSeeing(float bob)
    {
      seeing = bob;
      for (int i=0; i<nb_rayons; ++i) 
	{
	  rayon[i] = (float) (i+1) *seeing; 
	}
    }

#define default_idebut 12
  bool FillFlux(FitsImage& fitsimage, FitsImage& fitsweight, int, int);
  bool FillFlux_noweight(FitsImage& fitsimage, double, int i=0, int s=0);
  bool ComputeBarycenter(FitsImage& fitsimage, FitsImage& fitsweight);
  bool ComputeBarycenter_noweight(FitsImage& fitsimage);
  float compute_platitude(double &mean, double &var, int i_debut = default_idebut);
  float compute_platitude(int i_debut = default_idebut);
  float compute_indep_platitude(double &mean, double &var, int i_debut = default_idebut);
  float compute_indep_platitude(int i_debut = default_idebut);
  void write(ostream& s = cout, const bool write_header=false);
  void write_short(ostream& s = cout, const bool write_header=false);
  void write_profile(ostream& s = cout, const bool write_header=false);
};


typedef StarList<AperturePhotomBaseStar> AperturePhotomBaseStarList;
typedef AperturePhotomBaseStarList::const_iterator AperturePhotomBaseStarCIterator;
typedef AperturePhotomBaseStarList::iterator AperturePhotomBaseStarIterator;
typedef CountedRef<AperturePhotomBaseStar> AperturePhotomBaseStarRef;


// =================================================================
// =================================================================
//
// definition de la classe MultiMeasuredStar
//
// =================================================================
// =================================================================

class MultiMeasuredStar : public BaseStar
{
 public:

  struct mesure
    {
      float measured_aperture_flux[nb_rayons];
      float fluxmax;
      float fluxmax2;
      float measured_errors;
      float measurement_radius[nb_rayons];
      bool valid;
      bool is_final;
      float final_flux;
      float seeing;
  };
  mesure mesures[100];
  int nb_measurements;
  int nb_valid;

 public:
  MultiMeasuredStar() {return;};
  MultiMeasuredStar(const BaseStar bstar) : BaseStar(bstar)
    {nb_measurements=0; };   //! constructor
  ~MultiMeasuredStar(){};
  
  float compute_mean(int measurement);
  float compute_mean2(int measurement);
  void finalize(int measurement);
  float compute_platitude(int measurement);
  void AddInfo(float rayon_mesure[], float aperture_flux[], float eflux, bool val, float fluxmax, float fluxmax2);
  void AddInfo(float rayon_mesure[], bool must_be_false);
  void AddInfo(bool must_be_false);
  void AddInfo(AperturePhotomSEStar *aps);
  void write_valid(ostream& s, int star_id, bool write_header);
  void write_all(ostream& s, int star_id, bool write_header);
  bool match_aps(const SEStarList& laliste,  Gtransfo * pix2radec, AperturePhotomSEStar* aps);

};

typedef StarList<MultiMeasuredStar> MultiMeasuredStarList;
typedef MultiMeasuredStarList::const_iterator MultiMeasuredStarCIterator;
typedef MultiMeasuredStarList::iterator MultiMeasuredStarIterator;
typedef CountedRef<MultiMeasuredStar> MultiMeasuredStarRef;
