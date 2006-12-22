// This may look like C code, but it is really -*- C++ -*-

#ifndef SIMULATION__H
#define SIMULATION__H

#include "persistence.h"

#include "datacards.h"
#include "reducedimage.h"
#include "simsnstar.h"
#include "gtransfo.h"

class ImageGtransfo ;

typedef enum MethodeSim{Random =0, Damier, InHost, AdaptedToHost};


/*! \page effic Preparation of the fake supernova starlist.

1) DatSim and the method for generating a fake supernova SimSNStar,
that is : the generation of the supernova caracteristics (position,flux),
its eventual association to an host galaxy. 
No connexion with a model star is done at this point.

2)ForSim regroups what is necessary to generate the fake supernova SimSNStar:
It thus contains a DatSim, and a SEStarList of host galaxy + zero_point.
The creator reads the DatSim in the TOADSCARDS sub.datacard.
There is also a very simple creator provided, if one wants to do things
personnaly.

3)ForSimWModel inheriting from ForSim, contains what
 necessary to generate supernovae and marry them to a ModelStar.

It thus contains a SEStarList of selected stars to become model stars
(zero_point is that of ForSim), and a geometrical transfo.

Note that rightnow, both SEStarList for stars and galaxies 
are supposed to be in the same geom. frame, 
and the same photometric system : only one zero_point is given for both.

The geometrical transfo. is used when one wants to compute
the integer shift in position between the fake supernova and the 
model star : it can be computed either in the geom. frame of the given 
SEStarList (corresponding usually to the  reference frame ), either 
in the geom. frame of an image of the stack of new.
The reason for this was explained in the fake star section :
 the integer shift method relies
on the assumption that the geom. transfo between the images 
are close to a translation. This may holds better for the 
    images of the new stack (where the fake is added )
than for the new images + the ref image.

The method MakeListSn generate a list of fake supernovae associated 
to a model star.

4) routine(s), that could be put elsewhere for the selection
  of ModelStar.
   Rightnow there is one using a ReducedImageList i.e.
    mainly a list of catalogs. Note that it uses a few static routine
 for the selection of ModelStars, that, if proven useful elswhere in the code 
in the future, could be "un-static-fied".


*/

//! Simulation Data Cards. It regroups information
//! on the supernova itself (total number, position, magnitude)
//! and on a eventual host (magnitude cut to select host galaxies)
//! magnitude information depends are provided in filters g,r,i.



struct DatSim {
    bool Islcsim;
    MethodeSim Methode ;
    int numberOfFakes;
    int delta_x ;
    int delta_y ;
    int xsize_damier ;
    int ysize_damier ;
    //! sn mag generated in [minMag, maxMag]
  double minMag, maxMag;
    //! select galaxies above maxGalMag (ie mag < maxGalMag)
    double maxGalMag;
    //! objects with mag > LimGalMag are all galaxies
    double LimGalMag ;
    // in case of adapted generation: sn mag will be in [host mag +  deltaMagmin, host mag + deltaMagmax] AND in [minMag, maxMag]
    double deltaMagmin, deltaMagmax;
  bool hasZeroPoint;
  double zeroPoint;
  
  DatSim() { Islcsim= false ;numberOfFakes = 500; xsize_damier =0 ;ysize_damier=0; delta_x = 100 ; delta_y = 100 ; minMag = 22; maxMag = 26; maxGalMag=26 ;  LimGalMag = 24 ; deltaMagmin=0. ; deltaMagmax=1. ;}
    //! read parameters for simulation in datacards. magnitude are set according filter (default = i)
  DatSim(const string &FileName, string const & filter = "");
  void Read(const string &FileName, string const & filter = "");
  void Print();

//! generate a random flux according to datacard specification.
  SimSNStar *RandomFluxSN(double xsnref, double ysnref, double zerop_ref );
//! generate a random position according to datacard specification.
  SimSNStar *RandomSN(int Nx, int Ny,  double zerop_ref );
//! generate a random position and flux in a galaxy according to datacard specification.
  SimSNStar *RandomSNInGalaxy(SEStar & gal_onref ,double zerop_ref );

  bool HasZeroPoint() const { return hasZeroPoint;}
  double ZeroPoint() const { return zeroPoint;}

    private :
  void LitDataCards(DataCards &, string const & filter = "");

};


//! used to select model stars. 
//! stars are first chosen in stlref
//! and then are asked to be present and ok in all the
//! the catalogs of the ReducedImageList.
//! the reducedimages are supposed to be aligned on the stlref SEStarList.
//! Note that this doesn't work if the images don't overlap well : 
//! the intersection of the catalogs can then be empty !

// a changer de place (dans fakestar ? dans sestar ?)
// des que je serai fixee DH

void 
PreSelectModelStars(SEStarList const &starList, 
		    SEStarList &BellesEtoiles, const Image &model);



bool
SelectModelStars(ReducedImageList  & imglist, SEStarList const & stlref, 
		 SEStarList &BellesEtoiles, const Image &model);

#ifdef STORAGE
void 
SelectionOfModelStars(SEStarList const &starList,
		      SEStarList &BellesEtoiles, Image * model);
#endif


//! ForSim regroups all the necessary elements to perform the generation
class ForSim {

public :   

  //For the SN carateristics generation :
  // *************************************
  //! read from defaults or copied from a given DatSim 
  DatSim datsim ;
 //!  host gal in coor of ref.      
  SEStarList BellesGal;
  //! zero point for starlist in which model stars an host galaxies are selected. supposed to be the same, but this can be changed easily.
  double zero_point;
  //! size of the frame where SN will be simulated (=ref frame)
  int XrefCCD, YrefCCD ;
  //! facultatif
  //! size of the image frame where the SN are added (for grid simulation).
  int XCCD, YCCD ;
  CountedRef<Gtransfo> RaDecToPix;

  
 
 
  //! creator reads DatSim from default datacard, select model star
  //! in the list  ListForModelStars, finds the related zero-point etc. in 
  //! the RefImage header, select host gal in the ref image catalogue 
  //! for debug, a fits image with the used host gals is saved.
  ForSim(ReducedImage const & RefImage);
    
  //! this creator does nothing. ForSim can be filled afterwards as 
  //! done in above  creator with Fill method.
  ForSim(){zero_point=0.;XrefCCD=0; YrefCCD=0; XCCD=0; YCCD=0;}
  //! for debug, a fits image with the used host gals is saved.
  void Fill(ReducedImage const & RefImage);  
//! Does it.
  void MakeListSN(SimSNStarList & SNList);
  SimSNStarList *MakeListSN(){SimSNStarList *list = new SimSNStarList; MakeListSN(*list);return(list);}

private :
  //! as indicated in Name. 
  void Construct_SNList_Random(SimSNStarList & SNList) ;
  void Construct_SNList_Damier(SimSNStarList & SNList);  
  void Construct_SNList_WHost(SimSNStarList & SNList);

  void ConstructSNwithList(SimSNStarList & SNList);

};
    

//! ForSimWModel regroups all the necessary elements to perform the generation
// and the association to a model star.
class ForSimWModel : public ForSim {

public : 

  //For the mariage to a model star:
  //*********************************
  //! model stars in coor of ref.
  SEStarList BellesEtoiles; 
  //! geom transfo from ref to an image of the newstack, so that the integer shift is computed in the image coordinate system.
  GtransfoRef Transfo ; // transfo de ref vers image
  GtransfoRef TransfoInv ;

 
 
  //! creator reads DatSim from default datacard, select model star
  //! in the list  ListForModelStars, finds the related zero-point etc. in 
  //! the RefImage header, select host gal in the ref image catalogue
  //! for debug, a fits image with the used model stars is saved. 
  ForSimWModel(ReducedImage const & RefImage, 
	       SEStarList const & ListForModelStars, 
	       const ImageGtransfo* tf=NULL);

  ForSimWModel(ReducedImage const & RefImage, 
	       const SEStarList &ListForModelStars,
	       const SEStarList &ListForGalaxies,
	       const ImageGtransfo* tf=NULL);

    
  //! this creator does nothing. ForSim can be filled afterwards as 
  //! done in above  creator with Fill method.
  ForSimWModel(){}
  
  //! 
  ForSimWModel(ForSim const & aforsim):ForSim(aforsim){}
  //! 
  void Fill(ReducedImage const & RefImage, SEStarList const & ListForModelStars, 
	 const ImageGtransfo* tf=NULL);
    
  //If you need another creator, feel free to add yours.

  //! Does it.
  void MakeListSNWModel(SimSNWModelStarList & SNList);
  SimSNWModelStarList *MakeListSNWModel(){SimSNWModelStarList *list = new SimSNWModelStarList; MakeListSNWModel(*list);return(list);}
};
    



#endif
