#include "sestar.h"
#include "image.h"
#include "reducedimage.h"
#include "toadscards.h"
#include "senearstar.h"
#include "dodetection.h"
#include "fitsimage.h"
#include "datacards.h"

// All simulation stuff
// starFinder which uses an histo shape (log (flux/fluxmax)) vs fwhm to find stars
//

#include <stdlib.h> // for random
#ifdef __sun
extern "C" {
extern long random(void);
}
#endif

static bool IncreasingCStar(const SEStar *S1, const SEStar *S2)
{
  return (S1->Cstar() < S2->Cstar());
}

#ifdef STORAGE
static bool DecreasingCStar(const SEStar *S1, const SEStar *S2)
{
return (S1->Cstar() > S2->Cstar());
}
#endif


struct DatSim {
  int numberOfFakes;
  double minMag, maxMag;
  
  DatSim() { numberOfFakes = 100; minMag = 22; maxMag = 26;}
  void LitDataCards(DataCards &);
  DatSim(const string &FileName);
  void Print();
};



DatSim::DatSim(const string &FileName)
{
  if (FileExists(FileName))
    {
      DataCards cards(FileName);
      LitDataCards(cards);
    }
  else
    {
      cerr << " DatSim::DatSim Cannot open  FileName " << endl;
    }
}

// Reads the parameters for the simulation
void DatSim::LitDataCards(DataCards &Data)
{
  numberOfFakes = Data.IParam("NUMBER_OF_FAKES");
  minMag = Data.DParam("MIN_FAKE_MAG");
  maxMag = Data.DParam("MAX_FAKE_MAG");
  if (minMag > maxMag) swap (minMag,maxMag);
}

void DatSim::Print()
{
  cout << " FAKES : number  of fakes : " <<  numberOfFakes 
      << " mag range : ["<< minMag << ',' << maxMag << ']' << endl;
}





static void KillBadObjects(SEStarList &List)
{
  cout << "Nombre d'objects avant nettoyage " << List.size()<<endl;
  // On degage les objets a probleme
  for (SEStarIterator it= List.begin(); it != List.end();) 
    {
      SEStar *star = *it;
      if (star->Flag() != 0) 
	it = List.erase(it);
      else ++it;
    }
  cout << "Nombre d'objects apres nettoyage " << List.size()<<endl;
}

// selectionne datSim.numberOfFakes Belles Etoiles
// dans starList
static void 
Construct_Star(SEStarList &starList,SEStarList &BellesEtoiles, DatSim & datSim)
{


  // Vire les objets a problemes
  starList.sort(DecreasingFlux); 
  KillBadObjects(starList);


  // On cherche les belles etoiles dans l'image
  SEStarList etoiles;
  starList.sort(DecreasingFlux);
  StarFinder(starList,etoiles);

  // On vire les etoiles en dessous d'un certain signal sur bruit : 100
  cout << "Taille avant coupure en stnoise " << etoiles.size() << endl;
  for (SEStarIterator it= etoiles.begin(); it!= etoiles.end() ;)
    {
      // + coupure en signal sur bruit
      SEStar *mystar = *it;
      double STNoise = mystar->flux / mystar->EFlux();
      double x = mystar->x; 
      double y = mystar->y;

      // on verifie que l'objet n'est pas sur un bord, corrige un probleme
      // de la campagne CFHT01A
      if (STNoise < 100 || x<1.e-3 || y<1.e-3) it = etoiles.erase(it);
      else ++it;
    }
  cout << "Taille apres coupure en stnoise " << etoiles.size() << endl;
  cout << "Nombre de fakes " << datSim.numberOfFakes << endl;


  // On remplit une liste du nombre d'etoiles desire
  int starcont=0;
  do 
    {
      for (SEStarIterator it= etoiles.begin(); it!= etoiles.end() ; ++it)
	{
	  if (starcont == datSim.numberOfFakes) break;
	  SEStar *star = *it;
	  BellesEtoiles.push_back(new SEStar(*star));
	  ++starcont;
	  ++it; 
	  if (it == etoiles.end()) it = etoiles.begin();
	}
    } while (starcont != datSim.numberOfFakes);

  
}


static void 
Construct_Gal(SEStarList &galaxies, SEStarList & GalList,
	      double & expoNew, double & zerop, DatSim & datSim)
{
  cout << "Processing Fakes associated with galaxies " << endl;
  
  // Construction de la liste de galaxie a partir du stack des new
  
  GalList.sort(DecreasingFlux); 
  cout << " nombre d objets lus " << GalList.size() << endl;
  
  KillBadObjects(GalList);

  
  // keep only the (galaxies.size() - datSim.numberOfFakes most galaxy like.
  GalList.sort(IncreasingCStar);
  GalList.CutTail(datSim.numberOfFakes);
  
  cout << " nombre de galaxies " << GalList.size()  << endl;

  int galcont = 0;  
  
  for (SEStarIterator it= GalList.begin(); it != GalList.end();++it ) 
    { 
      SEStar *star = *it;
      
      // Caracteristiques des objets
      double star_flux=star->flux;
      double mag = zerop - 2.5 * log10(star_flux / expoNew);
      
      // construction de la liste des galaxies 
      // Cstar est pas terrible mais on n'a pas besoin d'objets tres beau ici
      double cstar = star->Cstar();
      if ( mag < datSim.maxMag && mag > datSim.minMag
	   && cstar<0.45) 
	{
	  ++galcont;
	  galaxies.push_back(new SEStar(*star));
	}
    }

}

#include "nbrandom.h"
//long int random(void);
double my_rndm()
{
  return ((1.0*random())/(RAND_MAX+1.0));
}

static void Construct_SnList_WHost(SEStarList &galaxies, 
				   SEStarList & BellesEtoiles, 
				   SENearStarList &NearStarRef, 
				   double expo, double zerop, DatSim & datSim,
				   int XCCD, int YCCD)
{
    
  // construction de la liste de supernovae
  for (SEStarIterator it= galaxies.begin(); it != galaxies.end();++it)
    {
      SEStar *star = *it;
      
      //Caracteristique de la galaxie consideree
      double gal_X=star->X();
      double gal_Y=star->Y();
      double GalFlux = star ->flux;
      double gal_width = star->Fwhm();
      
      // Cherche l'etoile la plus proche dans la liste generee
      SEStar NearestStar = *(BellesEtoiles.FindClosest(gal_X,gal_Y));

      // Caracteristiques de cette etoile
      double StarFlux = NearestStar.flux;
      double NearStar_X = NearestStar.X();
      double NearStar_Y = NearestStar.Y();
      
      // Determination de la position de la Sn par rapport a la galaxie
      double NSig_X = 6.0*my_rndm() - 3.0;
      double NSig_Y = 6.0*my_rndm() - 3.0;
      
      // Calcul du deplacement de l'etoile future Sn
      int dx = (int) floor(gal_X - NearStar_X + NSig_X * gal_width); 
      int dy = (int) floor(gal_Y - NearStar_Y + NSig_Y * gal_width);
      
      // Calcul du facteur photometrique
      double sn_mag = datSim.minMag+(datSim.maxMag - datSim.minMag)*my_rndm();
      double sn_flux = expo*pow(10,(zerop - sn_mag)*0.4);
      double phot_factor = sn_flux / StarFlux;
      
      // Image size : to avoid adding stars outside of the image
      
      // Construction of the star list
      if ( NearestStar.X() + dx >0 &&  NearestStar.X() + dx < XCCD 
	   &&  NearestStar.Y() + dy >0 &&  NearestStar.Y() + dy < YCCD)
	{
	  NearStarRef.push_back(new SENearStar(NearestStar, dx, dy, 
					       phot_factor, GalFlux));
	}
    }
}

static void Construct_SnList_NoHost( SEStarList &BellesEtoiles, 
				     SEStarList &starList, 
				     SENearStarList &NearStarRef, 
				     double & expo, double & zerop, 
				     DatSim & datSim, 
				     const ReducedImage *AnImage, 
				     Image & mask)
{

  Image img = FitsImage(AnImage->FitsImageName(Calibrated));  
  int compteur=0;
  while (compteur<datSim.numberOfFakes)
    {
      for (SEStarIterator it= BellesEtoiles.begin(); it != BellesEtoiles.end();++it)
	{
	  if (compteur>datSim.numberOfFakes) break;
	  SEStar *Star = *it;
	  double StarFlux = Star->flux;
	  
	  // Calcul du deplacement aleatoire de la fakes dans un rayon shift
	  double shift = 400; // HC
	  int dx = (int) floor( ((2 * shift * my_rndm() - shift))); 
	  int dy = (int) floor( ((2 * shift * my_rndm() - shift)));
	  
	  // Calcul du rapport photometrique
	  double sn_mag = datSim.minMag+(datSim.maxMag - datSim.minMag) 
	    * my_rndm();
	  double sn_flux = expo*pow(10,(zerop - sn_mag)*0.4);
	  double phot_factor = sn_flux / StarFlux;
	  
	  // Calcul de la position aleatoire de le fake
	  double X_fakes = Star->X() + dx;
	  double Y_fakes = Star->Y() + dy;
	  
	  // recherche de l'etoile la plus proche dans l'image
	  SEStar NStar = *(starList.FindClosest(X_fakes, 
						Y_fakes));
	  // recherche de l'etoile la plus proche dans la liste de fakes
	  SEStar NNStar;
	  if (compteur > 0)
	    NNStar = *(NearStarRef.FindClosest(X_fakes,Y_fakes));
	  
	  // On ajoute l'etoile dans la liste de fakes si elle est 
	  // suffisament loin des autres fakes, des objets de l'image,
	  // des etoiles saturees, des pixels morts et si ell est 
	  // dans l'image

	  int XCCD = AnImage->XSize();
	  int YCCD = AnImage->YSize();

	  
	  if ( 0 < X_fakes && X_fakes < XCCD &&
	       0 < Y_fakes && X_fakes < YCCD &&
	       fabs(X_fakes-NStar.X()) > 15  &&
	       fabs(Y_fakes-NStar.Y()) > 15  &&
	       fabs(X_fakes-NNStar.X()) > 15 &&
	       fabs(Y_fakes-NNStar.Y()) > 15 &&
	       IsNotBadPix(img, mask, int(X_fakes),int(Y_fakes),10))
	    {
	      NearStarRef.push_back(new SENearStar(*Star, dx, dy, 
						   phot_factor));
	      ++compteur;
	    }
	} 
    }
}


void MakeListSn( const ReducedImage *AnImage, ReducedImage *New, SENearStarList &NearStarRef, bool AssociateGal=false) // HC
{

  double zerop = AnImage->Zerop();
  
  DatSim datSim(DefaultDatacards());
  

  double expo = AnImage->Exposure();
  
  if (zerop==0)
    {
      cout << "Zero USNO missing!!" << endl;
    }
  cout << " using Ravel ZEROP " << zerop << " as zero point with expo of " << expo << endl;
  
  
  int XCCD = AnImage->XSize();
  int YCCD = AnImage->YSize();
  
  // On colle les fakes jusqu'a atteindre le nombre attendu
  // on genere un mask pour ne pas y coller les sne  
  Image img = FitsImage(AnImage->FitsImageName(Calibrated));

  FitsImage *satur = new FitsImage(AnImage->FitsSaturName());

  if (FitsImage(AnImage->FitsDeadName()).IsValid())
    {
      FitsImage *dead = new FitsImage(AnImage->FitsDeadName());
            
      Image mask = (*dead)+(*satur);
      delete dead, satur;      
    }
  Image mask = (*satur);
  delete satur; 
  // Construction de la liste a l'aide de CvDetection
  
  CandidateStarList stl;     
  double seeing=AnImage->Seeing();
  DatDetec datdet(DefaultDatacards());
  NewCvDetection(img, mask, stl, datdet, seeing, seeing);
  SEStarList *starList = (SEStarList*) &stl;

#ifdef STORAGE  //
  SEStarList stl(AnImage->CatalogName());
  SEStarList *starList = (SEStarList*) &stl;
#endif
  // Construction de la liste d'etoiles
  SEStarList BellesEtoiles;
  Construct_Star(*starList, BellesEtoiles, datSim);
  
  // Construction de la liste de fakes associes a des galaxies
  if (AssociateGal)
    {
      // Construction de la liste de galaxies a partir de la new
      SEStarList galaxies;
      SEStarList GalList(New->CatalogName());
      double expoNew = New->Exposure();
      Construct_Gal(galaxies, GalList,  expoNew, zerop, datSim); 
      Construct_SnList_WHost(BellesEtoiles, galaxies, NearStarRef, 
			     expo, zerop, datSim, XCCD, YCCD);
      
    }
  
  // Construction de la liste de fakes sans association de galaxie
  if(!AssociateGal)
    {
      cout << "Processing random Fakes " << endl;
      

      Construct_SnList_NoHost(BellesEtoiles, *starList, NearStarRef, 
			      expo, zerop, datSim, AnImage, mask);
    }

}

