#include "sestar.h"
#include "starlist.h"
#include "image.h"
#include "myrdm.h"
#include "simsnstar.h" 
#include "transformedimage.h"



#include "reducedimage.h"
#include "fitsimage.h"
#include "toadscards.h"
#include "simulation.h"


// ================ DatSim ========================== 

DatSim::DatSim(const string &FileName, string const & filter)
{
  Read(FileName, filter);
}
void 
DatSim::Read(const string &FileName, string const & filter)
{
  if (FileExists(FileName))
    {
      DataCards cards(FileName);
      LitDataCards(cards,filter);
    }
  else
    {
      cerr << " DatSim::Read Cannot open  FileName " << endl;
    }
}



// Reads the parameters for the simulation
void DatSim::LitDataCards(DataCards &Data, string const & filter)
{
  Methode = Random;
  string line = Data.SParam("METHOD");
  if (strstr(line.c_str(),"FAKES_INHOST"))
    Methode = InHost;

  if (strstr(line.c_str(),"FAKES_ADAPTED"))
    Methode= AdaptedToHost; 

  if (strstr(line.c_str(),"FAKES_RANDOM"))
    Methode = Random; 

  if (strstr(line.c_str(),"FAKES_ONGRID"))
    Methode = Damier ; 

  numberOfFakes = Data.IParam("NUMBER_OF_FAKES");
  double i_minMag = Data.DParam("MIN_FAKE_MAG_I");
  double i_maxMag = Data.DParam("MAX_FAKE_MAG_I");
  double r_minMag = Data.DParam("MIN_FAKE_MAG_R");
  double r_maxMag = Data.DParam("MAX_FAKE_MAG_R");
  double g_minMag = Data.DParam("MIN_FAKE_MAG_G");
  double g_maxMag = Data.DParam("MAX_FAKE_MAG_G");
  xsize_damier = Data.IParam("XSIZE_GRID");
  ysize_damier = Data.IParam("YSIZE_GRID");
  delta_x = Data.IParam("FAKE_DELTA_X");
  delta_y = Data.IParam("FAKE_DELTA_Y");
  deltaMagmin= Data.DParam("MIN_DELTA_MAG");
  deltaMagmax= Data.DParam("MAX_DELTA_MAG");
  if (minMag > maxMag) swap (minMag,maxMag);
  if ( deltaMagmin > deltaMagmax) swap(deltaMagmin,deltaMagmax);

  // take all galaxies with mag < maxGalMag
  double i_maxGalMag=26 ;  
  double r_maxGalMag=26 ; 
  double g_maxGalMag=26 ;
  // all objects with mag > LimGalMag are identified as galaxies
  double i_LimGalMag = 24 ;
  double r_LimGalMag = 24.5 ;
  double g_LimGalMag = 24. ;

  if ( (filter == "g") || (filter == "G"))
    {
      cerr << "Choosing g filter specifications for MC. " << endl ;
      minMag = g_minMag;
      maxMag = g_maxMag;
      LimGalMag  = g_LimGalMag ;
      maxGalMag = g_maxGalMag ;
    }
  else
    if ( (filter == "r") || (filter == "R"))
      {
	cerr << "Choosing r filter specifications for MC. " << endl ;
	minMag = r_minMag;
	maxMag = r_maxMag;
	LimGalMag  = r_LimGalMag ;
	maxGalMag = r_maxGalMag ;
      }
    else
      if ( (filter == "i") || (filter == "I"))
	{
	  cerr << "Choosing i filter specifications for MC. " << endl ;
	  minMag = i_minMag;
	  maxMag = i_maxMag;
	  LimGalMag  = i_LimGalMag ;
	  maxGalMag = i_maxGalMag ;
	}
      else
	{
	  minMag = i_minMag;
	  maxMag = i_maxMag;
	  LimGalMag  = i_LimGalMag ;
	  maxGalMag = i_maxGalMag ;
	  cerr << "Unrecognized filter : " << filter << ", choosing i filter specification." << endl ;
	}
  

}
void DatSim::Print()
{
  cout << " ** DatSim ** " << endl ;
  cout << " Method: " ;
  if (Methode == Random) cout << "Random" ;
  if (Methode == InHost ) cout << "InHost " ;
  if (Methode == AdaptedToHost) cout << "AdaptedToHost " ;
  if (Methode == Damier) cout << "Damier " ;
  cout << endl ;
  cout << " FAKES : number  of : " <<  numberOfFakes 
      << " mag range : ["<< minMag << ',' << maxMag << ']' << endl;
  cout << " GALAXIES: maximum mag : " << maxGalMag << endl;
  cout << " object fainter than mag : " << LimGalMag << " are selected as galaxies. " 
       << endl ;     
  cout << "           mag diff. between gal. and SN: ["<< deltaMagmin << ',' 
       << deltaMagmax << ']'  << endl;
  cout << " In case of GRID generation : " << endl
       << "         Delta X, Delta Y : " << delta_x << ", " << delta_y << endl
       << "         X-Dim, Y-Dim : " << xsize_damier << ", " << ysize_damier << endl ;
}




// *******FAKE Generation, according DatSim specification ******

/*SEStar DatSim::Ran(double xsnref, double ysnref, double zerop_ref )
{

  // tirage en coordonnees ref.
  double snmag = minMag+(maxMag - minMag) 
    * my_rndm();
  double sn_flux = pow(10,(zerop_ref - snmag)*0.4);
  SEStar sestar(xsnref,ysnref,sn_flux);
  return(sestar) ;
  } */
SimSNStar*
DatSim::RandomFluxSN(double xsnref, double ysnref, double zerop_ref )
{

  // tirage en coordonnees ref.
  double snmag = minMag+(maxMag - minMag) 
    * my_rndm();
  double sn_flux = pow(10,(zerop_ref - snmag)*0.4);
  SimSNStar *p = new SimSNStar(xsnref,ysnref,sn_flux);
  p->Mag_SN() = snmag ;
  return(p) ;
 } 

SimSNStar*
DatSim::RandomSN(int Nx, int Ny, double zerop_ref )
{

  // tirage en coordonnees ref.
  double xsnref = Nx * my_rndm();//coord image ref
  double ysnref = Ny * my_rndm();//coord image ref
  SimSNStar *p = RandomFluxSN(xsnref, ysnref , zerop_ref) ;
  return(p) ;
 } 

static
SEStar *RandomAndRemoveGalaxy(int Nx, int Ny, SEStarList & BG )
{
  // tirage en coordonnees ref. --> tirage au sort galaxie
  double xref = Nx * my_rndm();//coord image ref
  double yref = Ny * my_rndm();//coord image ref
  SEStar *galref = BG.FindAndRemoveClosest(xref,yref);
  return(galref);
}
static void
MagRange(DatSim  const & datsim, double galmag, 
	 double & minsnmag, double & maxsnmag)
{
  minsnmag = galmag + datsim.deltaMagmin; 
  maxsnmag = galmag + datsim.deltaMagmax;
  cerr << "Magnitude range from galaxy mag: " << minsnmag << " " << maxsnmag << endl ;
  cerr << "Magnitude range from datacards: " <<datsim.minMag << " " << 
datsim.maxMag << endl ;
  if ( minsnmag < datsim.minMag)
    minsnmag = datsim.minMag;
  if ( maxsnmag > datsim.maxMag)
    maxsnmag = datsim.maxMag ;
  cerr << "Final Magnitude range: " << minsnmag << " " << maxsnmag << endl ;
}
 
  
SimSNStar*
DatSim::RandomSNInGalaxy(SEStar & gal_onref ,double zerop_ref )
{

  double GalMag = -2.5*log10(gal_onref.flux)+ zerop_ref;
  // tirage SN dans galaxy
  double DeltaX  = my_rndm()*2.* gal_onref.A() - gal_onref.A();
  double DeltaY  = my_rndm()*2.* gal_onref.A() - gal_onref.A();
  // sn coord on image ref
  double xsnref = gal_onref.X() + DeltaX;
  double ysnref = gal_onref.Y() + DeltaY;

  double minmag = minMag;
  double maxmag = maxMag;
  if(Methode ==  AdaptedToHost)
    MagRange(*this, GalMag, minmag , maxmag);
  double snmag = minmag+(maxmag - minmag)*my_rndm();
  double sn_flux = pow(10,(zerop_ref - snmag)*0.4);
  SimSNStar *p = new SimSNStar(xsnref,ysnref,sn_flux);
  p->Mag_SN() = snmag;
  p->Mag_Gal() = GalMag;
  p->X_Gal() = gal_onref.X();
  p->Y_Gal() = gal_onref.Y();
  p->A_Gal() = gal_onref.A();
  p->Fluxmax_Gal() = gal_onref.Fluxmax();
  return(p) ;
 } 











//  ********** Model Star Selection



static
void 
PreSelectModelStars(SEStarList const &starList, 
		    SEStarList &BellesEtoiles, Image * model)
{

  // Vire les objets a problemes 
  SEStarList stl;
  starList.CopyTo(stl);
  cout << "List size for stars selection: " << stl.size() << endl;
  RemoveNonOKObjects(stl);
 cout << "List size after good objetcs selection: " << stl.size() << endl;
  // On cherche les belles etoiles dans l'image
  SEStarList etoiles;
  StarFinder(stl,etoiles);
  // On vire les etoiles en dessous d'un certain signal sur bruit 
  cout << "List size after stars selection: " << etoiles.size() << endl;

  //etoiles.write("etoiles.list");

  int n = 0 ;
  BaseStarList *lb = ( BaseStarList *) &etoiles;
  for (SEStarIterator it= etoiles.begin(); it!= etoiles.end(); ++it )
    {

      SEStar *mystar = *it;
      double STNoise = mystar->flux / mystar->EFlux();
      bool no_neigbor = ! (mystar->HasBigCloseNeighbor(*lb, 1.5*mystar->StampSize(), 1.e-3)) ;
      bool checkflux = true ;
      if ( model != NULL )
	{
	  double sumflux = model->SumPixels(mystar->x, mystar->y, mystar->StampSize());
	  double error = abs( (mystar->flux - sumflux) / mystar->flux ) ;
	  // ok if | error | < 0.2 on model image, supposed to be very clean. 
	  checkflux = ( error < 0.2 ) ;
	}
      if ( checkflux && no_neigbor && (STNoise > 75 )) // convient pour SNLS
	{
	  BellesEtoiles.push_back(new SEStar(*mystar));
	  n++ ;
	}
    } 
  cout << "After S/N cut, number of model stars " 
       << BellesEtoiles.size() << endl;
}

// short utilitaries used in PrefIll Constructor
static
void 
SelectionOfModelStars(SEStarList const &starList,
		 SEStarList &BellesEtoiles, Image * model)
{
  unsigned int MaxNumber = 100 ;
  PreSelectModelStars(starList, BellesEtoiles, model);
  BellesEtoiles.FluxSort(); // first of the list has the biggest flux

  if (BellesEtoiles.size() > MaxNumber) 
    BellesEtoiles.CutTail(MaxNumber); // keeping MaxNumber stars.


  cerr << "Final number of model stars (datacard : " 
       <<  MaxNumber << ") : " << BellesEtoiles.size() << endl;
}


// removes in stl the stars that are not present
// or not OK in the catalogs of a ReducedImageList.
// all  the images and stl are supposed to be ALIGNED.

static
bool SelectStarsInList(ReducedImageList & imglist, SEStarList & stl)
{
  for (ReducedImageIterator c = imglist.begin();
       c != imglist.end(); ++c)
    {
      
      ReducedImage *ti = *c ;
      if (! ti->HasCatalog())
	ti->MakeCatalog();
      if (! ti->HasCatalog())
	{
	  cerr << " Failling producing catalog for " << ti->Name() 
	       << ": exiting " << endl ;
	  return false;
	}
      SEStarList stli(ti->CatalogName());
      cerr << "Reading: " << stli.size() << " stars in " << ti->CatalogName() << endl ;
      double min_dist = 4. ;// no need to be too tough
      // as the distances are big
      min_dist *= min_dist;
      int n = 0, nout=0 ;
      if (! ti->HasWeight() )
	ti->MakeWeight();
      if (! ti->HasWeight())
	{
	  cerr << " Failling producing weight image for " << ti->Name() 
	       << ": exiting " << endl ;
	  return false;
	}	
      FitsImage weight(ti->FitsWeightName());
      for (SEStarIterator it = stl.begin(); it !=  stl.end() ; )
	{
	  SEStar *Stari = stli.FindClosest((*it)->x,(*it)->y ) ;
	  double dist2 = Stari->Dist2(*(*it));
	  if ( dist2 < min_dist) // ref star recovered in list i
	    if ( !(Stari->IsOK(weight))  )
	      {
		it = stl.erase(it); n++ ;
	      }
	    else
	      it++;
	  else
	    {
	      cerr << "Star not present: " << (*it)->x << " " << (*it)->y << endl ;
	     it = stl.erase(it); nout++ ;
	    } 
	}
      cerr << " Removing from StarList " << n << " objects that are bad on " << ti->Name() << endl ;
      cerr << " Removing from StarList " << nout << " objects that are not present on " << ti->Name() << endl ;
    }
  return true;
}



bool 
SelectModelStars(ReducedImageList & imglist, SEStarList const &stlref,
		 SEStarList &BellesEtoiles, Image * model)
{
  PreSelectModelStars(stlref,BellesEtoiles,model);
  cerr << "Number of selected stars : " << BellesEtoiles.size() << endl ;
  bool ok = SelectStarsInList(imglist, BellesEtoiles);
  if (!ok)  return false ;
  cerr << "Final Number of selected stars : " << BellesEtoiles.size() << endl ;

  return( (BellesEtoiles.size() > 0) ) ;
}








// ********** selection de galaxies hotes.**********


static bool
GalIsOK(double galmag, DatSim const & datsim)
{
  double minsnmag = galmag + datsim.deltaMagmin; 
  double maxsnmag = galmag + datsim.deltaMagmax;
  if ( minsnmag > datsim.maxMag)
    return false;
  if ( maxsnmag < datsim.minMag)
    return false;
  return true ;
}
 
static void 
SelectHostGalaxies(DatSim  const & datsim, 
	      SEStarList &galaxies, SEStarList & liste, 
	      double zerop_gal)
{
  cout << "Selecting host galaxies " << endl;
  RemoveNonOKObjects(liste); 
  SEStarList lgal;
  // all objects FAINTER than datsim.LimGalMag are galaxies.
  GalaxyFinder(liste, lgal, datsim.LimGalMag, zerop_gal);
  cout << "Selecting " << lgal.size() << " objects as galaxies." <<  endl;
  for (SEStarIterator it= lgal.begin(); it != lgal.end();++it ) 
    { 
      SEStar *star = *it;
      double GalMag = -2.5*log10(star->flux)+zerop_gal;
      if ( GalMag< datsim.maxGalMag )
	{
	  bool ok = true ;
	  if (datsim.Methode ==  AdaptedToHost) 
	    ok = GalIsOK(GalMag,datsim) ; 
	  if (ok)
	    galaxies.push_back(new SEStar(*star));
	}
    }
   
  cout << "After magnitude cut: " << galaxies.size()  << endl;
}

// ================ ForSim ========================== 


// une routine static de debug
//verification du seeing des etoiles modeles.
static void
CheckSeeing(SEStarList const &BellesEtoiles, 
	    ReducedImage const & ForModelStars)
{
  double S = 0., S2=0. ; 
  for (SEStarCIterator it= BellesEtoiles.begin(); it != BellesEtoiles.end(); it++)
    {
      double fwhm = (*it)->Fwhm();
      S += fwhm ;
      S2 += fwhm*fwhm ;
    }
  S  = S/BellesEtoiles.size()/(2.*sqrt(2.*log(2.)));
  S2 = S2/BellesEtoiles.size()/(8.*log(2.));
  cerr << "Seeing computed on model stars: " << S << " +/- " << sqrt(S2-S*S) << endl ;
  cerr << "Seeing as indicated on model image: " << ForModelStars.Seeing() << endl ;
}



// une routine static de debug
//sauvegarde d'une image a zero sauf autour des objets de la stl:
// pour verifier d'un coup d'eil les etoile smodeles ou les galaxies hotes.
#include "fitsimage.h"

static
void SaveModel(SEStarList & stl, 
	       ReducedImage const & r, string nom )
{

 string catmodel = r.Dir()+ "/" + nom + ".list" ;
 stl.write(catmodel);
 cerr << "Writing " << catmodel << endl ; 
 string name = r.Dir()+ "/" + nom + ".fits" ;
 // test on the existence
 // so that when simulations are called in sequence, the
 // model stars/galaxies are not saved each time. 
 if ( ! FileExists(name) )
   {
     Image test(r.XSize(),r.YSize() );
     {
       FitsImage ref(r.FitsName());
       GtransfoIdentity identity;
       for (SEStarCIterator c = stl.begin() ; c != stl.end() ; c++)
	 {
	   ModelStar model(*(*c), (*c)->StampSize(),0, 0,1.) ;
	   model.AddToImage(ref, test,(Gtransfo *) &identity);
	 }
     }
     
     FitsHeader head(r.FitsName());
     FitsImage tt(name,head,test);
   }
}  

// ****************** ForSim ***************

void
ForSim::Fill(ReducedImage const & RefImage)
{
 
  zero_point = RefImage.AnyZeroPoint() ;
  cerr << "Image de selection des Modeles :" << RefImage.Name() 
       << " Zero Point : " << zero_point << endl ;
  XCCD = 0 ;
  YCCD = 0 ;
  XrefCCD = RefImage.XSize();
  YrefCCD = RefImage.YSize(); // on suppose ref et image i meme dim.
  if (datsim.Methode == Damier)
    {
      if (datsim.xsize_damier > 0)
	XCCD = datsim.xsize_damier;
      else
	XCCD = RefImage.XSize();
      if (datsim.ysize_damier > 0)
	YCCD = datsim.ysize_damier; 
      else
	  YCCD = RefImage.YSize();
      cout << "Working on grid " << XCCD << " x " 
	   << YCCD << endl;
    }
 
 
  if (datsim.Methode == InHost || datsim.Methode == AdaptedToHost)
    {
      SEStarList liste(RefImage.CatalogName());
      cerr << "Image de selection des Galaxies :" << RefImage.Name() 
       << " Zero Point : " << zero_point << endl ;
      SelectHostGalaxies(datsim, BellesGal, liste, zero_point );
      SaveModel(BellesGal, RefImage, "galaxies" );
    }
}

ForSim::ForSim(ReducedImage const & RefImage)
{
  string band = RefImage.Band();
  datsim.Read(DefaultDatacards(), band);
  datsim.Print();
  Fill(RefImage);
}






// ****************** ForSimWModel ***************

// partie commune aux 2 createurs adaptes a mon utilisation
// personnelle (DH)
// prend une DatSim et la recopie, recupere taille ref, point zero,
// son catalogue pour la selection des hotes 


void
ForSimWModel::Fill(ReducedImage const & RefImage,
		   SEStarList const & ListForModelStars, 
		   const ImageGtransfo* tf)
{
   ForSim::Fill(RefImage);

  ListForModelStars.CopyTo(BellesEtoiles);
  CheckSeeing(BellesEtoiles, RefImage);
  SaveModel(BellesEtoiles, RefImage, "model");

  Transfo = tf->TransfoFromRef();
  TransfoInv = tf->TransfoToRef() ;
  cerr << "Transfos For FAKE computing: " << endl ;
  Transfo->dump();
  TransfoInv->dump();
}


ForSimWModel::ForSimWModel(ReducedImage const & RefImage,
	       SEStarList const & ListForModelStars, 
	       const ImageGtransfo* tf)
{
 
  // datacard lue dans le default
  string band = RefImage.Band();
  datsim.Read(DefaultDatacards(), band);
  datsim.Print();
  Fill(RefImage,ListForModelStars,tf);
}

void 
ForSimWModel::Construct_Sn_Random(SimSNWModelStarList & SNList)
{
  // imprimer les 5 premieres etoiles modeles pour debug .......
  int nnn = 0 ;
  cerr << "This are the first 5 ModelStar " << endl ;
  for (SEStarCIterator it= BellesEtoiles.begin(); it != BellesEtoiles.end() && nnn < 5 ; it++)
    {
      double x1,y1 ;
      Transfo->apply((*it)->x,(*it)->y,x1, y1);
      cerr << (*it)->x << " " << (*it)->y << " "
	   <<x1 << " "<< y1 << endl ;
      nnn++;
    }
 
  int compteur=0 ;
  while (compteur<datsim.numberOfFakes)
    {
      SimSNStar *psim = datsim.RandomSN(XrefCCD, YrefCCD,zero_point);
      SimSNWModelStar *p = new SimSNWModelStar(*psim);
      delete psim ;
      p->MariageToAModelStar(BellesEtoiles,Transfo,TransfoInv);
      // apres mariage, les coordonnees de la fake ont peut etre un 
      // peu bouge.  verif que est toujours dans la ref.
      if ( 0 < p->x  && p->x  < XrefCCD &&
	   0 < p->y  && p->y < YrefCCD)
	{
	  SNList.push_back(p);
	  ++compteur;
	}
      else
	delete p ;
    } 
}





//** position en Damier
void 
ForSimWModel::Construct_Sn_Damier(SimSNWModelStarList & SNList)
{

  // position en damier sur image 1 et non ref !
  int xsn1 = datsim.delta_x ;
  int ysn1 = datsim.delta_y ;
  while (ysn1 < YCCD )
    {

      
      // coord ref sn.
      double xsnref, ysnref ;
      TransfoInv->apply(xsn1, ysn1,xsnref,ysnref);
      SimSNStar *psim = datsim.RandomFluxSN(xsnref, ysnref, zero_point);
      SimSNWModelStar *p = new SimSNWModelStar(*psim);
      delete psim ;
      p->MariageToAModelStar(BellesEtoiles,Transfo,TransfoInv);
      SNList.push_back(p);
	
      xsn1 += datsim.delta_x ;
      if (xsn1 > XCCD )
	{xsn1 = datsim.delta_x ; ysn1 += datsim.delta_y ;}
    }
  cerr << SNList.size() << " stars generated " << endl ;
}

//** Dans Hote.


void 
ForSimWModel::Construct_SnList_WHost(SimSNWModelStarList & SNList)
{
  SEStarList BG;
  (BellesGal).CopyTo(BG); //coord as on image ref.
  if (BG.size() == 0)
    {
      cerr << "0 Model galaxies selected, no SNe generation" << endl ;
      exit(1);
    }
 
  
  // construction de la liste de supernovae ds gal
  // on tire une galaxie au hasard,on l'enleve de la liste
  // (pas 2 sn ds la meme galaxie !) puis on tire
  // 1 sn au hasard
  int NSN= 0 ;
  while (NSN <datsim.numberOfFakes && BG.size()>0 )
    {
      SEStar *galref = RandomAndRemoveGalaxy(XrefCCD,YrefCCD , BG ) ;
      SimSNStar *psim = datsim.RandomSNInGalaxy(*galref , zero_point);
      SimSNWModelStar *p = new SimSNWModelStar(*psim);
      delete psim ;
      p->MariageToAModelStar(BellesEtoiles, Transfo,TransfoInv);
      SNList.push_back(p);
      // a cause du mariage, la position peut bouger.  verif que est toujours dans la ref.
      if ( 0 < p->x  && p->x  < XrefCCD &&
	   0 < p->y  && p->y < YrefCCD)
	{
	  SNList.push_back(p);
	  ++NSN;
	}
      else
	delete p ;
      delete galref;
      NSN++;
    }
  cerr << NSN<< " supernovae generated in host." << endl ;
}

// ********** Enrobage de tout ca. **********




  
void 
ForSimWModel::MakeListSN(SimSNWModelStarList & SNList)
{    
  switch (datsim.Methode)
    { 
    case  Random :
      cout << "Generating random Fakes list." << endl;
      Construct_Sn_Random(SNList);
      break;
    case Damier :
      cout << "Generating Fakes list on grid " << XCCD << " x " << YCCD << endl;
      Construct_Sn_Damier(SNList);
      break;
    case InHost : 
    case AdaptedToHost :  
      cout << "Generating Fakes in Host list." << endl;
      Construct_SnList_WHost(SNList);
      break;
    default:
      cerr << "Unknown Method " << datsim.Methode  << " for Fake SN Generation " << endl ;
    }
    

}


