#include "simsnstar.h"
#include "image.h"
#include "gtransfo.h"
#include "agaussian.h"
#include "daophotpsf.h"
#include "myrdm.h"
#include <fstream>



/********************* Model  Star ********************/

ModelStar::ModelStar()
: BaseStar(0.,0.,0.)
{
  
  stampsize = 0. ;
  xShift =  0 ;
  yShift = 0 ;
  photFactor = 0. ;
}


ModelStar::ModelStar(const BaseStar &AStar, const double size, const int Dx, const int Dy, const double PhotFactor) 
  : BaseStar(AStar)
{
  stampsize = size;
  xShift =  Dx;
  yShift = Dy;
  photFactor = PhotFactor;

}

string ModelStar::WriteHeader_(ostream &pr, const char *i) const
{
  if (i == NULL) i = "";
  string BaseStarFormat = BaseStar::WriteHeader_(pr,i);
  pr << "# dx"<< i << " : offset to real position" << endl;
  pr << "# dy"<< i << " : offset to real position" << endl;
  pr << "# size"<< i << " : size of the stamp to be added" << endl;
  pr << "# photfac"<< i << " : phot factor from model to fake. " << endl;
  return BaseStarFormat +" ModelStar 1 ";
}


void
ModelStar::dumpn(ostream& s) const
{
 BaseStar::dumpn(s);
 s << " xShift : "   << xShift
   << " yShift : "   << yShift 
   << " stampsize : "   << stampsize
   << " photFactor : "   << photFactor ;
}

void
ModelStar::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}
void ModelStar::writen(ostream & pr) const
{
  BaseStar::writen(pr);
  pr << ' ' << xShift << ' ' << yShift  << ' ' << stampsize << ' ' <<  photFactor ;
}

/*
void ModelStar::read_it(istream& r, const char* Format)
{
  BaseStar::read_it(r, Format);
  r >>  xShift ;
  r >> yShift  ;
  r >> stampsize ;
  r >> photFactor  ;
}
ModelStar*  ModelStar::read(istream& r, const char *Format)
{
  ModelStar *pstar = new ModelStar();  
  pstar->read_it(r, Format);
  return(pstar);
}
*/


static void 
AddModelToImage(double xmod, double ymod, int xShift , int yShift,
		double photFactor, double stampsize, const Image &image, 
		Image & dest,  Image * psat, double satlevel, double gain) 
{
  int flagsaturated = 0 ;
  int xstart = int(xmod - stampsize + 0.5);
  int ystart = int(ymod - stampsize + 0.5);

  int taille =  int(2*stampsize + 1. + 0.5);

  int xend = min(xstart + taille,image.Nx()) ;
  int yend = min(ystart + taille,image.Ny()) ;

  xstart = max(xstart,0);
  ystart = max(ystart,0);

 
  double FFlux=0. ;
  int Aire = 0 ;  
  double  skylevel = 0. ;
  for (int j = ystart; j < yend; ++j)
    {
      int targetj = j+ yShift;
      if (targetj < 0) continue;
      if (targetj >= dest.Ny()) break;
      for (int i = xstart; i < xend; ++i)
	{
          int targeti = i + xShift;
	  if (targeti < 0 ) continue;
          if (targeti >= dest.Nx()) break;
	  
	  // the CCD gain is in units of (e-/ADU)
	  double val = ((double)image(i,j)- skylevel)*photFactor;
	  // now set a random value for this image
	  val = Poisson(val*gain)/gain;	  
	  
	  // total flux computed for debug.
	  FFlux += val ;

	  if ( (satlevel > 0 ) && 
	       (dest(targeti , targetj)+ val > satlevel ) )
	    {
	      //cerr << "dest: " << targeti << " " << targetj << " " << satlevel << endl ;
	      dest(targeti , targetj) = satlevel ;
	      if (psat != NULL)
		{
		  if (flagsaturated == 0 )
		    if ((*psat)(targeti , targetj) == 0 ) // not saturated before, saturated because of SNe
		      {
			cerr << "SN is saturated (i= " 
			     << targeti << ", j= " 
			     << targetj << "), saturation level = " 
			     << satlevel << ", image = " 
			     << dest(targeti , targetj) 
			     << ", new image val = " 
			     << dest(targeti , targetj)+ val  
			     << " : " ; 
			flagsaturated = 1;
		      }	
		    
		  (*psat)(targeti , targetj) = 1 ;
		}
	      else
		if (flagsaturated == 0 )
		  {
		    cerr << "SN or SN zone  is saturated: " ; 
		    flagsaturated = 1;
		  }
	      
	    }
	  else
	    {
	      dest(targeti , targetj) += val ;
	      //cerr << "dest: " << targeti << " " << targetj << " " << dest(targeti , targetj) << endl ;
	    }
	  // area computed for debug.
	  Aire++;
	}
    }
  
  //cout << "in AddModelToImage: photFactor, FFLux, satlevel = " << photFactor << " " << FFlux << " " << satlevel << endl;
  
} 


void ModelStar::AddToImage(const Image &image, Image & dest, 
			   const Gtransfo *Transfo, Image * psat,
			   double satlevel, double gain) const
{
  double x_model, y_model ;
  Transfo->apply(x,y,x_model, y_model);
  AddModelToImage(x_model, y_model, xShift , 
		  yShift, photFactor,	
		  stampsize, image, dest, psat, satlevel, gain);
}




/* ******************** Simulated SN Star ******************* */



SimSNStar::SimSNStar()
: BaseStar(0.,0.,0.)
{
  mag_sn = 0.  ;
  mag_gal= 0.  ;
  x_gal= 0.  ;
  y_gal= 0.  ;
  a_gal= 0.  ;
  fluxmax_gal = 0.;
}

SimSNStar::SimSNStar(double xx, double yy, double ff)
: BaseStar(xx,yy,ff)
{
  mag_sn = 0.  ;
  mag_gal= 0.  ;
  x_gal= 0.  ;
  y_gal= 0.  ;
  a_gal= 0.  ;
  fluxmax_gal = 0.;
}

void
SimSNStar::dumpn(ostream& s) const
{
 BaseStar::dumpn(s);
 s << " mag_sn : "   << mag_sn
   << " mag_gal : " << mag_gal
   << " x_gal : " << x_gal
   << " y_gal : " << y_gal
   << " a_gal : " << a_gal
   << " fluxmax_gal : " << fluxmax_gal ;
}

void
SimSNStar::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}

string SimSNStar::WriteHeader_(ostream &pr, const char *i) const
{
  if (i == NULL) i = "";
  string BaseStarFormat = BaseStar::WriteHeader_(pr,i);
  pr << "# snmag"<< i << " : generated magnitude" << endl;
  pr << "# galmag"<< i << " : galaxy magnitude " << endl; 
  pr << "# x_gal"<< i << " : host pos on ref " << endl;   
  pr << "# y_gal"<< i << " : host pos on ref" << endl; 
  pr << "# a_gal"<< i << " : galaxie size in pix." << endl; 
  pr << "# fmax_gal"<< i << " : galaxie pix. max value on ref" << endl;
  return BaseStarFormat + " SimSNStar 1 ";
}


void SimSNStar::writen(ostream & pr) const
{
  BaseStar::writen(pr);
  pr << ' ' << mag_sn << ' ' << mag_gal << ' ' << x_gal<< ' ' << y_gal 
     << ' ' << a_gal << ' ' << fluxmax_gal ;
}


/*
void SimSNStar::read_it(istream& r, const char* Format)
{
  BaseStar::read_it(r, Format);
  r >> mag_sn  ;
  r >> mag_gal  ;
  r >> x_gal  ;
  r >> y_gal  ;
  r >> a_gal  ;
  r >>  fluxmax_gal ;
  return;
}


SimSNStar*  SimSNStar::read(istream& r, const char *Format)
{
  SimSNStar *pstar = new SimSNStar();  
  pstar->read_it(r, Format);
  return(pstar);
}
*/


void 
SimSNStarList::NewFlux(double NewZeroPoint)
{
  for (SimSNStarIterator i = begin(); i != end(); ++i)
    {
      (*i)->NewFlux(NewZeroPoint) ;
    }
}



void 
SimSNStarList::AddWGaussianToImage(double sigmax, double sigmay, double rho,
				   Image & dest, 
				   const Gtransfo *Transfo,   
				   Image * psat, 
				   double satlevel) const 
{  
  AGaussian gauss;
  gauss.sigma_x = sigmax ;
  gauss.sigma_y = sigmay ;
  gauss.rho = rho  ;
  for (SimSNStarCIterator i = begin(); i != end(); ++i)
    {
      const SimSNStar *s = *i;
      // position of SN on image i
      double Xsn, Ysn ; 
      Transfo->apply(s->x,s->y,Xsn, Ysn );
      gauss.xc =  Xsn ;
      gauss.yc =  Ysn ;
      gauss.flux = s->flux ;
      gauss.AddToImage(dest, psat, satlevel);
    }
}

void 
SimSNStarList::AddWGaussianToImage(double sigmax, double sigmay, double rho,
				   Image & dest,    
				   Image * psat, 
				   double satlevel) const 
{  
  BaseStarList *list = (BaseStarList *) this ;
  AddListWGaussianToImage(sigmax, sigmay, rho, list,
			  dest, psat, satlevel);
}

// Converter :
BaseStarList* SimSN2Base(SimSNStarList * This)
{ return (BaseStarList*) This;} 

const BaseStarList* SimSN2Base(const SimSNStarList * This)
{ return (BaseStarList*) This;} 


#include "starlist.cc" /* since starlist is a template class */
template void StarList<SimSNStar>::CopyTo(StarList<SimSNStar>&) const;  /* to force instanciation */


/*     ************** SimSNWModelStar ************** */

string SimSNWModelStar::WriteHeader_(ostream &pr, const char *i) const
{
  if (i == NULL) i = "";
  string SimSNFormat =  SimSNStar::WriteHeader_(pr,i);
  string ii = i ;
  string uu = "mod" + ii ;
  string ModelStarFormat = model_on_ref.WriteHeader_(pr,uu.c_str());
  return SimSNFormat + ModelStarFormat + " SimSNWModelStar 1 ";
  }

void SimSNWModelStar::dumpn(ostream & pr) const
{
  SimSNStar::dumpn(pr);
  pr << " Modele : " ;
  model_on_ref.dumpn(pr);
}
void
SimSNWModelStar::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}

void SimSNWModelStar::writen(ostream & pr) const
{
  SimSNStar::writen(pr);
  pr << " " ;
  model_on_ref.writen(pr);
}

/*
void SimSNWModelStar::read_it(istream& r, const char* Format)
{
  SimSNStar::read_it(r, Format);
  model_on_ref.read_it(r, Format);
  return;
}

SimSNWModelStar*  SimSNWModelStar::read(istream& r, const char *Format)
{
  SimSNWModelStar *pstar = new SimSNWModelStar();  
  pstar->read_it(r, Format);
  return(pstar);
}
*/


int integer_delta(double xsn, double xmodel)
{
  int dx = 0 ;
  double delta = xsn - xmodel ;
  if (delta > 0 )
    dx = (int) (delta +0.5) ;
  else
    dx = (int) (delta -0.5) ;
  return(dx);

}



void 
SimSNWModelStar::MariageToAModelStar(SEStarList const & BellesEtoiles, 
			       const Gtransfo *Transfo, const Gtransfo *TransfoInv)
{
  double xsn1, ysn1; //coord sn image 1.
  Transfo->apply(x,y,xsn1, ysn1);

  // nearest model star in coord image ref
  SEStar model = *(BellesEtoiles.FindClosest(x,y)) ;
  // its  coordinates on image 1
  double xmodel1, ymodel1;
  Transfo->apply(model.X(), model.Y(),xmodel1, ymodel1 );

  // calcul dx sur image1
  int dx = integer_delta(xsn1, xmodel1);
  int dy = integer_delta(ysn1, ymodel1);
	  
      
  // Calcul du rapport photometrique
  double phot_factor = flux / model.flux;
	  
  // actual position fake sur image 1
  double XSN1 = xmodel1 + dx;
  double YSN1 = ymodel1 + dy;
  // position fake sur image ref (ne doit pas etre tres differente): 
  double XSN_ref, YSN_ref ;
  TransfoInv->apply(XSN1,YSN1,XSN_ref,YSN_ref);

  x = XSN_ref ;
  y = YSN_ref ;
  ModelStar modelstar(model,model.StampSize(),dx,dy,phot_factor);
  model_on_ref = modelstar ;
}

void 
SimSNWModelStar::MariageToAModelStar(SEStarList const & BellesEtoiles)
{
  
  GtransfoIdentity id;
  MariageToAModelStar(BellesEtoiles, &id, &id);
}








//static routine for debug, to check if the model star method is sound.
static string dirdebug = "mc";
static string debug= dirdebug + "/addmodel.debug";
// prdebugstate is incremented when sticking a whole list of fake --> corresponds to an image number.
static int prdebugstate=0;
static void PrintDebugHead(ofstream & pr)
{
  pr << "#n : image number" << endl 
     <<"#xsn_r : x sn on ref " << endl
     <<"#ysn_r : y sn on ref " << endl
     <<"#xmod_r : x model on ref" << endl
     <<"#ymod_r : y model on ref" << endl
     <<"#xsn_th : x sn on image computed with tf" << endl
     <<"#ysn_th : y sn on image computed with tf " << endl
     <<"#xmod : x model on image computed with tf" << endl
     <<"#ymod : y model on image computed with tf" << endl
     <<"#xsn : x fake = xmod+dx" << endl
     <<"#ysn : y fake = xmod+dy" << endl
     <<"#dx : " << endl
     <<"#dy : " << endl
     <<"#fracx : dif de xsn_th -xmod " << endl
     <<"#fracy : " << endl
     << "# end"  << endl ;

} 


#include "fileutils.h"
static void PrepareDebug()
{
  if (prdebugstate == 0)
    {
      MKDir(dirdebug.c_str());
      remove(debug.c_str());
      ofstream pr(debug.c_str());
      PrintDebugHead(pr);
      pr.close();
      prdebugstate=1;
    }
  return ;
}
static
void DebugModelMethod(const SimSNWModelStar  & sim, const Gtransfo *Transfo)
{
  // position of SN on image i
  double Xsn, Ysn ; 
  Transfo->apply(sim.x,sim.y,Xsn, Ysn );  
  // position of Model Star on image i
  double Xmodel, Ymodel ;
  Transfo->apply((sim.model_on_ref).x,(sim.model_on_ref).y,Xmodel, Ymodel);
  // on etudie la fraction entiere de Xsn - Xmodel, qui doit etre la plus petite possible.
  double fracx = Xsn - Xmodel - integer_delta(Xsn , Xmodel ) ;
  double fracy = Ysn - Ymodel - integer_delta(Ysn , Ymodel ) ;
  // print de debug
  if (prdebugstate > 0) // le print a deja commence, on append
    {
      ofstream prdebug(debug.c_str(), ios::app);
      prdebug << prdebugstate << " "
	      << sim.x << " " 
	      << sim.y << " "
	      << (sim.model_on_ref).x << " " 
	      << (sim.model_on_ref).y << " "
	      << Xsn << " " 
	      << Ysn << " " 
	      << Xmodel << " " 
	      << Ymodel << " " 
	      << Xmodel + (sim.model_on_ref).XShift() << " "
	      << Xmodel + (sim.model_on_ref).YShift() << " " 
	      << (sim.model_on_ref).XShift() << " "
	      << (sim.model_on_ref).YShift() << " " 
	      << fracx << " " 
	      << fracy << " " << endl ;
      prdebug.close();
    }

}



void 
SimSNWModelStarList::NewFlux(double NewZeroPoint)
{
  for (SimSNWModelStarIterator i = begin(); i != end(); ++i)
    {
      (*i)->NewFlux(NewZeroPoint) ;
    }
}

void 
SimSNWModelStarList::AddWModelToImage(const Image &image, Image & dest,  
				const Gtransfo *Transfo,  
				Image * psat,
				double satlevel, bool print_debug) const
{
  if (print_debug) PrepareDebug();
  for (SimSNWModelStarCIterator i = begin(); i != end(); ++i)
    {
      const SimSNWModelStar *s = *i;
      if (print_debug) DebugModelMethod(*s,Transfo);
      s->AddWModelToImage(image, dest, Transfo, psat, satlevel);
    }
  if (print_debug) prdebugstate += 1 ;
}


 



 




// Converter :
BaseStarList* SimSNWModel2Base(SimSNWModelStarList * This)
{ return (BaseStarList*) This;}

const BaseStarList* SimSNWModel2Base(const SimSNWModelStarList * This)
{ return (BaseStarList*) This;} 

#include "starlist.cc" /* since starlist is a template class */
//template class StarList<SimSNWModelStar>;  /* to force instanciation */

/***************** Utilitaires Dao **************/
//The place for this is to be discussed with Seb



// a comparer a
// void Daophot::Addstar(const string& PsfFileName, const string& AddFileName, const string& AddPicName, const int InSeed)
// a mettre ds daophoutils ?
void AddWDaoPsfToImage(DaoPsf const & daopsf, double xc, double yc, 
		       double flux,Image & image,
		       Image * psat,
		       double saturation)
{
  double dpdx, dpdy ;
  int NN = (int) (daopsf.Radius() + 1.) ; // au-dela de Radius, la 
  // DaoPsf.Value retourne 0.
  // cerr <<"taille fenetre " <<  NN << endl ;
 int flagsaturated = 0 ;
  for (int j = -NN; j <=NN ; j++)
    {
      int jj = (int) (j +yc + 0.5) ;
      if ( jj < 0 ) continue ;
      if ( jj >= image.Ny() ) break  ;
      for (int i = -NN; i <=NN ; i++)
	{
	  int ii = (int) (i +xc + 0.5) ;
	  if ( ii < 0 ) continue ;
	  if ( ii >= image.Nx() ) break  ;
	  double val = daopsf.Value(ii,jj,xc,yc,dpdx,dpdy) * flux;
	  if ( (saturation > 0 ) && 
	       (image(i,j)+ val > saturation ) )
	    {
	      image(ii,jj) = saturation ;
	      if (psat != NULL)
		{
		  if (flagsaturated == 0 ) // la 1ere fois qu'on s'en apercoit
		    if ((*psat)(ii,jj) == 0 ) // not saturated before, saturated because of SNe
		      {
			cerr << "SN is saturated (i= " 
			     << ii << ", j= " 
			     << jj << "), saturation = " 
			     << saturation << ", image = " 
			     << image(ii,jj) << ", new image val = " 
			     << image(ii,jj)+ val  << endl ;
			flagsaturated = 1; // pour ne pas le dire 36 fois
		      }	
		    
		  (*psat)(ii,jj) = 1 ;
		}
	      else
		if (flagsaturated == 0 )
		  {
		    cerr << "SN or SN zone  is saturated: " << endl ;
		    flagsaturated = 1;
		  }	      
	    }
	  else
	    image(ii,jj) += val ; // il faudrait faire un tirage poissonien 
	  // sur val.
	}
    }
} 

void AddListWDaoPsfToImage(DaoPsf const & daopsf, BaseStarList *List,
		       Image & img,  Image * psat,
		       double saturation)
{
  for (BaseStarIterator it= List->begin(); it!= List->end(); ++it )
    {     
      AddWDaoPsfToImage(daopsf, (*it)->x, (*it)->y,  (*it)->flux,img,psat,saturation);
    }
}

#include "imagepsf.h"
void AddWPsfToImage(ImagePSF &psf ,double xc, double yc, 
		       double flux,Image & image,
		       Image * psat,
		       double saturation)
{
int xmin,xmax,ymin,ymax;
psf.StampLimits( xc,yc,xmin,xmax,ymin,ymax);

 int flagsaturated = 0 ;
  for ( int j =  ymin ;j < ymax ; j++)
    {

      for (int i = xmin; i < xmax ; i++)
	{

	  double val = psf.PSFValue(xc,yc,i,j) * flux;
	  if ( (saturation > 0 ) && 
	       (image(i,j)+ val > saturation ) )
	    {
	      image(i,j) = saturation ;
	      if (psat != NULL)
		{
		  if (flagsaturated == 0 ) // la 1ere fois qu'on s'en apercoit
		    if ((*psat)(i,j) == 0 ) // not saturated before, saturated because of SNe
		      {
			cerr << "SN is saturated (i= " 
			     << i << ", j= " 
			     << j << "), saturation = " 
			     << saturation << ", image = " 
			     << image(i,j) << ", new image val = " 
			     << image(i,j)+ val  << endl ;
			flagsaturated = 1; // pour ne pas le dire 36 fois
		      }	
		    
		  (*psat)(i,j) = 1 ;
		}
	      else
		if (flagsaturated == 0 )
		  {
		    cerr << "SN or SN zone  is saturated: " << endl ;
		    flagsaturated = 1;
		  }	      
	    }
	  else
	    image(i,j) += val ; // il faudrait faire un tirage poissonien 
	  // sur val.
	}
    }
} 

void AddListWPsfToImage(ImagePSF &psf, BaseStarList *List,
		       Image & img,  Image * psat,
		       double saturation)
{
  for (BaseStarIterator it= List->begin(); it!= List->end(); ++it )
    {     
      AddWPsfToImage( psf, (*it)->x, (*it)->y,  (*it)->flux,img,psat,saturation);
    }
}


