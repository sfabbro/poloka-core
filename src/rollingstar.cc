#include <fstream>

#include "rollingstar.h"
#include "nstarmatch.h"
#include "candidatestar.h"
#include "candidatestar.h"
#include "fastfinder.h"
#include "dodetection.h"
#include "reducedimage.h"
#include "fitsimage.h"
//#include "datdetect.h"
#include "toadscards.h"

#ifdef STORAGE
RollingStar::RollingStar()
  : BaseStar(0.,0.,0.)
{
  //  Set_to_Zero();
}
#endif

RollingStar::RollingStar(int & n)
{
  Nlists = n;
  flux = new double[Nlists];
  eflux = new double[Nlists];
  stnoise = new double[Nlists];
  prcti = new double[Nlists];
  chi2 = new double[Nlists];
  
};

RollingStarList::RollingStarList(const NStarMatchList & nstml, const ReducedImage &imRef)
{

  Nlists = nstml.Nlists();
  DatDetec datdet(DefaultDatacards()); 
  double seeing_ref = imRef.Seeing();
  double rayon = datdet.n_rayon * seeing_ref; 
  
  FitsImage imgref(imRef.FitsName());
  double fd_ref = 0 ;
  SEStarList stlref(imRef.CatalogName()); 
  FastFinder finder(*((BaseStarList*) &stlref));

  for (NStarMatchCIterator si = nstml.begin(); si != nstml.end(); ++si)
    {
      double xbar = 0;
      double ybar = 0;
      double sumx = 0;
      double sumy = 0;
      RollingStar *rst = new RollingStar(Nlists);
      
      
      for (int i =0; i< Nlists; i++)
	{
	  
	  if(si->StarExist(i))
	    {
	      const CandidateStar *current = (CandidateStar*) si->star[i];
	      double vx = current->Sigx();
	      if ( current->x < 0 || vx < 0)
		{
		  cout << current->x << " " << vx << endl;	      
		  continue;
		}
	      double vy = current->Sigy();
	      if ( current->y < 0 || vy < 0)
		{
		  cout << current->y << " " << vy << endl;
		  continue;
		}
	      vx *= vx;
	      vy *= vy;
	      xbar += current->x / vx;
	      ybar += current->y / vy;
	      sumx += 1/vx;
	      sumy += 1/vy;

	      rst->flux[i] = current->flux;
	      rst->eflux[i] = current->EFlux();
	      rst->stnoise[i] = current->SigToNoise();
	      rst->chi2[i] = current->Chi2();
	      
	    }
	  else 
	    {
	      rst->flux[i] = 0;
	      rst->eflux[i] = 0;
	      rst->stnoise[i] = 0;
	      rst->chi2[i] = 100000;
	    }
	  rst->x = xbar/sumx;
	  rst->y = ybar/sumy;
	  
	  //cout << rst->x << " " << rst->y << endl;
          int npix;
	  double faper = Flux_Aperture(imgref, rst->x, rst->y, rayon, fd_ref,npix  );
	  
	  for (int i =0; i< Nlists; i++)
	    {
	      if(si->StarExist(i))
		rst->prcti[i] = rst->flux[i]/faper;
	      else
		rst->prcti[i]=0;
	    }
	  BaseStar *pstar =  rst ;
	  
	  // la plus proche sur ref
	  const BaseStar *pr = finder.FindClosest(*pstar,100);
	  
      
	  if (pr != NULL)
	    {
	      rst->host = *pr ;
	      
	      rst->dass_ref = pr->Distance(*pstar);
	    }
	  else
	    rst->dass_ref = -1;
	}
      this->push_back(rst);
      
    } 
}
  
void RollingStar::WriteHeader(ostream & pr)
{
  pr << "# x : x position of the canidate (pixels)" << endl 
     << "# y : y position of the canidate (pixels)" << endl ;
  
  BaseStar::WriteHeader_(pr,"ref");
  pr << "# dass_ref : dist. with nearest object on ref." << endl ;
  
  for (int i=1; i<Nlists+1;i++)
  {
    pr << "# flux"<< i <<" : flux in pixel units on sub"<<i << endl 
       << "# eflux"<< i <<" : error on the flux" <<i<< endl;
    pr << "# stnoise"<< i <<" : S/N for obj. detecetd on sub" << i << endl;
    pr << "# prcti"<< i <<" : percentage increase ie flux" << i << "/fap_r " << endl; 
    pr << "# chi2"<< i <<" : chi2 for the candidate " << i  << endl; 
    

  }
#ifdef STORAGE  
  for (int i =0; i< Nlists ; i++)
    for (int j = i+1; j < NLists ; j++)
      pr << "# dass"<< i << j <<" : association distance for "
	 <<i<<" and "<<j<< endl ; 
#endif  
  pr << "# end " << endl ;  static char format[256];
  sprintf(format,"RollingStar %d",1 );
  //return format; 
}


void RollingStar::write(ostream &pr) const
{
  
       
  pr << x << " "; 
  pr << y << " "; 
  host.writen(pr);
  pr << " " << dass_ref << " " ;
  
  for (int i =0; i< Nlists ; i++)
    {
      pr << flux[i] << " ";
      pr << eflux[i] << " ";
      pr << stnoise[i]   << " ";
      pr << prcti[i]   << " "; 
      pr << chi2[i]   << " "; 
    }
  pr << endl ;
#ifdef STORAGE 
  for (int i =0; i<this->nlists ; i++)
     
    
    for (int j = i+1; j<this->nlists ; j++)
      pr << Nstarm.Distance(i,j, *Transfo(i,j))<< " " ;
#endif
  
}

void RollingStarList::write(const string &FileName) 
{
  ofstream pr(FileName.c_str());
  
  
  RollingStar toto(Nlists);
  toto.WriteHeader(pr);

  for(RollingStarCIterator it= this->begin(); it!= this->end(); it++)
    {
      (*it)->write(pr);
    }
  pr.close();
}  
