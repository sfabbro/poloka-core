#include <fstream>
#include <iomanip>

#include "sestar.h"
#include "image.h"
#include "candstar.h"






//######################## ROUTINES DE READ WRITE ##############



void CandStar::write_scan(ostream &pr) const
{
  pr << star.x << " " 
     << star.y << " " 
     << dass_1 << " " 
     << dass_2 << " " 
     << dass_12 << " "
     << star1.Fluxmax() << " "
     << star2.Fluxmax() << " "
     << star.Fluxmax() << " "
     << star1.flux << " "
     << star2.flux << " "
     << star.flux << " "
     << star.Noise()   << ' ' 
     << star1.Noise() << ' ' 
     << star2.Noise() << ' ' 
     << star1.SigToNoise()<< ' ' 
     << star2.SigToNoise() << ' ' 
     << star.SigToNoise()   << ' ' 
     << star.Faper_ref()   << ' ' 
     << star.PrctIncrease()   << ' ' 
     << star_ref.flux << endl ;
}



void CandStar::WriteHeader_Scan(ostream &pr)
{
  pr << "# x : x position (pixels)" << endl 
     << "# y : y position (pixels)" << endl 
     << "# dass_1 : dist between object on sub and obj. on sub1 " << endl
     << "# dass_2 : dist between object on sub and obj. on sub2 " << endl
     << "# dass_12 : dist between object on sub2 and obj. on sub1 " << endl
     << "# amp1 : Peak pixel value above background on sub1" << endl 
     << "# amp2 : Peak pixel value above background on sub2" << endl 
     << "# amp : Peak pixel value above background on sub" << endl 
     << "# flux1 : flux in pixel units on sub1 " << endl 
     << "# flux2 : flux in pixel units on sub2 " << endl 
     << "# flux : flux in pixel units on sub " << endl ;
  pr << "# noise1 : noise (for sub1) " << endl;  
  pr << "# noise2 : noise (for sub2) " << endl;  
  pr << "# noise : noise (for sub) " << endl;   
  pr << "# stnoise1 : S/N for obj. detected on sub1 " << endl;
  pr << "# stnoise2 : S/N for obj. detected on sub2 " << endl;
  pr << "# stnoise : S/N for obj. detected on sub " << endl;
  pr << "# fap_r : aperture flux computed on ref at candidate position " << endl; 
  pr << "# prcti : percentage increase ie flux/fap_r " << endl; 
  pr << "# hostflux : Flux of the nearest object on the ref " << endl; 
  pr << "# end" << endl;
}


void CandStar::WriteHeader(ostream &pr)
{
  star.WriteHeader_(pr,"");
  star1.WriteHeader_(pr,"1");
  string format = star2.WriteHeader_(pr,"2");
  string formatse = star_ref.WriteHeader_(pr,"r");
  pr << "# dass_r : dist. with nearest object on ref. " << endl;
  pr << "# dass_1 : idem for  obj. detected on sub1 " << endl;
  pr << "# dass_2 : idem for  obj. detected on sub2 " << endl;
  pr << "# dass_12 : distance between detection on sub1 and sub2 " << endl;
  pr << "# end" << endl;
  // if you change the format above, change the format tag here,
  // read it back in CandStar::Read and use the value to "Read" correctly
  pr << "# format " << format << " " << formatse << " CandStar 5 " << endl;
}
  
int CandStar::write(ostream & pr) const
{
  star.writen(pr);
  star1.writen(pr);
  star2.writen(pr);
  star_ref.writen(pr);
  pr << dass_ref << ' ' 
     << dass_1 << ' ' 
     << dass_2 << ' '  
     << dass_12 << ' ' 
     << endl;
  return 1;
}


int CandStar::write_nice(ostream & pr) const
{
  pr << resetiosflags(ios::scientific) ;
  pr << setiosflags(ios::fixed) ;
  pr  << setprecision(3) ;

  pr << endl << "Candidat " << numero << endl ;
  pr << "--------------------------" << endl ;
  pr << "x,y,flux, noise, S/N : " << star.x << " " << star.y << " " 
     << star.flux
     << " " << star.Noise() << " " << star.SigToNoise() << endl ;
  pr << "dass_1, flux1, noise1, S/N1 : " << dass_1 << " " << star1.flux
     << " " << star1.Noise()  << " "  << star1.SigToNoise()  << endl ;
  pr << "dass_2, flux2, noise2 : " << dass_2 << " " << star2.flux
     << " " << star2.Noise()  << " "  << star2.SigToNoise()  << endl ;
  pr << "dass_12 " << dass_12 << endl ;


  double u1 =0., u3=0., u4=0., u5=0., u6=0. ;
  if ( fabs(star_ref.Fluxmax()) > 1.e-10)
    u1 =  star.Fluxmax()/star_ref.Fluxmax() ;
  if ( fabs(star_ref.flux ) > 1.e-10)
    u3 =  star.flux/star_ref.flux ;

  

  pr << endl <<  "dass_r, fluxmax, fluxmaxr, fluxmax/fluxmaxr : " 
     << dass_ref << " | " 
     << star.Fluxmax() << " "  << star_ref.Fluxmax() << " " 
     << u1 << endl ;
  pr <<  "flux, flux_ref ouverture , flux/flux_ref ouverture: " 
     << star.flux << " "  << star.Faper_ref()  << " " 
     << star.PrctIncrease()  << endl ;
  pr <<  "flux, flux_ref  , flux/flux_ref: " 
     << star.flux << " "  << star_ref.flux << " " 
     << u3  << endl ;


  pr << endl <<  "( F(4*seeing) - F(3*seeing) ), (F4-F3/F3) , Noise, (F4-F3)/Noise : "  
     << u4 <<  " "  
     << u5  <<  " " 
     << u6  << endl ;


  pr << endl << "candidat : " << endl;
  star.dumpn(pr);
  pr << endl ;

  pr << "objet associe sur ref: " << endl;
  star_ref.dumpn(pr);
  pr << endl ;
  return 1;
}


void CandStar::read_it(istream& r, const char *Format)
{
  
  star.read_it(r,Format);
  star1.read_it(r,Format);
  star2.read_it(r,Format);
  star_ref.read_it(r,Format);
  r >> dass_ref;
  r >> dass_1;
  r >> dass_2;
  r >> dass_12;
}


int   
CandStarList::write(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  CandStar().WriteHeader(pr);
  for(CandStarCIterator it= begin(); it!= end(); it++)
    {
      (*it).write(pr);
    }
  pr.close();
  return 1 ;
}  

int   
CandStarList::write_nice(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  for(CandStarCIterator it= begin(); it!= end(); it++)
    {
      (*it).write_nice(pr);
    }
  pr.close();
  return 1 ;
}  

int   
CandStarList::write_scan(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  CandStar::WriteHeader_Scan(pr);
  for(CandStarCIterator it= begin(); it!= end(); it++)
    {
      (*it).write_scan(pr);
    }
  pr.close();
  return 1 ;
}  


int CandStarList::read(const string &FileName)
{ 
  ifstream r(FileName.c_str());
  char buff[1024];
  char *format = 0;
  char c;
  while( r >> c ) // pour tester end-of-file
    {
      r.unget() ;
      if ( (c == '#') ) // on saute (presque toujours) la ligne
        {
	  r.getline(buff,1024);
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
	    {
	      format = p + strlen("format");
	    }
        }
      else
        {
	  push_back(CandStar()); /* insert "empty" CandStar at the end */
	  back().read_it(r, format); /* fill it */
	}
    }
  r.close();
  return 1;
}



//##################### END OF CANDSTAR DEFINITIONS #######################################


//########################    ROUTINES DE CUTS    ############################

void
CandStar::ComputeFaperRef(Image & imgref, double rayon, double fd_ref)
{
  star.ComputeFaperRef(imgref,rayon,fd_ref );
  //cout << "Flux sur ref : " <<  star.x << " " << star.y << " " << rayon << " " <<  fd_ref << " " 
  //      << star.Faper_ref() << " " << star.flux << endl ;
  return;
}



void
CandStarList::ComputeFAperRef(Image & imgref, double rayonref, double fd_ref)
{
  
  for (CandStarIterator it= this->begin(); it!= this->end(); ++it )
    it->ComputeFaperRef(imgref,rayonref,fd_ref);
  
}


bool
CandStar::IsGood(DatDetec & dat)
{
  // percentage increase 
  // on compare le flux du cand et le flux ds une ouverture
  // sur la reference
  // si PrctIncrease<0 et Faper_ref <0, 
  // flux > dat.frac_faper_ref*Faper_ref et c'est ok
  if (!star.IsGood(dat.sigtonoise , dat.frac_faper_ref))
    return false;
  
  // test for signal to noise on sub1
  if (!star1.IsGood(dat.sigtonoise_1))
    return false ; 
  
  // test for signal to noise on sub2  
  if (!star2.IsGood(dat.sigtonoise_2))
    return false ;

  // saturated host  
  if (( dass_ref > 0. ) && ( dass_ref < dat.dass_ref_min ) 
      && (star_ref.IsSaturated() ) )
    return false ;  

  return true;  
}


// ROUTINES DE CONSTRUCIONS DE LISTES
void
CandStarList::Construct_List(CandidateStarList & stl,  CandidateStarList & stl1,  
	       CandidateStarList & stl2, SEStarList & stlref, double dist_min)

{
  // ttes les etoiles sont dans le systeme Candidate ie SE
  int n = 0 ;
  for (CandidateStarIterator it= stl.begin(); it!= stl.end(); it++ )
    { 
      CandidateStar *pstar =  *it ;
      // la plus proche sur 1
      CandidateStar *p1 = stl1.FindClosest(*pstar);
      // la plus proche sur 2
      CandidateStar *p2 = stl2.FindClosest(*pstar);
      // la plus proche sur ref
      SEStar *pr = stlref.FindClosest(*pstar);
      if ((p1 != NULL) && (p2 != NULL))
	{
	  double dist1 = p1->Distance(*pstar);
	  double dist2 = p2->Distance(*pstar);
	  double dist = p2->Distance(*p1);
	  if ( ( dist1 < dist_min ) && ( dist2 < dist_min ))
	    {
	      CandStar candstar;
	      candstar.star = *pstar;
	      candstar.star1 = *p1;
	      candstar.star2 = *p2;
	      candstar.dass_1 = dist1;
	      candstar.dass_2 = dist2;
	      candstar.dass_12 = dist;

	      if (pr != NULL)
		{
		  candstar.star_ref  = *pr ;
		  candstar.dass_ref = pr->Distance(*pstar);
		}
	      else
		candstar.dass_ref = -1;
	      (candstar.star).Numero() = n;
	      (candstar.star1).Numero() = n;
	      (candstar.star2).Numero() = n;
	      this->push_back(candstar); n++;
	    }
	}
    }
    
}

// utilise lors de la detection sur 3 sub
void
CandStarList::Construct_Simple(CandidateStarList & stl,  CandidateStarList & stl1,  
			       CandidateStarList & stl2, SEStarList & stlref, Image & imgref,
			       double dass_min, double rayon_ref, double fd_ref)
{
  this->Construct_List(stl, stl1,  stl2, stlref, dass_min);
  this->ComputeFAperRef(imgref, rayon_ref, fd_ref);
}



void
CandStarList::Cut(CandStarList & stlcut, DatDetec & dat)
{  
  int n = 0 ;
  for (CandStarIterator it= this->begin(); it!= this->end(); ++it )
    {

      CandStar &cstar =  *it ;
      if ( cstar.IsGood(dat) )
	{
	  cstar.numero = n;
	  stlcut.push_back(cstar); n++;
	}
    }
}

