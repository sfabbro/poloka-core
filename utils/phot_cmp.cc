#include <string>
#include <vector>
#include <cmath>

#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <math.h>


#include "dictfile.h"
#include "fileutils.h"


using namespace std;
void FastComputeMeans(DictFile & l);
void ComputeMeans(DictFile & l);
void WriteShortHeader_(DictFile & l, ofstream & pr, string suffixe);
double LocalJulDate(const int day, const int month, const int year, 
	       const int hour, const int min, const double sec);

/*************** ajouter ds dictfile.cc *************/
/*
void DictFileEntry::writen(ofstream & pr) const
{
  for (int ii = 0 ; ii < file.Dict().size() ; ii++)
    pr << elements[ii] << " " ;

}
*/

/*************** ajouter ds dictfile.h *************/
/* dans class DictFileEntry

  void writen(ofstream & pr) const;

*/





int main(int argc, char **argv) {

  char cc ;
  string nomc, nomn, nomo ;
  while ((cc = getopt(argc, argv, "hc:n:o:")) != -1) 
    {
      switch (cc)
	{
	case 'h' :
	  //usage();
	  break;
	case 'c' :
	  nomc = optarg ;
	  break;
	case 'n' :
	  nomn = optarg ;
	  break;
	case 'o' :
	  nomo = optarg ;
	  break;
	default:
	  //usage();
	  cerr << "bad option " << endl ;
	}
    }
  double mjdref = 52640.0 ;

 
  DictFile lc(nomc.c_str()); 
  DictFile ln(nomn.c_str()); 

  // calcul du fluxmoyen
  cerr << "Computing means for " << nomc << " list" << endl ;
  FastComputeMeans(lc) ;
  cerr << "Computing means for " << nomn << " list" << endl ;
  FastComputeMeans(ln) ;

  cerr << "Joining lists " << endl ;
  ofstream pr(nomo.c_str());
  WriteShortHeader_(lc,pr,  "c");
  WriteShortHeader_(ln,pr,  "n");
  pr << "#end" << endl ;


  int *c_seen = new int[lc.size()] ;
  for(int i = 0 ; i < lc.size() ; i++) c_seen[i] = -1 ;
  int count = 0 ;
  for(DictFileIterator it =  ln.begin(); it !=ln.end(); it++, count++)
    {
      double nstar = it->Value("star");
      int nnstar = int(nstar);
      string name="";
      if (it->HasKey("name")) name = (string)  it->Value("name");
      bool is_seen = false ;
      int nc = 0 ;
      for(DictFileIterator itc =  lc.begin(); itc !=lc.end(); itc++, nc++)
	{	  
	  double nstarc = itc->Value("star");
	  int nnstarc = int(nstarc);

	  string namec="" ;
	  if (itc->HasKey("name")) 
	    {
	      namec = (string) itc->Value("name");
	    }
	  bool is_same = false ;
	  if (( name != "" ) && ( namec != "" ))
	    {
	      if ( name == namec ) is_same = true ;
	    }
	  else
	    {
	      if ( !it->HasKey("mjd")) cerr << "PAS DE CLEF MJD" << endl ;
	      if ( !itc->HasKey("mjd")) cerr << "PAS DE CLEF MJD C" << endl ;
	      double mjd = it->Value("mjd");
	      double mjdc = itc->Value("mjd");
	      if (count < 1) cerr << setprecision(10) << mjd << " " << mjdc << endl ;
	      if (fabs(mjd-mjdc)<0.0000001)
		 is_same = true ;

	    }
	  if (( nnstarc == nnstar ) && is_same)
	    {
	      itc->writen(pr);
	      it->writen(pr);
	      pr << endl ;
	      if(c_seen[nc] > 0 ) cerr << " ERROR : starc " << nnstarc << " + date " 
				       << namec << "seen twice !!! " << endl ;
	      c_seen[nc] = 1 ;
	      is_seen = true ;
	      break ;
	    } 	  
	}
      if (! is_seen )
	{
	  cerr << " starn  " << nnstar << " + date " 
	       << name << " not seen " << endl ;
	   for(int ie = 0 ; ie < lc.Dict().size() ; ie++)
	    pr << " -1 " ;
	   it->writen(pr);
	   pr << endl ;
	}
    }
  int nc = 0 ;
  for(DictFileIterator itc =  lc.begin(); itc !=lc.end(); itc++, nc++)
    {
     if(c_seen[nc] < 0 ) 
       {
	 double nstarc = itc->Value("star");
	  int nnstarc = int(nstarc);
	  string namec  = "" ;
	  if ( itc->HasKey("name")) namec = (string) itc->Value("name");
	  cerr << " starc " << nnstarc << " + date " 
	       << namec << " not seen " << endl ;
	  itc->writen(pr);
	  for(int ie = 0 ; ie < ln.Dict().size() ; ie++)
	   pr << " -1 " ;
	  pr << endl ;
	}
    }

  pr.close();
	  



}

void WriteShortHeader_(DictFile & l, ofstream & pr, string suffixe) {

  unsigned presentSize = l.Dict().size();


  // do not write global keys and values
 

 // invert the dictionnary:
  map<int,string> tags;
  for (Dictionnary::const_iterator it = l.Dict().begin(); it != l.Dict().end(); ++it)
    tags[it->second] = it->first;

  //write the header
  for (unsigned i = 0; i < presentSize; ++i)
    pr << "# " <<  tags[i] << suffixe << " : " << endl ;
}


void ComputeMeans(DictFile & l)
{
  
  int *star_numbers = new int[l.size()] ;  
  double *flux = new double[l.size()]; 
  double *err_flux = new double[l.size()];

  l.AddKey("npoints");
  l.AddKey("mnflux");
  l.AddKey("mnerr");
  l.AddKey("rmsflux");
  int i = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++, i++)
    {
      double nstar = line->Value("star");
      flux[i] = line->Value("flux");
      err_flux[i] = line->Value("error");
      star_numbers[i] = int(nstar);
      line->AddKey("npoints", " -10 " );
      line->AddKey("mnflux", " -10 " );
      line->AddKey("mnerr", " -10 " );
      line->AddKey("rmsflux", " -10 " );
    }
  int index = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      cerr << "Computing mean for star " << index << endl ;
      index++;
      double mflux = line->Value("mnflux");
      if (mflux > -10) continue ;
      double nstar = line->Value("star");
      int nnstar = int(nstar);
      int cc=0 ;
      int ntot = 0 ;
      double S = 0, S2=0 ;
      double Se = 0 ;
      double mean=-1, emean=-1, rms=-1 ;
      for(int ic = 0 ; ic < l.size() ; ic++)
	{
	  if (star_numbers[ic] != nnstar) continue ;
	  double f = flux[ic] ;
	  double ef = err_flux[ic] ;
	  S += f ;
	  S2 += f*f ;
	  Se += ef ;
	  ntot++;
	}
      if ( ntot > 1 )
	{
	  mean = S/ntot ;
	  emean = Se/ntot ;
	  rms = (S2 - ntot*mean*mean)/(1.*(ntot-1));
	  if (rms > 0 )
	    rms = sqrt(rms);
	}
      // on update non seulement cette etoile mais ttes celles qui ont le meme numero
      for(DictFileIterator line2 = l.begin(); line2 != l.end(); line2++)
	{
	  double nstar2 = line2->Value("star");
	  int nnstar2 = int(nstar2);
	  if ( nnstar2 == nnstar )
	    {
	      line2->ModKey("npoints", ntot);
	      line2->ModKey("mnflux", mean);
	      line2->ModKey("rmsflux", rms);
	      line2->ModKey("mnerr", emean);
	    }
	}
    }

}


class LC
{
public :
  vector<double> vec_flux;
  vector<double>vec_eflux;
  int star_number ;
  int npoint ;
  double mflux ;
  double meflux ;
  double rmsflux; 
  LC() : star_number(-1), npoint(0), mflux(-1), meflux(-1), rmsflux(-1) {};
};

void FastComputeMeans(DictFile & l)
{
  
 
  l.AddKey("npoints");
  l.AddKey("mnflux");
  l.AddKey("mnerr");
  l.AddKey("rmsflux");


  map<int,LC> llcc ;

  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      line->AddKey("npoints", " -10 " );
      line->AddKey("mnflux", " -10 " );
      line->AddKey("mnerr", " -10 " );
      line->AddKey("rmsflux", " -10 " );

      if (line->HasKey("name") )
	{
	  string thename = line->Value("name");
	  string pattern = "p" ;
	  RemovePattern(thename, pattern);
	  line->ModKey("name", thename);
	}
      
      double nstar = line->Value("star");
      int nnstar = int(nstar) ;

      double flux = line->Value("flux");
      double err_flux = line->Value("error");

      llcc[nnstar].star_number = nnstar;
      llcc[nnstar].npoint += 1 ;
      llcc[nnstar].vec_flux.push_back(flux);
      llcc[nnstar].vec_eflux.push_back(err_flux);
    }

  // on clacul les moyennes
  for(map<int,LC>::iterator it = llcc.begin();  it!=llcc.end(); ++it) 
    {
      int nstar = it->first ;
      cerr << "Computing mean for star " << nstar << endl ;
      LC & mylc = it->second ;
      int ntot = 0 ;
      double S = 0, S2=0 ;
      double Se = 0 ;
      double mean=-1, emean=-1, rms=-1 ;
      for(int ic = 0 ; ic < mylc.vec_flux.size() ; ic++)
	{
	  double f = mylc.vec_flux[ic] ;
	  double ef = mylc.vec_eflux[ic] ;
	  S += f ;
	  S2 += f*f ;
	  Se += ef ;
	  ntot++;
	}
      if ( ntot > 1 )
	{
	  mean = S/ntot ;
	  emean = Se/ntot ;
	  rms = (S2 - ntot*mean*mean)/(1.*(ntot-1));
	  if (rms > 0 )
	    rms = sqrt(rms);
	}
      mylc.npoint = ntot ;
      mylc.mflux = mean ;
      mylc.rmsflux = rms ;
      mylc.meflux = emean ;
    }
  // on met a jour 
  for(DictFileIterator line2 = l.begin(); line2 != l.end(); line2++)
    {
      double nstar2 = line2->Value("star");
      int nnstar2 = int(nstar2);
      LC & mylc = llcc[nnstar2];
      line2->ModKey("npoints", mylc.npoint);
      line2->ModKey("mnflux", mylc.mflux);
      line2->ModKey("rmsflux",  mylc.rmsflux);
      line2->ModKey("mnerr", mylc.meflux);
    }
}



