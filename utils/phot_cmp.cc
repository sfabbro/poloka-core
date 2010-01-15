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
  //double mjdref = LocalJulDate(01, 01, 2003, 0, 0, 0)-2400000.5;
  //cerr << setprecision(12) << setiosflags(ios::fixed) << mjdref << endl ;

  /*DictFile l1(nomc.c_str());
  ComputeMeans(l1) ;
  l1.Write("essai.list");
  DictFile l2(nomc.c_str());
  FastComputeMeans(l2) ;
  l2.Write("essai_fast.list");

  exit(0);*/

 
  DictFile lc(nomc.c_str()); 
  DictFile ln(nomn.c_str()); 

  // calcul du fluxmoyen
  cerr << "Computing means for " << nomc << endl ;
  FastComputeMeans(lc) ;
  cerr << "Computing means for " << nomn << endl ;
  FastComputeMeans(ln) ;

  cerr << "Joining lists " << endl ;
  ofstream pr(nomo.c_str());
  WriteShortHeader_(lc,pr,  "c");
  WriteShortHeader_(ln,pr,  "n");
  pr << "#end" << endl ;


  int *c_seen = new int[lc.size()] ;
  for(int i = 0 ; i < lc.size() ; i++) c_seen[i] = -1 ;
  int count = 0 ;
  for(DictFileIterator it =  ln.begin(); it !=ln.end(); it++)
    {
      double nstar = it->Value("star");
      int nnstar = int(nstar);
      // double mmjd  = it->Value("mmjd");
      string name  = it->Value("name");
      bool is_seen = false ;
      int nc = 0 ;
      for(DictFileIterator itc =  lc.begin(); itc !=lc.end(); itc++, nc++)
	{	  
	  double nstarc = itc->Value("star");
	  int nnstarc = int(nstarc);
	  //double mjdc  = itc->Value("mjd");
	  //double mmjdc  = mjdc - mjdref ;
	  string namec  = itc->Value("name");
	  //if (( nnstarc == nnstar ) && ( fabs(mmjd-mmjdc)<1.e-5))
	  if (( nnstarc == nnstar ) && ( name == namec ))
	    {
	      itc->writen(pr);
	      it->writen(pr);
	      pr << endl ;
	      if(c_seen[nc] > 0 ) cerr << " ERROR : starc " << nnstarc << " + date " 
				       << namec << "seen twice !!! " << endl ;
	      //			       << mmjd << "seen twice !!! " << endl ;
	      c_seen[nc] = 1 ;
	      is_seen = true ;
	      break ;
	    } 	  
	}
      if (! is_seen )
	{
	  cerr << " starn  " << nnstar << " + date " 
	       << name << " not seen " << endl ;
	  //	       << mmjd << " not seen " << endl ;
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
	  //double mjdc  = itc->Value("mjd");
	  //double mmjdc = mjdc-mjdref ;
	  string namec  = itc->Value("name");
	  cerr << " starc " << nnstarc << " + date " 
	       << namec << " not seen " << endl ;
	  //	       << mmjdc << " not seen " << endl ;
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
  l.AddKey("mean_flux");
  l.AddKey("mean_err");
  l.AddKey("rms_flux");
  int i = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++, i++)
    {
      double nstar = line->Value("star");
      flux[i] = line->Value("flux");
      err_flux[i] = line->Value("error");
      star_numbers[i] = int(nstar);
      line->AddKey("npoints", " -10 " );
      line->AddKey("mean_flux", " -10 " );
      line->AddKey("mean_err", " -10 " );
      line->AddKey("rms_flux", " -10 " );
    }
  int index = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      cerr << "Computing mean for star " << index << endl ;
      index++;
      double mflux = line->Value("mean_flux");
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
	      line2->ModKey("mean_flux", mean);
	      line2->ModKey("rms_flux", rms);
	      line2->ModKey("mean_err", emean);
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
  l.AddKey("mean_flux");
  l.AddKey("mean_err");
  l.AddKey("rms_flux");


  map<int,LC> llcc ;

  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      line->AddKey("npoints", " -10 " );
      line->AddKey("mean_flux", " -10 " );
      line->AddKey("mean_err", " -10 " );
      line->AddKey("rms_flux", " -10 " );
      
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
      line2->ModKey("mean_flux", mylc.mflux);
      line2->ModKey("rms_flux",  mylc.rmsflux);
      line2->ModKey("mean_err", mylc.meflux);
    }
}


double LocalJulDate(const int day, const int month, const int year, 
	       const int hour, const int min, const double sec)
{

  /* decimal day fraction	*/
  double frac = (( double)hour/ 24.0)
    + ((double) min / 1440.0)
    + (sec / 86400.0);

  /* convert date to format YYYY.MMDDdd	*/
  double gyr = (double) year
    + (0.01 * (double) month)
    + (0.0001 * (double) day)
    + (0.0001 * frac) + 1.0e-9;

  /* conversion factors */
  long iy0, im0;
  if ( month <= 2 )
    {
      iy0 = year - 1L;
      im0 = month + 12;
    }
  else
    {
      iy0 = year;
      im0 = month;
    }
  long ia = iy0 / 100L;
  long ib = 2L - ia + (ia >> 2);

  /* calculate julian date	*/
  long ljd;
  if ( year <= 0L )
    ljd = (long) ((365.25 * (double) iy0) - 0.75)
      + (long) (30.6001 * (im0 + 1L) )
      + (long) day + 1720994L;
  else
    ljd = (long) (365.25 * (double) iy0)
      + (long) (30.6001 * (double) (im0 + 1L))
      + (long) day + 1720994L;

  /* on or after 15 October 1582	*/
  if ( gyr >= 1582.1015 )ljd += ib;

  return (double) ljd + frac + 0.5;
}	

