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
#include "vutils.h"
#include "gaussianfit.h"


using namespace std;

class LC
{
public :
  vector<double> vec_flux;
  vector<double> vec_eflux;
  int star_number ;
  int npoint ;
  double mag ;
  double mflux ;
  double gmflux ;
  double rmsflux; 
  double grmsflux; 
  LC() : star_number(-1), npoint(0), mag(-1), mflux(-1),  gmflux(-1), rmsflux(-1), grmsflux(-1) {};
  void write(){cerr << "star " << star_number << " : " << mflux << " " << rmsflux << " " 
       << gmflux << " " << grmsflux << " " << endl ;};
  void ComputeMeans() ;
};









void FastComputeMeans(DictFile & l, map<int,LC> & llcc);
void WriteShortHeader_(DictFile & l, ofstream & pr, string suffixe);

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
  bool only_star = false ;
  while ((cc = getopt(argc, argv, "hc:n:o:S")) != -1) 
    {
      switch (cc)
	{
	case 'h' :
	  //usage();
	  break;
	case 'c' :
	  nomc = optarg ;
	  break;
	case 'S' :
	  only_star = true  ;
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

 
  DictFile lc(nomc.c_str()); 
  DictFile ln(nomn.c_str()); 

  // calcul du fluxmoyen
  cerr << "Computing means for " << nomc << " list" << endl ;
  map<int,LC> llcc ;
  FastComputeMeans(lc, llcc) ;
  cerr << "Computing means for " << nomn << " list" << endl ;
  map<int,LC> llcn ;
  FastComputeMeans(ln, llcn) ;

  
  cerr << "Writing Star List " << endl ;
  string nomo_star = CutExtension(nomo) ;
  nomo_star = nomo_star+"_star.list";
  ofstream prls(nomo_star.c_str());
  prls << "#n : " << endl 
       << "#mag : " << endl
       << "#mfc : " << endl
       << "#rmsfc : " << endl
       << "#gmfc : " << endl
       << "#grmsfc : " << endl
       << "#mfn : " << endl
       << "#rmsfn : " << endl
       << "#gmfn : " << endl
       << "#grmsfn : " << endl
       << "#n1 : " << endl
       << "#f1c : " << endl
       << "#f1n : " << endl
       << "#gf1c : " << endl
       << "#gf1n : " << endl
       << "#end  " << endl ;
  double f1c=-1, gf1c=-1, f1n=-1, gf1n=-1 ;
  int n1 = -1 ;
  for(map<int,LC>::iterator it = llcc.begin();  it!=llcc.end(); ++it) 
    {
      int nstarc = it->first ;
      LC & mylcc = it->second ;
      double mag = mylcc.mag ;
      if(mag > 18.5 || mag < 18.) continue ;
      if(llcn.find(nstarc) != llcn.end() )
	{
	  LC & mylcn = llcn[nstarc] ;
	  n1 = nstarc ;
	  f1c = mylcc.mflux ;
	  gf1c = mylcc.gmflux ;
	  f1n = mylcn.mflux ;
	  gf1n = mylcn.gmflux ;
	}
      if(n1 > 0 ) break ;
    }

  for(map<int,LC>::iterator it = llcc.begin();  it!=llcc.end(); ++it) 
    {
      int nstarc = it->first ;
      LC & mylcc = it->second ;
      if(llcn.find(nstarc) == llcn.end() ) continue ;
      LC & mylcn = llcn[nstarc] ;
      prls << setprecision(12);
      prls << nstarc << " " ;
      prls << mylcc.mag << " " ;
      prls << mylcc.mflux << " " << mylcc.rmsflux << " " 
	   << mylcc.gmflux << " " << mylcc.grmsflux << " " 
	   << mylcn.mflux << " " << mylcn.rmsflux << " " 
	   << mylcn.gmflux << " " << mylcn.grmsflux << " " ;
      prls << n1 << " " << f1c << " " << gf1c << " " 
	   << f1n << " " << gf1n << endl ;
    }
  prls.close(); 
 
  if (only_star) return(0) ;


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
	  

  delete [] c_seen ;

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




void LC::ComputeMeans()
{
  double *fflux = new double[vec_flux.size()] ; 
  int ntf=0 ; 
  for(int ii = 0 ; ii < vec_flux.size(); ii++) 
  //Penser a proteger contre les flux nuls ou >0 (done)
    if (vec_flux[ii]>1) {
      fflux[ntf]=vec_flux[ii] ; 
      ntf++;
    }
  npoint = ntf ;
  int nval = ntf ;
  double k = 3.5 ;
  mflux = clipmean(fflux, nval, rmsflux, k);
  int ngval = ntf ;
  k = 3. ;
  gmflux = gaussianfit(fflux, ngval,gmflux, grmsflux, k, /* first_evalutation = */ true) ;
  delete  [] fflux  ;
  write();
}


void FastComputeMeans(DictFile & l,map<int,LC> & llcc )
{
  
 
  l.AddKey("npoints");
  l.AddKey("mnflx");
  l.AddKey("gmnflx");
  l.AddKey("rmsflx");
  l.AddKey("grmsflx");



  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      line->AddKey("npoints", " -10 " );
      line->AddKey("mnflx", " -10 " );
      line->AddKey("gmnflx", " -10 " );
      line->AddKey("rmsflx", " -10 " );
      line->AddKey("grmsflx", " -10 " );

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

      double mag = line->Value("mag");

      double err_flux = line->Value("error");

      llcc[nnstar].mag = mag;
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
      /*int ntot = 0 ;
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
      mylc.meflux = emean ;*/
      mylc.ComputeMeans();
    }
  // on met a jour 
  cerr << "Mise a jour " << endl ;
  for(DictFileIterator line2 = l.begin(); line2 != l.end(); line2++)
    {
      double nstar2 = line2->Value("star");
      int nnstar2 = int(nstar2);
      LC & mylc = llcc[nnstar2];
      //mylc.write() ;
      line2->ModKey("npoints", mylc.npoint);
      line2->ModKey("mnflx", mylc.mflux);
      line2->ModKey("gmnflx", mylc.gmflux);
      line2->ModKey("rmsflx",  mylc.rmsflux);
      line2->ModKey("grmsflx",  mylc.grmsflux);
    }
}











//  old routine
void ComputeMeans(DictFile & l)
{
  
  int *star_numbers = new int[l.size()] ;  
  double *flux = new double[l.size()]; 
  double *err_flux = new double[l.size()];

  l.AddKey("npoints");
  l.AddKey("mnflx");
  l.AddKey("gmnflx");
  l.AddKey("rmsflx");
  l.AddKey("grmsflx");
  int i = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++, i++)
    {
      double nstar = line->Value("star");
      flux[i] = line->Value("flux");
      err_flux[i] = line->Value("error");
      star_numbers[i] = int(nstar);
      line->AddKey("npoints", " -10 " );
      line->AddKey("mnflx", " -10 " );
      line->AddKey("gmnflx", " -10 " );
      line->AddKey("rmsflx", " -10 " );
      line->AddKey("grmsflx", " -10 " );
    }
  int index = 0 ;
  for(DictFileIterator line = l.begin(); line != l.end(); line++)
    {
      cerr << "Computing mean for star " << index << endl ;
      index++;
      double mflux = line->Value("mnflx");
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
	      line2->ModKey("mnflx", mean);
	      line2->ModKey("rmsflx", rms);
	    }
	}
    }

}
