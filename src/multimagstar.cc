#include <iostream>
#include <cmath> // asin, sqrt
#include <iomanip>

#include "multimagstar.h"

#include "image.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include "globalval.h"

#include "fastifstream.h"

#include "listmatch.h"


void
MultiMagSEStar::SetToZero()
{

  alpha = 0. ;
  delta = 0. ;
  x_orig = 0. ;
  y_orig = 0. ;
  gx=0. ;
  gy=0. ;
  gmxx = 0. ;
  gmyy = 0. ;
  gmxy = 0. ;
  gmxx_loc = 0. ;
  gmyy_loc = 0. ;
  gmxy_loc = 0. ;
  ell_dist = 0. ;
  norm_dist = 0. ;
  dist = 0. ;

  star_dist = -1 ;
  //star = NULL ;
}



// CONSTRUCTORS
MultiMagSEStar::MultiMagSEStar(const SEStar &sestar, double *phot_autoaper) : SEStar(sestar)
{
  SetToZero();
  x_orig = x ;
  y_orig = y;
  ell_aper.SetParameters(sestar,1.,0.5*phot_autoaper[0],0.5*phot_autoaper[1]);
  // les magbox sont crees et initialisees a la demande
}






std::string MultiMagSEStar::WriteHeader_(ostream & pr, 
				     const char* i) const
{
  std::string sestarFormat = SEStar::WriteHeader_(pr, i);
  if (!i) i="";
  pr << "#ra" << i  << " : " << endl;
  pr << "#dec" << i  << " : " << endl;
  pr << "#x_orig" << i  << " : " << endl;
  pr << "#y_orig" << i  << " : " << endl;

  pr << "#gx" << i  << " : " << endl;
  pr << "#gy" << i  << " : " << endl;
  pr << "#gmxx" << i  << " : " << endl;
  pr << "#gmyy" << i  << " : " << endl;
  pr << "#gmxy" << i  << " : " << endl;
  pr << "#gmxxl" << i  << " : " << endl;
  pr << "#gmyyl" << i  << " : " << endl;
  pr << "#gmxyl" << i  << " : " << endl;

  char cc1[500];
  sprintf(cc1,"es%s",i);
  char cc2[500];
  sprintf(cc2,"eg%s",i);
  ell_aper.WriteHeader_(pr,cc1);
  g_ell_aper.WriteHeader_(pr,cc2);



  pr << "#nm"<<i<< " : number of magboxes " << endl;



  for (unsigned k=0; k < magboxes.size(); ++k)
    {
      string kk = magboxes[k].calib.band ;
      //pr << "#band" << kk << i  << " : " << endl;
      pr << "#zp" << kk << i  << " : " << endl;
      pr << "#sigzp" << kk << i  << " : " << endl;
      pr << "#seeing" << kk << i  << " : " << endl;
      pr << "#fc" << kk << i  << " : " << endl;
      pr << "#efc" << kk << i  << " : " << endl;
      pr << "#f" << kk << i  << " : " << endl;
      pr << "#ef" << kk << i  << " : " << endl;
      pr << "#m_a" << kk << i  << " : " << endl;
      pr << "#em_a" << kk << i  << " : " << endl;
      pr << "#m_c" << kk << i  << " : " << endl;
      pr << "#em_c" << kk << i  << " : " << endl;
      pr << "#m" << kk << i  << " : " << endl;
      pr << "#em" << kk << i  << " : " << endl;
      pr << "#fmx" << kk << i  << " : fluxmax " << endl;
      pr << "#flag" << kk << i  << " : flag " << endl;
      pr << "#flagb" << kk << i  << " : flagbad " << endl;

    }

  return sestarFormat+ "MultiMagSEStar 7"; 
}



void
MultiMagSEStar::dumpn(ostream& s) const
{
  SEStar::dumpn(s);
  s << "alpha : " << alpha << " " 
    << "delta : " << delta << " " 
    << "x_orig : " << x_orig << " " 
    << "y_orig : " << y_orig << " " 
    << "gx : " << gx << " "    
    << "gy : " << gy << " "    
    << "gmxx : " << gmxx << " "   
    << "gmyy : " << gmyy << " "   
    << "gmxy : " << gmxy << " "  
    << "gmxx_loc : " << gmxx_loc << " "  
    << "gmyy_loc : " << gmyy_loc << " "   
    << "gmxy_loc : " << gmxy_loc << " "   ;

   s << "size magboxes : " << magboxes.size() << ' ';
  for (unsigned k=0; k < magboxes.size(); ++k)
    {
      const CalibBox &cal =  magboxes[k].calib;
      s << "band : " << cal.band << " ZP : " << cal.ZP 
	<< " sigZP :  " << cal.sigZP << " seeing : " << cal.seeing << " " ;
      s <<  " flux circ : " << magboxes[k].f_circ 
	<< " err. flux circ : " <<  magboxes[k].ef_circ 
       <<  " flux auto : " << magboxes[k].f_auto 
	<< " err. flux auto : " <<  magboxes[k].ef_auto 
	<< " mag auto : " <<  magboxes[k].m_auto 
	<< " err. mag auto : " <<  magboxes[k].em_auto << " " 
	<< " mag circ : " <<  magboxes[k].m_circ 
	<< " err. mag circ : " <<  magboxes[k].em_circ << " " ; 
      s << " fluxmax : " << magboxes[k].fluxmax << " " ;
      s << " flag : " << magboxes[k].flag << " " ;
      s << " flagbad : " << magboxes[k].flagbad << " " ;
    }
}

void MultiMagSEStar::writen(ostream& s) const
{
  SEStar::writen(s); 
  ios::fmtflags  old_flags =  s.flags();
  int oldprec = s.precision(); 
  s<< setiosflags(ios::fixed) ;
  s  << setiosflags(ios::scientific) ;
  s << setprecision(14);
  s << alpha << " " ;
  s << delta << " " ;
  s << x_orig << " " ;
  s << y_orig << " " ; 
  s << gx << " " ;   
  s << gy << " " ;   
  s << gmxx << " " ;   
  s << gmyy << " " ;   
  s << gmxy << " " ;  
  s << gmxx_loc << " " ;   
  s << gmyy_loc << " " ;   
  s << gmxy_loc << " " ;  

  ell_aper.writen(s);
  g_ell_aper.writen(s);

    
 s << magboxes.size() << ' ';
  for (unsigned k=0; k < magboxes.size(); ++k)
    {
     
       const CalibBox &cal =  magboxes[k].calib;
       //s << cal.band << " " ; //!! attention, pas lisible par paw du coup
       s << cal.ZP << " " ;
       s << cal.sigZP << " " ;
       s << cal.seeing << " " ;
       s << magboxes[k].f_circ << " " ;
       s << magboxes[k].ef_circ << " " ;
       s << magboxes[k].f_auto << " " ;
       s << magboxes[k].ef_auto << " " ;
       s << magboxes[k].m_auto << " " ;
       s << magboxes[k].em_auto << " " ; 
       s << magboxes[k].m_circ << " " ;
       s << magboxes[k].em_circ << " " ;  
       s << magboxes[k].m << " " ;
       s << magboxes[k].em << " " ; 
       s << magboxes[k].fluxmax << " " ;
       s << magboxes[k].flag << " " ;
       s << magboxes[k].flagbad << " " ;

    }
 s.flags(old_flags);
 s << setprecision(oldprec);
    
}



void MultiMagSEStar::read_it(fastifstream& r, const char* Format)
{
  SEStar::read_it(r, Format);
  int format = DecodeFormat(Format, "MultiMagSEStar");
  r >> alpha  ;
  r >> delta  ;
  r >> x_orig  ;
  r >> y_orig  ;
  
  if (format >= 4)
    {
      r >> gx ;
      r >> gy ;
      r >> gmxx ;
      r >> gmyy ;
      r >> gmxy ;
      r >> gmxx_loc ;
      r >> gmyy_loc ;
      r >> gmxy_loc ;
    } 
  if (format >= 5)
    {
      ell_aper.read(r);
      g_ell_aper.read(r);
    }
  if (format <= 1)
    {
      double undouble ;
      for(int i = 0 ; i < 35 ; i++)
	r >> undouble ;
    }

  unsigned nmag;
  r >> nmag; 
  for (unsigned k=0; k < nmag; ++k)
    {
      magboxes.push_back(ShortMagBox());
      ShortMagBox &mb = magboxes.back();
      CalibBox &cal = mb.calib;
      if (format <= 2)
	r >> cal.band ;
      r >> cal.ZP  ;
      r >> cal.sigZP  ;
      r >> cal.seeing  ;
      r >> mb.f_circ ;
      r >> mb.ef_circ ;
      r >> mb.f_auto ;
      r >> mb.ef_auto ;
      if (format >= 7)
	{  	  
	  r >> mb.m_auto ;
	  r >> mb.em_auto ; 
	  r >> mb.m_circ ;
	  r >> mb.em_circ ; 
	  r >> mb.m ;
	  r >> mb.em ;
	}
      else
	{
	  r >> mb.m_auto ;
	  r >> mb.em_auto ; 
	}
      if (format >= 6)
	{
	  r >> mb.fluxmax ;
	  r >> mb.flag ;
	  r >> mb.flagbad ;
	}
   }
}

MultiMagSEStar* MultiMagSEStar::read(fastifstream& r, const char* Format)
{
  MultiMagSEStar *pstar = new MultiMagSEStar();  
  pstar->read_it(r, Format);
  return(pstar);
}



void  MultiMagSEStar::ComputeMag(int kbox, string band, double ZP, double eZP, double & old_ZP)

{
  ShortMagBox &mb = magboxes[kbox];
  CalibBox &cal = mb.calib;  
  old_ZP = cal.ZP ;
  //double zz =  cal.ZP - ZP ;
  //cerr << "band : " << zz  << endl ;

  cal.band = band ;
  cal.ZP = ZP;
  cal.sigZP = eZP;
  if ( (mb.f_auto > 0) && (mb.ef_auto > 0) )
    {			 			  
      mb.m_auto = -2.5*log10(mb.f_auto)+ZP ;
      mb.em_auto = 2.5/log(10.)*mb.ef_auto/mb.f_auto; 
    }
  else
    {
      mb.m_auto = 99 ;
      mb.em_auto = 99 ;
    }
  if ((mb.f_circ > 0) && (mb.ef_circ > 0) )
    {
      mb.m_circ = -2.5*log10(mb.f_circ)+ZP ;
      mb.em_circ = 2.5/log(10.)*mb.ef_circ/mb.f_circ;
    }
  else
    {
      mb.m_circ = 99 ;
      mb.em_circ = 99 ;
    }
  mb.m = mb.m_auto ;
  mb.em = mb.em_auto ;
}



#include "ellipticaper.h"

// equation de l'ellipse : cxx . X^2 + cyy . Y^2 + cxy XY = 1, ellipse de 1/2 grd axe A, de 1/2 petit axe B,  et d'inclinaison theta
// d^2 = cxx x^2 + cyy y^2 + cxy y^2 = distance elliptique en nombre de A() ou B(), ca vaut 1 sur le pourtour de l'ellipse.
// le flux dans SExtractor est calcule sur l'ellipse:
// cxx . X^2 + cyy . Y^2 + cxy XY = R^2
// R = rayon d'integration = kron FACTOR  = 2.5 * kron RADIUS
// qui est ce qu'on recupere de SExtractor et malhencontreusement appele >Kronradius() chez nous.
// donc pour retranscrire d en "rayon de l'ellipse de photom en pixels", il faut x par R*sqrt(AB)
// si R * sqrt(AB) < RadMin= PHOT_AUTOAPER[0]*0.5 ds la datacard, alors on suit la procedure de SExtractor : reset les cxx,cyy,cxy a 1,1,0 (cercle) et R =  Radius = = PHOT_AUTOAPER[1]*0.5, en pixel donc.
double MultiMagSEStar::SqEllipticDistance(double xx, double yy) const
{
  double dist2 =  ell_aper.SqEllipticDistance(xx,yy);
  return dist2;
}
// est -on ou non dans l'ellipse (ou le cercle) de photometrie ?
// d/Radius si cercle, d/kron factor ie KronRadius() si ellipse.

double MultiMagSEStar::NormalizedDistance(double xx, double yy) const
{
  double dist =  ell_aper.NormalizedEllipticDistance(xx,yy);
  return(dist);
}




bool IncSqEllipticDist(const MultiMagSEStar *S1, const MultiMagSEStar *S2)
{
  return (S1->ell_dist < S2->ell_dist );
}

bool IncNormalizedDistance(const MultiMagSEStar *S1, const MultiMagSEStar *S2)
{
  return (S1->norm_dist < S2->norm_dist );
}




#include "starlist.cc"

//instantiate all template routines from StarList
template class StarList<MultiMagSEStar>;
MultiMagSEStarList::MultiMagSEStarList(const string &FileName)
{
  read(FileName);
  int nbox = front()->magboxes.size() ;
  // 2 cas :
  // 1) old list : pas de globals, mais etiquette cal.band lue

  // 2) new list : globals contenant les etiquettes a mettre dans cal.band
  // verification des globals
   if (nbox >  0 )
    {
      bool new_version = ( (front()->magboxes[0].calib).band == "" );
      // verification homogeneite
      for (int ii = 0 ; ii < nbox ; ii++)
	{
	  CalibBox &cal =  front()->magboxes[ii].calib ;
	  if ( (( cal.band == "") && ! new_version ) || ( ( cal.band != "") && new_version ) )
	    cerr << "Error in reading filename : " << FileName << "  cal bands partially empty" << endl ;
	}
      if ( ! new_version ) // on mets les etiquettes de chaque magbox en globals
	{
	  for (int ii = 0 ; ii < nbox ; ii++)
	    {
	      CalibBox &cal =  front()->magboxes[ii].calib ;
	      char  kkey[20];
	      sprintf(kkey,"MAG_%d",ii);
	      GlobVal().AddKey(kkey,cal.band) ;
	    }
	}
      else // on mets les etiquettes des globals ds les magbox
	{
	  // verification autant de magboxes que d'etiquettes
	  for(int ii = 0 ; ii < nbox ; ii++)
	    {
	      char  kkey[20];
	      sprintf(kkey,"MAG_%d",ii);
	      if ( ! GlobVal().HasKey(kkey) )
		cerr << "Error in reading filename : " << FileName 
		     << "  no " << kkey << endl;
	    }
	  // on remets les etiquettes
	  for (MultiMagSEStarIterator mit = this->begin();mit != this->end();mit++)
	    {
	      for (unsigned int ii = 0 ; ii < (*mit)->magboxes.size(); ii++ )
		{
		  CalibBox &cal =  (*mit)->magboxes[ii].calib ;
		  char  kkey[20];
		  sprintf(kkey,"MAG_%d",ii);
		  cal.band = GlobVal().getStringValue(kkey) ;
		}
	    }
	}
    }
}


// just copy and initialize
MultiMagSEStarList::MultiMagSEStarList(const SEStarList &L)
{ 
  double phot_autoaper[2] = {16., 16.}  ;
  if ( L.GlobVal().HasKey("PHOT_AUTOAPER_0") )
    phot_autoaper[0] = L.GlobVal().getDoubleValue("PHOT_AUTOAPER_0") ;
  if ( L.GlobVal().HasKey("PHOT_AUTOAPER_1") )
    phot_autoaper[1] = L.GlobVal().getDoubleValue("PHOT_AUTOAPER_1") ;
  cerr << "PHOT_AUTOAPERS : " << phot_autoaper[0] << " " << phot_autoaper[1]<< endl ;
  GlobVal().AddKey("PHOT_AUTOAPER_0",phot_autoaper[0]) ;
  GlobVal().AddKey("PHOT_AUTOAPER_1",phot_autoaper[1]) ;


  for (SEStarCIterator i = L.begin(); i != L.end(); ++i)
    push_back(new MultiMagSEStar(**i, phot_autoaper));
}



void
MultiMagSEStarList::CopySEStarList(const SEStarList &L)
{ 
  double phot_autoaper[2] = {16., 16.}  ;
  if ( L.GlobVal().HasKey("PHOT_AUTOAPER_0") )
    phot_autoaper[0] = L.GlobVal().getDoubleValue("PHOT_AUTOAPER_0") ;
  if ( L.GlobVal().HasKey("PHOT_AUTOAPER_1") )
    phot_autoaper[1] = L.GlobVal().getDoubleValue("PHOT_AUTOAPER_1") ;
  cerr << "PHOT_AUTOAPERS : " << phot_autoaper[0] << " " << phot_autoaper[1]<< endl ;
  GlobVal().AddKey("PHOT_AUTOAPER_0",phot_autoaper[0]) ;
  GlobVal().AddKey("PHOT_AUTOAPER_1",phot_autoaper[1]) ;


  for (SEStarCIterator i = L.begin(); i != L.end(); ++i)
    push_back(new MultiMagSEStar(**i,phot_autoaper ));   

 
}



void Check(MultiMagSEStarList &LM, SEStarList &L)
{
  cerr << "Dans liste SE : " << endl ;
  for (SEStarIterator i = L.begin(); i != L.end() ; ++i)
    {
      SEStar *star  = *i;
      MultiMagSEStar *magstar = LM.FindClosest(star->x, star->y);
      if (magstar != NULL)
	{
	  double d2 = (star->x-magstar->x)*(star->x-magstar->x) + (star->y-magstar->y)*(star->y-magstar->y) ;
	  if (d2 > 1)
	    {
	      cerr << "MAL RETROUVEE " << endl ;
	      cerr << "SESTAR: " ;
	      star->dump();
	      cerr << "MAGSTAR: " ;
	      magstar->dump();
	    }
	}
      else
	{ 
	  star->dump();
	 cerr << "PAS RETROUVEE " << endl ; 
	}
    }


  cerr << "Dans liste MULTIMAGSE : " << endl ;
  for (MultiMagSEStarIterator i = LM.begin(); i != LM.end(); ++i)
    {
      MultiMagSEStar *magstar  = *i;
      SEStar *star = L.FindClosest(magstar->x, magstar->y);
      if (star != NULL)
	{
	  double d2 = (star->x-magstar->x)*(star->x-magstar->x) + (star->y-magstar->y)*(star->y-magstar->y) ;
	  if (d2 > 1)
	    {
	      cerr << "MAL RETROUVEE " << endl ;
	      cerr << "MAGSTAR: " ;
	      magstar->dump();
	      cerr << "SESTAR: " ;
	      star->dump();
	    }
	}
      else
	{
	  magstar->dump();
	 cerr << "PAS RETROUVEE " << endl ; 
	}
    }

}

void
MultiMagSEStarList::ComputeMag(int kbox, string band, double ZP, double eZP)
{
  for (MultiMagSEStarIterator mit = this->begin();mit != this->end();mit++)
    {
      (*mit)->ComputeMag(kbox, band, ZP, eZP);
    }
  return ;
}

bool
MultiMagSEStarList::UpDate_Assoc(SEStarList &L, string band)
{ 
  int nband_tot = GetNBand();
  char kkey[50] ;
  sprintf(kkey,"MAG_%d",nband_tot);
  GlobVal().AddKey(kkey,band);
    
  SetNBand(nband_tot + 1);
  // on rajoute une boite a chacun
  for (MultiMagSEStarIterator mit = this->begin(); mit != this->end(); ++mit)
    {
      MultiMagSEStar *magstar = *mit ;
      magstar->magboxes.push_back(ShortMagBox());
      ShortMagBox &mb = magstar->magboxes.back();
      CalibBox &cal = mb.calib;
      cal.band = band ;// pour info de provenance, peut changer apres avec la calib effectivement utilisee
	  // NB dans le cas ou le seeing est calcule par ailleurs et 
	  // stocke dans Fwhm()
      cal.seeing =  -1 ; 
      mb.flag = -1 ;
      mb.flagbad = -1 ;
      mb.fluxmax = -1 ;
      mb.f_auto = -1; // toujours cette vieille "convention" avec sextractorbox
      mb.ef_auto = -1; 
      mb.f_circ = -1 ;// doit etre calcule par ailleurs
      mb.ef_circ = -1;  
      mb.seeing = -1 ;
      mb.f_aper = -1;
      mb.ef_aper =  -1 ;
      mb.f_aper_other =  -1 ;
    }
  CountedRef<Gtransfo> identity = new GtransfoLin(); 
  BaseStarList *L1 = (BaseStarList*) this ; 
  BaseStarList *L2 = (BaseStarList*)  &L;
  double d_max =  0.01; // il faudra faire un cut dessus par la suite !
  StarMatchList *lmatch = ListMatchCollect(*L1, *L2, identity, d_max);
  double frac = lmatch->size()/(1.*size()) ;
  cerr << "Matching : " << lmatch->size() << " in " << size() << " " << frac << "%" << endl ;
  if ( frac < 0.90 )
    {
      cerr << "Not enough match " << endl ;
      return false ;
    }
  int nn = 0;
  for(StarMatchIterator it = lmatch->begin() ; it != lmatch->end(); it++)
    {
      MultiMagSEStar *mstar  = (*it).s1 ;
      SEStar *sestar = (*it).s2 ;
      ShortMagBox &mb = mstar->magboxes.back();
      CalibBox &cal = mb.calib;
      if ( nn < 5 )
	{
	  cerr << "Updating with assoc in band: " << cal.band << endl ;
	  nn++;
	}
      if (cal.band != band )
	{
	  cerr << "Mismatch in  band " << band << " " << cal.band << endl ;
	  return false ;
	}
      // NB dans le cas ou le seeing est calcule par ailleurs et 
      // stocke dans Fwhm()
      cal.seeing = sestar->Fwhm() ;
      mb.flag = sestar->Flag() ; 
      mb.flagbad = sestar->FlagBad() ; 
      mb.fluxmax = sestar->Fluxmax() ; 
      mb.f_auto = sestar->Flux_auto(); // toujours cette vieille "convention" avec sextractorbox   
      mb.ef_auto = sestar->Eflux_auto(); 

      mb.f_circ = sestar->Flux_circ_aper();// doit etre calcule par ailleurs
      mb.ef_circ = sestar->Eflux_circ_aper(); 

      mb.seeing = sestar->local_seeing ;
      mb.f_aper = sestar->aper_flux ;
      mb.ef_aper =  sestar->err_aper_flux ;
      mb.f_aper_other =  sestar->aper_flux_other ;
    }
  return true ;
}


bool
MultiMagSEStarList::UpDate(SEStarList &L, string band)
{
  if (size() != L.size())
    {
      cerr << "Taille des listes differentes dans Update " << " " << size() << " " << L.size() << endl ;
      return(UpDate_Assoc(L,band));
    } 
  int nband_tot = GetNBand();
  char kkey[50] ;
  sprintf(kkey,"MAG_%d",nband_tot);
  GlobVal().AddKey(kkey,band);
    
  SetNBand(nband_tot + 1);

 
  int nn = 0 ;
  MultiMagSEStarIterator mit = this->begin();
  for (SEStarCIterator i = L.begin(); i != L.end() && mit != this->end(); ++i, ++mit)
    {
      MultiMagSEStar *magstar = *mit ;
      const SEStar *star = *i ;
      double d2 = (magstar->x -star->x)*(magstar->x -star->x)+(magstar->y -star->y)*(magstar->y -star->y);
      if (d2 > 1.e-5)
	{	 
	  cerr << "etoiles differentes dans les listes " << band << endl ;
	  return(false);
	} 

      if ( nn < 5 )
	{
	  cerr << "Nbre de MagBoxes: " << magstar->magboxes.size() << endl ;
	}

      magstar->magboxes.push_back(ShortMagBox());
      ShortMagBox &mb = magstar->magboxes.back();
      CalibBox &cal = mb.calib;
      cal.band = band ;// pour info de provenance, peut changer apres avec la calib effectivement utilisee
	  // NB dans le cas ou le seeing est calcule par ailleurs et 
	  // stocke dans Fwhm()
      cal.seeing =  star->Fwhm(); 
      mb.flag = star->Flag() ; 
      mb.flagbad = star->FlagBad() ; 
      mb.fluxmax = star->Fluxmax() ; 
      mb.f_auto = star->Flux_auto(); // toujours cette vieille "convention" avec sextractorbox
      mb.ef_auto = star->Eflux_auto(); 
      mb.f_circ = star->Flux_circ_aper();// doit etre calcule par ailleurs
      mb.ef_circ = star->Eflux_circ_aper();

      mb.seeing = star->local_seeing ;
      mb.f_aper = star->aper_flux ;
      mb.ef_aper =  star->err_aper_flux ;
      mb.f_aper_other =  star->aper_flux_other ;

      if ( nn < 5 )
	{
	  cerr << "New Nbre de MagBoxes: " << magstar->magboxes.size() << endl ;
	  nn++;
	}
    }
  return(true);
}



MultiMagSEStar*  MultiMagSEStarList::FindEllipticClosest(double xx, double yy) const
{

  double min_dist2 = 1e30;
  const MultiMagSEStar *minstar = NULL;
  double dist2;
  for (MultiMagSEStarCIterator si = begin(); si!= end(); ++si) 
    { 
      const MultiMagSEStar *s = *si;
      dist2 =  s->SqEllipticDistance(xx, yy);
      if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;}
    }
  return (MultiMagSEStar *) minstar; // violates constness


}




void MultiMagSEStarList::check() const
{
  if (!empty())
    {
      unsigned nmag = front()->magboxes.size();
      for (MultiMagSEStarCIterator i = begin(); i !=end(); ++i)
	{
	  if ((*i)->magboxes.size() != nmag)
	    {
	      cerr << " writing a MultiMagSEStarList  with different number of apertures" << endl;
	      cerr << " hope you know what you are doing " << endl;
	      (*i)->dump(cerr);
	    }
	}
    }
}


int MultiMagSEStarList::write(const std::string &FileName) const
{
  check();
  return StarList<MultiMagSEStar>::write(FileName);
}



void 
MultiMagSEStarList::ComputeAlphaDelta(const FitsHeader & head)
{
  Gtransfo *wcs;
  WCSFromHeader(head, wcs);
  for (MultiMagSEStarIterator i = this->begin(); i != this->end(); ++i)
    {
      MultiMagSEStar *star = *i ;
      wcs->apply(star->x,star->y,star->alpha ,star->delta); 
    }
}

int MultiMagSEStarList::GetBandNumber(string band)
{
  if (GlobVal().HasKey(band) )
    return((int)(GlobVal().getDoubleValue(band)));
  else
    return(-1);
}

void MultiMagSEStarList::SetBandNumber(string band, int n)
{
  GlobVal().AddKey(band,n);
  return;
}

int MultiMagSEStarList::GetNBand()
{
  if (GlobVal().HasKey("NBAND") )
    return((int)(GlobVal().getDoubleValue("NBAND")));
  else
    return(0);
}
void MultiMagSEStarList::SetNBand(int n)
{
  if ( ! GlobVal().HasKey("NBAND") )
    GlobVal().AddKey("NBAND",n);
  else
    GlobVal().setDoubleValue("NBAND",n);
  return ;
}

int MultiMagSEStarList::MatchToOtherList(BaseStarList *l)
{
  CountedRef<Gtransfo> identity = new GtransfoLin(); 
  BaseStarList *L1 = (BaseStarList*) this ;
  double d_max =  1000.; // il faudra faire un cut dessus par la suite !
  StarMatchList *lmatch = ListMatchCollect(*L1, *l, identity, d_max);
  for(StarMatchIterator it = lmatch->begin() ; it != lmatch->end(); it++)
    {
      MultiMagSEStar *mstar  = (*it).s1 ;
      mstar->star = (*it).s2 ;
      mstar->star_dist = sqrt( (mstar->x-((*it).s2)->x)*(mstar->x-((*it).s2)->x)+(mstar->y-((*it).s2)->y)*(mstar->y-((*it).s2)->y) ) ;
    }
  return(lmatch->size());
}




/************ converters *************/
BaseStarList* MultiMagSE2Base(MultiMagSEStarList * This)
{ return (BaseStarList *) This;}

const BaseStarList* MultiMagSE2Base(const MultiMagSEStarList * This)
{ return (const BaseStarList *) This;}

BaseStarList& MultiMagSE2Base(MultiMagSEStarList &This)
{ return (BaseStarList &) This;}

const BaseStarList& MultiMagSE2Base(const MultiMagSEStarList &This)
{ return (const BaseStarList &) This;}

/* ================================== */

SEStarList* MultiMagSE2SE(MultiMagSEStarList * This)
{ return (SEStarList *) This;}

const SEStarList* MultiMagSE2SE(const MultiMagSEStarList * This)
{ return (const SEStarList *) This;}

SEStarList& MultiMagSE2SE(MultiMagSEStarList &This)
{ return (SEStarList &) This;}

const SEStarList& MultiMagSE2SE(const MultiMagSEStarList &This)
{ return (const SEStarList &) This;}

/* ================================== */
AperSEStarList* MultiMagSE2AperSE(MultiMagSEStarList * This)
{ return (AperSEStarList *) This;}

const AperSEStarList* MultiMagSE2AperSE(const MultiMagSEStarList * This)
{ return (const AperSEStarList *) This;}

AperSEStarList& MultiMagSE2AperSE(MultiMagSEStarList &This)
{ return (AperSEStarList &) This;}

const AperSEStarList& MultiMagSE2AperSE(const MultiMagSEStarList &This)
{ return (const AperSEStarList &) This;}


// on prend les neighbors plus proche que dist (pixels), mais on les ordonne
// par leur distance elliptique.
// ell_dist : distance elliptique en nombre de a 
// dist : distance en arc sec
void FindEllipticNeighb(double xsn, double ysn, double dist,
			MultiMagSEStarList & stlin, 
			MultiMagSEStarList & stl_neighb)
{
  //double n2lim = nrad_lim*nrad_lim;
  double d2lim = dist*dist;
  for(MultiMagSEStarIterator itse = stlin.begin(); itse!= stlin.end();itse++)
    {
      MultiMagSEStar *sestar = *itse;
      double n2 = sestar->SqEllipticDistance(xsn,ysn);
      double d2 = (xsn - sestar->x)*(xsn - sestar->x) +
		(ysn  - sestar->y)*(ysn  - sestar->y);
      if (d2 < d2lim)
	{
	  sestar->ell_dist = sqrt(n2) ; 
	  sestar->dist = sqrt(d2) ;
	  stl_neighb.push_back(sestar);
	}
    }
  stl_neighb.sort(&IncSqEllipticDist);
  return;
}

void CheckNeighbFinders(string name, double xsn, double ysn, double dist,
			MultiMagSEStarList & stlin)
{
  //double n2lim = nrad_lim*nrad_lim;
  double d2lim = dist*dist;
  MultiMagSEStarList  stl_neighb;
  for(MultiMagSEStarIterator itse = stlin.begin(); itse!= stlin.end();itse++)
    {
      MultiMagSEStar *sestar = *itse;
      double n2 = sestar->SqEllipticDistance(xsn,ysn);
      double d2 = (xsn - sestar->x)*(xsn - sestar->x) +
		(ysn  - sestar->y)*(ysn  - sestar->y);
      if (d2 < d2lim)
	{
	  sestar->ell_dist = sqrt(n2) ; 
	  sestar->dist = sqrt(d2) ;
	  stl_neighb.push_back(sestar);
	}
    }
  stl_neighb.sort(&IncSqEllipticDist);
  int ii =0;
  int ifirst=-1;
  double d_min_norm=1000.;
  for(MultiMagSEStarIterator itse = stl_neighb.begin(); itse!= stl_neighb.end();itse++)
    {
      MultiMagSEStar *neighb = *itse;
      double d_norm = neighb->NormalizedDistance(xsn,ysn);
      if (d_norm  < d_min_norm)
	{
	  ifirst=ii;
	  d_min_norm = d_norm  ;
	}
      cerr << name << "_" << ii << "  " << neighb->dist << " " << neighb->ell_dist << " " << d_norm << endl ;
      ii++;
    }
  if (ifirst != 0)
    cerr << "Normalized Distance : l'hote de " << name << " a change : " << ifirst << endl ;
  return;
}

void FindNormalizedDistNeighb(double xsn, double ysn, double dist,
			      MultiMagSEStarList & stlin, 
			      MultiMagSEStarList & stl_neighb)
{
  //double n2lim = nrad_lim*nrad_lim;
  double d2lim = dist*dist;
  for(MultiMagSEStarIterator itse = stlin.begin(); itse!= stlin.end();itse++)
    {
      MultiMagSEStar *sestar = *itse;
      double d2_ell = sestar->SqEllipticDistance(xsn,ysn);
      double d_norm = sestar->NormalizedDistance(xsn,ysn);
      if (d_norm < 0)
	cerr << "##Pb dans calcul NormalizedDistance" << endl ;
      double d2 = (xsn - sestar->x)*(xsn - sestar->x) +
		(ysn  - sestar->y)*(ysn  - sestar->y);
      if (d2 < d2lim)
	{
	  sestar->ell_dist = sqrt(d2_ell) ; 
	  sestar->norm_dist = d_norm ;
	  sestar->dist = sqrt(d2) ;
	  stl_neighb.push_back(sestar);
	}
    }
  stl_neighb.sort(&IncNormalizedDistance);
  return;
}






