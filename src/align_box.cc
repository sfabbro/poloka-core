#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>


#include "sestar.h"
#include "starlist.h"
#include "align_box.h"
#include "datacards.h"
#include "frame.h"
#include "gtransfo.h"
#include "listmatch.h"


#include "starmatch.h"

#ifdef STORAGE

DatAlign::DatAlign() {
#ifdef STORAGE
    prlevel = 0 ;
  marge = 0 ;
  //  flag = 0 ;
  saturation  = 0. ;
  b_min = 0 ;
  degre = 0 ;
  type_align_phot = 0 ;
  ph_fluxmax_prct = 0 ;
  ph_fluxmax_N = 0 ;
  n_sigma = 0. ;
  lim_n1n2 = 0. ;
#endif
}

void
DatAlign::Saturation(const double saturlevel)
{
  saturation = saturlevel;
}

void
DatAlign::Default(const double saturlevel)
{
  saturation  = saturlevel  ;
  ph_fluxmax_prct = 0.3 ;
  ph_fluxmax_N = 200 ;
}


void
DatAlign::Print(ostream & s) const 
{
  s  << "*****DatAlign******" << endl ;
  s << "saturation: " << saturation<< endl ;
  s << "fluxmax_prct: " <<  ph_fluxmax_prct << endl ;
  s << "fluxmax_N: " <<   ph_fluxmax_N<< endl ;
  s  << endl ;
}




void 
DatAlign::LitDataCard(DataCards & data)
{
  ph_fluxmax_prct =  data.DParam("ALIGN_P_PIXMAX_PRCT") ;
  ph_fluxmax_N =  data.IParam("ALIGN_P_PIXMAX_N") ;
}

bool 
SM_DecFluxMax(const StarMatch & S1, const StarMatch & S2)
{
  const SEStar * p1 = S1.s1;
  const SEStar * p2 = S2.s1;
  bool test = DecFluxMax(p1,p2);
  return(test) ;
}


int
GoodForPhotomAlign(DatAlign const & dat, 
		   StarMatchList & liste)
{
  return GoodForPhotomAlign(liste, dat.saturation, 
			    dat.ph_fluxmax_prct, dat.ph_fluxmax_N);
}




int
GoodForPhotomAlign(StarMatchList & liste, double saturation, 
		   double prctage, int N_max)
{
  // on enleve les saturees, on prend les xx% plus brillantes
  // en brillance de surface.
  int Ngardees = 0 , Nok = 0;
  if (!getenv("NOFILTERING")) 
    for (StarMatchIterator it= liste.begin(); it!= liste.end(); )
      {
	StarMatch starm = *it ;
	
	
	SEStar * pstar1 = (SEStar *) starm.s1;
	SEStar * pstar2 = (SEStar *) starm.s2;
	
	if ((pstar1->FlagBad() == 0 ) &&
	    ( pstar1->Flag() == 0 ) &&
	    ( pstar1->Fluxmax() < saturation ) && 
	    (pstar2->FlagBad() == 0 ) &&
	    ( pstar2->Flag() == 0 ) &&
	    ( pstar2->Fluxmax() < saturation ))
	  {
	    it++ ; Nok++ ; ;
	  }
	else
	  {
	    it = liste.erase(it) ;
	  }
      }
  else
    Nok = liste.size();
  liste.sort(&SM_DecFluxMax);

  Ngardees = (int) (prctage*1.0*Nok);

  cout << Ngardees << " plus brillantes pour alignement photom."<< endl ;

  // si dat.fluxmax_N > 0, 
  // on prend le min de dat.fluxmax_N et  Ngardees (on est alors sur que < Nok1)
  if (( N_max > 0 ) && (Ngardees > N_max))
      Ngardees = N_max ;

  int count = 0 ;
  StarMatchIterator it ;
  for (it= liste.begin(); (it!= liste.end()) && ( count < Ngardees) ;++it )
    {      	
      count++ ;
    }
  liste.erase(it,liste.end()) ;
  
  return(Ngardees);

}


#ifdef IS_IT_USEFUL

// cut des out lier a n sigma une fois une solution trouvee

int
CutForPhotomAlign(Poly & pl, double n_sigma, 
		   StarMatchList & liste)
{
  // on enleve les saturees, on prend les xx% plus brillantes
  // en brillance de surface.
  int Ngardees = 0 ;
  double n2sig = n_sigma*n_sigma ;
  for (StarMatchIterator it= liste.begin(); it!= liste.end(); )
    {
      StarMatch starm = *it ;


      SEStar * pstar1 = (SEStar *) starm.s1;
      SEStar * pstar2 = (SEStar *) starm.s2;

      double dif =  pstar2->flux - pl(pstar1->flux) ;
      dif = dif *dif ;
      double err2 =  pstar1->EFlux()*pstar1->EFlux()+pstar2->EFlux()*pstar2->EFlux() ;

      if ((err2 > 1.e-10) && ((dif/err2)<n2sig))
	{
	  it++ ; Ngardees++ ; ;
	}
      else
	{
	  it = liste.erase(it) ;
	}
    }
  return(Ngardees);

}



// version simplifiee pour seb

void
SearchPhotomAlign(double & A, double & B, const SEStarList & stl1,
		  const SEStarList & stl2, const Gtransfo * tf,
		  const double saturation, const double Nmax)
{

  BaseStarList *l1 = (BaseStarList*)  &stl1;
  BaseStarList *l2 = (BaseStarList*)  &stl2;
  
  double dist_max = 2. ;
  StarMatchList *matchlist = ListMatchCollect(*l1, *l2, tf, dist_max) ;
 
  double prctage = 1. ;
  int Ngood = GoodForPhotomAlign(*matchlist, saturation,  prctage, Nmax);
  cout << " Selected " << Ngood << " stars for photometric alignment " << endl ;
  int degre = 1 ;
  double chi2 ;

  Poly pl = PhotomAlign_(*matchlist, degre, chi2);
  if (Ngood > 0 ) cout << " First pass Chi2 Fit : " << chi2/Ngood ;

  int n_sigma = 3 ;
  int N1  = CutForPhotomAlign(pl, n_sigma, *matchlist) ;
  cout << " Left with " << N1 << " stars for a second pass " << endl ;
  pl = PhotomAlign_(*matchlist, degre, chi2);
  if (N1 > 0 ) cout << " Second pass Chi2 Fit : " << chi2/N1 << endl;

  A = pl[1] ;
  B = pl[0] ;

}


Poly
SearchPhotomAlign(DatAlign & datal, SortieAlign & sortie, 
		  SEStarList & stl1,
		  SEStarList & stl2, Gtransfo * tf)
{
  BaseStarList *l1 = (BaseStarList*)  &stl1;
  BaseStarList *l2 = (BaseStarList*)  &stl2;

  double dist_max = 2. ;
  StarMatchList *liste = ListMatchCollect(*l1, *l2, tf, 
					       dist_max) ;

  Poly pl(datal.degre); // il est a zero
  // 1ere passe
  sortie.n_1 = GoodForPhotomAlign(datal, *liste);
  cout << sortie.n_1 << " pour alignement photometrique " << endl ;

  if ( sortie.n_1 < NLIM_ALIGN ) 
    {
      cerr << "Number of selected stars for photom alignment too small: " 
	   << sortie.n_1 << endl ;

      return(pl);
    }

  if (datal.prlevel > 0 )
    {
      ofstream pr("match_phot1.list");
      liste->write(pr);
      pr.close();

    }

  // on fitte Y = a X 
  if ( datal.type_align_phot == 1)
    {
      cout << "###### test phot " << endl ;
      Poly  plu = PhotomAlign2_(*liste) ;
 
      return(plu);
    } 
 
  // default : procedure ou on fitte un polynome


  Poly pl1 = PhotomAlign_(*liste, datal.degre, sortie.chi2_1);
  if ( datal.prlevel > 0 )
    {
      cout << "Poly 1: " << pl1 << endl ;
      cout << "chi2 1: " << sortie.chi2_1 << endl ;
    }
      

  if ( sortie.chi2_1 < 0 ) 
    {
      cerr << "fit failed for photometric alignement " << endl ;
      return(pl);
    }
  

  pl = pl1 ;


  sortie.chi2_1 =  sortie.chi2_1 /  (1.*sortie.n_1) ;

  //2eme passe
  sortie.n_sigma = datal.n_sigma ;
  sortie.n_2 = CutForPhotomAlign(pl1, datal.n_sigma, *liste) ;

  if ( ( sortie.n_2 > 0 ) && (sortie.n_2 > datal.lim_n1n2 * sortie.n_1) )
    {
      if (datal.prlevel > 0 )
	{
	  ofstream pr("match_phot2.list");
	  liste->write(pr);
	  pr.close();

	}
      Poly pl2 = PhotomAlign_(*liste, datal.degre, sortie.chi2_2);
      if ( datal.prlevel > 0 )
	{
	  cout << "Poly 2: " << pl2 << endl ;
	  cout << "chi2 2: " << sortie.chi2_2 << endl ;
	}
      sortie.chi2_2 =  sortie.chi2_2 /  (1.*sortie.n_2) ;

      if ( (sortie.chi2_2>0) && ( sortie.chi2_2 < sortie.chi2_1 ) )
	{	 
	  pl = pl2;

	  // 3eme passe

	  sortie.n_3 = CutForPhotomAlign(pl2, datal.n_sigma, *liste) ;


	  if ( ( sortie.n_3 > 0 ) && (sortie.n_3 > datal.lim_n1n2 * sortie.n_3) )
	    {
	      if (datal.prlevel > 0 )
		{
		  ofstream pr("match_phot3.list");
		  liste->write(pr);
		  pr.close();

		}
	      Poly pl3 = PhotomAlign_(*liste, datal.degre, sortie.chi2_3);
	      if ( datal.prlevel > 0 )
		{
		  cout << "Poly 3: " << pl3 << endl ;
		  cout << "chi2 3: " << sortie.chi2_3 << endl ;
		}
	      sortie.chi2_3 =  sortie.chi2_3 /  (1.*sortie.n_3) ;

	      if ( (sortie.chi2_3>0) && ( sortie.chi2_3 < sortie.chi2_2 ) )
		{	 
		  pl = pl3;
		}
      

	    }
	}

    }

  return(pl);
}

int MeanSigAl(double *mag,double *err,int *ndata,
	      double nsigma,double *mean,double *sigma)
{
int deb= 1;
int pass,passmx,i,n=0;
double m,s2,v,scut;

*mean = *sigma = 0.;
passmx = ( nsigma <= 0. ) ? 1 : 2 ;
for (pass=1;pass<=passmx;pass++) 
  {
    m = s2 = 0.;
    n=0;
    scut = 1.e+35;
    if( pass == 2 ) 
      scut=nsigma* *sigma;
    for (i=0;i<*ndata;i++) 
      {
	v = *(mag+i);
	if( *(err+i) > 0. && fabs(v-*mean) < scut ) 
	  {
	    n++;
	    m += v;
	    s2 += v * v;
	  } 
      }
    if ( n >= 2 ) 
      {
	*mean = m / n;
	*sigma = sqrt( s2 / n - m/n * m/n );
      } 
    else 
      {
	*mean = *sigma = 0.;
	*ndata=n;
	return(-1);
      }
    if ( deb>0 ) printf("MeanSig: pass=%d mean=%f sigma=%f n=%d\n"
			,pass,*mean,*sigma,n);
  }
 *ndata=n;
 return(0);
}


Poly 
PhotomAlign2_(StarMatchList & liste)
{
  int N = liste.size() ;
  cout << "Nbre etoiles pour fit " << N << endl ;
  Poly pl(1); 
 double *rapp = new double[N];
 double *err = new double[N];
 int ii = 0 ;


 for(StarMatchCIterator it = liste.begin(); it != liste.end() ; it++)
   {    
     StarMatch starm = *it ;
     SEStar * pstar1 = (SEStar *) starm.s1;
     SEStar * pstar2 = (SEStar *) starm.s2;

     double flux1 = pstar1->flux ;
     double flux2 = pstar2->flux ;
     if ( flux1 >1.e-10 ) 
       {
	 rapp[ii] = flux2/flux1 ;   
	 err[ii] = 1. ;
	 ii++ ;
       }
   }

 double coeff = 1., sigma = 1.  ;
 int ndata = ii ;
 MeanSigAl(rapp, err, &ndata, 2. , &coeff, &sigma ) ;


 delete [] err;
 delete [] rapp;




 pl[1] = coeff ;
 pl[0] = 0. ;
 return(pl);
}








Poly 
PhotomAlign_(StarMatchList & liste, int degre, double & chi2)
{
  int N = liste.size() ;
  cout << "Nbre etoiles pour fit " << N << endl ;
  Poly pl(degre); 
 double *flux1 = new double[N];
 double *flux2 = new double[N];
 double * eflux2 = new double[N] ;
 int ii = 0 ;


 for(StarMatchCIterator it = liste.begin(); it != liste.end() ; it++)
   {    
     StarMatch starm = *it ;
     SEStar * pstar1 = (SEStar *) starm.s1;
     SEStar * pstar2 = (SEStar *) starm.s2;

     flux1[ii] = pstar1->flux ;
     flux2[ii] = pstar2->flux ;
     //eflux2[ii] = pstar2->EFlux() ;
     eflux2[ii] = pstar1->EFlux()*pstar1->EFlux()+pstar2->EFlux()*pstar2->EFlux() ;
     //eflux2[ii] = sqrt(flux2[ii]) ;
     
     

     ii++ ;
   }
 cout << " degre " << degre << endl ;
 Vector VFlux1(N,flux1 );
 Vector VFlux2(N,flux2 );
 Vector EVFlux2(N,eflux2 );
 Vector ErrCoef(degre);  
 chi2 = -1 ;     
 TRY{
   chi2 = pl.Fit(VFlux1,VFlux2,EVFlux2,degre,ErrCoef);}
   CATCH(i){
     cerr << "Erreur Fit " << endl ;
     if ( i == sizeMismatchErr)
       {
	 cerr << " Erreur Fit: sizeMismatchErr " << endl ;
       }
     if ( i == singMatxErr)
       {
	 cerr << " Erreur Fit: singMatxErr " << endl ;
       }
     THROW(i);
   }
   ENDTRY 

  cout << "chi2 : " << chi2 << endl ;

 delete [] flux1 ;
 delete [] flux2 ;
 delete [] eflux2 ;
 return(pl);
}

// On applique une transfo photom

void 
Poly_Apply(Poly const & pl, SEStarList & stl)
{
 
  for (SEStarCIterator it= stl.begin(); it!=stl.end(); it++)
    {
      SEStar * pstar = (SEStar *) *it ;
      pstar->flux *= pl[1];
      /*pstar->Flux_aper() *= pl[1];
      pstar->Flux_iso() *= pl[1];
      pstar->Flux_isocor() *= pl[1];*/
    }


}
// On applique une transfo photom

void 
Poly_Apply(Poly const & pl, Image & img)
{
 int tx = img.Nx();
 int  ty = img.Ny();

  Pixel  *p = img.begin() ;

  for(int i=0; i < ty ; i++ )
      { 
        for(int j=0; j < tx  ; j++ )
	  { 
	    *p = pl(*p);
	    p++;
	  }
      }
}

void 
Permutation_XY(SEStarList & stl)
{
 
  for (SEStarCIterator it= stl.begin(); it!=stl.end(); it++)
    {
      SEStar * pstar = (SEStar *) *it ;
      double xi = pstar->x ;
      pstar->x = pstar->y ;
      pstar->y = xi ;
    }
return;
      
}




StarMatchList *
SimpleAlignment( BaseStarList *lr1, BaseStarList *lr2, double distance)
{
  GtransfoLin ident;
  StarMatchList *lma = ListMatchCollect(*lr1, *lr2, &ident, 
					distance) ; 
  cout << lma->size() << " etoiles matchees " << endl ;
  return(lma);
}

#endif /* IS_IT_USEFUL */


/* discards stars which are too close from sides of the frame */
void SEStarListCutEdges(SEStarList *L, const Frame &F, const double MinDist)
{
  for (SEStarIterator si = L->begin(); si != L->end(); )
    {
      if (F.MinDistToEdges(**si) < MinDist) 
	{
	  si = L->erase(si); 
	}
      else ++si;
    }
}

#endif


