#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
 
#include "hostabsorption.h"
#include "sampledfunction.h"
#include "fakesnia.h"
//#include "basestar.h"
#include "fakelist.h"
#define c 299792.458
#ifndef  WAVELENGTH_STEP
/* 10 Angst */
#define WAVELENGTH_STEP 10 
#endif

#ifdef STORAGE
//! object aperture data is incomplete or corrupted
#define APER_CORRUPT(m)    ((bool) ((m) & APER_COR_FLAG) )
// object isophotal data are incomplete or corrupted
#define ISO_CORRUPT(m)    ((bool) ((m) & ISO_COR_FLAG ) )
// memory overflow during deblending
#define MEM_OVFLOW1(m)    ((bool) ((m) & MEM_OV1_FLAG) )
// memory overflow during extraction
#define MEM_OVFLOW2(m)    ((bool) ((m) & MEM_OV2_FLAG) )
#endif


static string DefaultLCModel() 
{

  char *where = getenv("LCMODEL");
  if (where) return (where);
  else  return "";
}

int DecodeFormat(const char *FormatLine, const char *StarName)
{
if (!FormatLine || !StarName) return 0;
const char *p= strstr(FormatLine, StarName);
 if (!p) return  0;
return atoi( p + strlen(StarName));
}




FakeSNIa::FakeSNIa()
: Point(0.0,0.0)
{
  Set_to_Zero();
}
FakeSNIa::FakeSNIa(double xx, double yy)
{
 Set_to_Zero();
 x=xx;
 y=yy;
}



FakeSNIa::FakeSNIa(double xx, double yy, double RS,double day,double istretch,double iextinction, double
icolor, double idispersion, double ialpha, double ibeta, double irv, double iMW_rv, double iMW_extinction)
  : Point(xx,yy)
{  
redshift=RS;
d0=day;
stretch=istretch;
h_extinction=iextinction;
color= icolor;
dispersion=idispersion;
alpha=ialpha;
beta=ibeta;
h_rv=irv;
MW_rv = iMW_rv;
MW_extinction = iMW_extinction;

h_ra=0.0;
h_dec=0.0;
h_imag =0.0;



}

void
FakeSNIa::Set_to_Zero()
{
  x=0.0;
  y=0.0;
  d0=0.0;
  stretch=1.0;
  h_extinction=0.0;
  redshift=0.0;
  color= 0.0;
  dispersion=0.0;
alpha=0.0;
beta=0.0;
h_rv=0.0;
h_ra=0.0;
h_dec=0.0;
h_imag =0.0;
MW_rv = 0.0;
MW_extinction = 0.0;


}

double
FakeSNIa::Magnitude(const GeneralCosmo &cosmo_model,Instrument &instrum, const
string &band,const double day, const string magsys, const double H0)
{


 

  double d0 = 10; // pc, distance for absolute magnitude
  



 
  
  SaltModel* sn_model = new SaltModel();
  LcParam* m_pday = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("DayMax"));
  LcParam* m_stretch = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("Stretch"));
  LcParam* m_dflux = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("Dflux"));
  LcParam* m_ext = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("Tau"));
  LcParam* m_color = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("Color"));
  LcParam* m_rv = const_cast<LcParam *>(sn_model->params.LocateOrAddParam("Rv"));
  
  	m_ext->val = 0.0; 
	m_pday->val = .0;
	m_stretch->val = 1;  
	m_color->val =0.0;
	m_rv->val=0.0;
  


  
  const Filter instrumental_filter = instrum.EffectiveFilterByBand(band);
  
  // calcul du point zero
  const Sampled1DFunction &refstar4magsystem = RefSpectrumForMags(magsys);
  double ref_flux = IntegPhotons(refstar4magsystem,instrumental_filter,0.);
  double zp = 2.5*log10(ref_flux);
  
  double flux_at_dl,mag,dl;
    
	sn_model->SetRedshift(Redshift());	
	m_ext->val = H_Extinction();
	m_pday->val = D0();
	m_stretch->val = Stretch();  
	m_color->val = Color();		 
	m_rv->val=H_Rv();
    
    dl = cosmo_model.Dl(redshift); // (en unites de c/H0)
    dl  = (c/H0)*dl; // Mpc 
    dl *= 1.0E6; // pc
    
    m_dflux->val = pow(d0/dl,2);

    /********************************************compute filter with MW absorption begin**************************************/
    
    
  Filter EffectiveFilterWithMilkyWay = instrumental_filter;
  if(fabs(MW_Extinction())>1.e-10) 
  {
    double data[2];
    data[0] = MW_Extinction();
    data[1] = MW_Rv(); 
    // mean Rv for the Galaxy (see E.L. Fitzpatrick astro-ph/9809387 , J.A. Cardelli et al., ApJ 425, 245 (1989) )
    Computed1DFunction HostGalacticExtinction(&FluxExtinctionCardelli,data);
    (General1DFunction &) EffectiveFilterWithMilkyWay = Product(WAVELENGTH_STEP,2,&instrumental_filter,&HostGalacticExtinction);   
   } 
    
    /********************************************compute filter with MW absorption end**************************************/
        // facteur 1/(1+z) ????
    flux_at_dl = sn_model->Flux(&EffectiveFilterWithMilkyWay,day)/(1+Redshift());
    mag = zp-2.5*log10(flux_at_dl);
 
	
	mag = mag + dispersion -alpha*(stretch-1) + beta*color;
	return(mag);	
	
}





void
FakeSNIa::dumpn(ostream& s) const
{
  s << " Ra : " << x ;
  s << " Dec : " << y ;
  s << " redshift : " << Redshift ();
  s << " stretch : " << Stretch() ;
  s << " day : " << D0() ;
  s << " color_term : " << Color() ;
  s << " Dispersion_term : " << Dispersion() ;
  s << " alpha : " << Alpha();
  s << " beta : " << Beta();
  
  s << " MW Rv : " << MW_Rv();
  s << " MW Extinction : " << MW_Extinction() ;    
  
  s << " Host Rv : " << H_Rv();
  s << " Host Extinction : " << H_Extinction() ;
  s << " Host Ra : " << H_Ra();
  s << " Host Dec : " << H_Dec(); 
  s << " Host Imag : " << H_Imag();


}
 


void
FakeSNIa::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}

void
FakeSNIa::write(ostream& s)  const
{
	writen(s);
	s << endl;
}

void
FakeSNIa::writen(ostream& s)  const
{

  s << Ra() << " ";
  s << Dec() << " ";
  s << Redshift() << " ";
  s << Stretch() <<" ";
  s <<  D0() <<" " ;
  s << Color() <<" ";
  s << Dispersion() <<" ";
  s << Alpha() <<" ";
  s << Beta() <<" ";
  
  s << MW_Rv() <<" ";
  s << MW_Extinction() <<" ";  
  
  s << H_Rv() <<" ";
  s << H_Extinction() <<" ";  
  s << H_Ra() <<" ";
  s << H_Dec() <<" "; 
  s << H_Imag() <<" ";


  
}


void
FakeSNIa::read_it(istream& s, const char * Format)
{
  int format = DecodeFormat(Format, "FakeSNIa");
    s >> Ra();
    s >> Dec();
    s >> Redshift();
    s >> Stretch();
    s >>  D0();
    s >> Color();
    s >> Dispersion();
    s >> Alpha();
    s >>  Beta();

if (format >2)
{
  s >> MW_Rv();
  s >> MW_Extinction();

  s >> H_Rv();
  s >> H_Extinction();
  s >> H_Ra();
  s >> H_Dec(); 
  s >> H_Imag();

}
 
  return ;
}

FakeSNIa*  FakeSNIa::read(istream& r, const char *Format)
{
  FakeSNIa *pstar = new FakeSNIa();  
  pstar->read_it(r, Format);
  return(pstar);
}


void FakeSNIa::WriteHeader(ostream & stream) const
{
	string format =WriteHeader_(stream);
	stream <<"# format " << format << endl;
	stream << "# end" << endl;
	
}

std::string FakeSNIa::WriteHeader_(ostream & pr, const char *i) const
{
  if (i== NULL) i= "";
  
  pr    << "# Ra"<< i <<" : " << endl 
  	<< "# Dec"<< i <<" : " << endl 
	<< "# Redshift"<< i <<" : Redshift of the fake object" << endl
	<< "# Stretch"<< i <<" : Stretch " << endl
	<< "# D0"<< i <<" : Day of max luminosity" << endl
    	<< "# Color"<< i <<" : Use to compute the polynome correction" << endl
	<< "# Dispersion_term"<< i <<" : mag effective dispersion" << endl
	<< "# Alpha"<< i <<" : strectch mag correction factor (-aplha*(s-1))" << endl
	<< "# Beta"<< i <<" : color mag correction factor (+Beta*Color)" << endl
	
	
	<< "# MW_Rv"<< i <<" :  MW Cardelli factor" << endl
	<< "# MW_Extinction"<< i << " : MW E(B-V) use to compute MW Cardelli'extinction"<< endl	
	
	
	<< "# H_Rv"<< i <<" :  host Cardelli factor" << endl
	<< "# H_Extinction"<< i << " : Host E(B-V) use to compute host Cardelli'extinction"<< endl	
	<< "# H_Ra"<< i <<" : Host Ra" << endl
	<< "# H_Dec"<< i <<" : Host Dec" << endl	
	<< "# H_Imag"<< i <<" : Host i mag" << endl;	
	
/* 1 is the current format id for FakeSNIas (when being written) it must correspond
to the right behaviour of the read routine ( and match what write does ! ) */
return " FakeSNIa 3 ";
}






#ifdef USE_ROOT
ClassImp(FakeSNIa)

/* To Generate the sestardict.cc file :
LINKDEF_CONTENT : #pragma link C++ class FakeSNIa+;
*/
#endif /* USE_ROOT */


//********************   FINDEFINITION FakeSNIa   *********************


#include "fakelist.cc" /* since starlist is a template class */


#ifdef USE_ROOT
template class FakeListWithRoot<FakeSNIa>;
ClassImpT(FakeListWithRoot,FakeSNIa);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class CountedRef<FakeSNIa>-;
LINKDEF_CONTENT : #pragma link C++ class list<CountedRef<FakeSNIa> >;
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeSNIa> >::unique();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeSNIa> >::sort();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeSNIa> >::merge(list <CountedRef<FakeSNIa> >)&;
LINKDEF_CONTENT : #pragma link C++ class StarList<FakeSNIa>-;
LINKDEF_CONTENT : ostream& operator << (ostream&, const StarList<FakeSNIa>&);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<FakeSNIa>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<FakeSNIa>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<FakeSNIa>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef FakeSNIaIterator;
*/
#include "root_dict/sestardict.cc"
#endif /* USE_ROOT */

template class FakeList<FakeSNIa>; // because StarListWithRoot<> derives from StarList<>


