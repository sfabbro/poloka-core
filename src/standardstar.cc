#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#include "standardstar.h"
#include "usnoutils.h"
#include "wcsutils.h"
#include "frame.h"
#include "starmatch.h"
#include "listmatch.h"

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif


// Converter :
BaseStarList* Standard2Base(StandardStarList * This)
{ return (BaseStarList*) This;}

const BaseStarList* Standard2Base(const StandardStarList * This)
{ return (BaseStarList*) This;}

StandardStar::StandardStar()
: BaseStar(0.,0.,0.)
{
  Set_to_Zero();
}

StandardStar::StandardStar(double xx, double yy, double ff)
: BaseStar(xx,yy,ff)
{
  Set_to_Zero();
}

void StandardStar::Set_to_Zero()
{
  name="";
  airmass=0.0;
  fluxpersec=0.0;
  efluxpersec=0.0;
  color=NONE;
  ra="";
  dec="";
  vmag=0.0;
  bvmag=0.0;
  ubmag=0.0;
  vrmag=0.0;
  rimag=0.0;
  vimag=0.0;
  n=0;
  m=0;
  dvmag=0.0;
  dbvmag=0.0;
  dubmag=0.0;
  dvrmag=0.0;
  drimag=0.0;
  dvimag=0.0;
}

void StandardStar::Read(istream& r, const char *Format)
{
  //  char Name[16];
  //  char Ra[12];
  //  char Dec[12];
  int format = DecodeFormat(Format, "StandardStar");

  if (format > 0)
    {
      r >> x;
      r >> y;
      r >> flux;
      //      r.get(Name,13);
      //      for (int u=0; u<13; u++) Name[u]=Name[u+1];
      //      name = string(Name);
      r >> name;
      r >> ra;
      //      r.get(Dec,11);
      //      dec = string(Dec);
      r >> dec;
    }
  else
    {
      //      r.get(Name,8);
      //      name = string(Name);

      r >> name;

      //      r.get(Ra,13);
      //      ra = string(Ra);
      //      ra[2]=':';ra[5]=':';
      r >> ra;

      //      r.get(Dec,13);
      //      dec = string(Dec);
      r >> dec;
    }
  r >> vmag;
  r >> bvmag;
  r >> ubmag;
  r >> vrmag;
  r >> rimag;
  r >> vimag;
  r >> n;
  r >> m;
  r >> dvmag;
  r >> dbvmag;
  r >> dubmag;
  r >> dvrmag;
  r >> drimag;
  r >> dvimag;
  if (format == 0)
    {
      x = RaStringToDeg(ra);
      y = DecStringToDeg(dec);
      flux = 0.0;
    }
}


StandardStar* StandardStar::read(istream& r, const char *Format)
{
  StandardStar *pstar = new StandardStar();
  pstar->Read(r,Format);
  return(pstar);
}

void StandardStar::dumpn(ostream& s) const
{
  s << " x : " << x;
  s << " y : " << y;
  s << " flux : " << flux;
  s << " name : " << name;
  s << " ra : " << ra;
  s << " dec : " << dec;
  s << " magV : " << vmag;
  s << " magBV : " << bvmag;
  s << " magUB : " << ubmag;
  s << " magVR : " << vrmag;
  s << " magRI : " << rimag;
  s << " magVI : " << vimag;
  s << " n : " << n;
  s << " m : " << m;
  s << " dmagV : " << dvmag;
  s << " dmagBV : " << dbvmag;
  s << " dmagUB : " << dubmag;
  s << " dmagVR : " << dvrmag;
  s << " dmagRI : " << drimag;
  s << " dmagVI : " << dvimag;
}

void StandardStar::dump(ostream& s) const
{
  dumpn(s);
  s << endl;
}

void StandardStar::writen(ostream& s) const
{
  s << x << " " ;
  s << y << " " ;
  s << flux << " " ;
  s << name << " " ;
  s << ra << " " ;
  s << dec << " " ;
  s << vmag << " " ;
  s << bvmag << " " ;
  s << ubmag << " " ;
  s << vrmag << " " ;
  s << rimag << " " ;
  s << vimag << " " ;
  s << n << " " ;
  s << m << " " ;
  s << dvmag << " " ;
  s << dbvmag << " " ;
  s << dubmag << " " ;
  s << dvrmag << " " ;
  s << drimag << " " ;
  s << dvimag << " " ;
}

void StandardStar::write(ostream& s) const
{
  writen(s);
  s << endl;
}

string StandardStar::WriteHeader_(ostream & pr, const char *i) const
{
  if (i==NULL) i="";
  string baseStarFormat =  BaseStar::WriteHeader_(pr, i);
  pr << "# name"<< i << " : name of the standard star " << endl
     << "# ra"<< i << " : ra... " << endl
     << "# dec"<< i << " : dec... " << endl
     << "# magV"<< i << " : std mag in V " << endl
     << "# magBV"<< i << " : std(magB-magV) " << endl
     << "# magUB"<< i << " : std(magU-magB) " << endl
     << "# magVR"<< i << " : std(magV-magR) " << endl
     << "# magRI"<< i << " : std(magR-magI) " << endl
     << "# magVI"<< i << " : std(magV-magI) " << endl
     << "# n"<< i << " : cf Landolt paper " << endl
     << "# m"<< i << " : cf Landolt paper " << endl
     << "# dmagV"<< i << " : error on magV " << endl
     << "# dmagBV"<< i << " : error on magBV " << endl
     << "# dmagUB"<< i << " : error on magUB " << endl
     << "# dmagVR"<< i << " : error on magVR " << endl
     << "# dmagRI"<< i << " : error on magRI " << endl
     << "# dmagVI"<< i << " : error on magVI " << endl;

  //static char format[256];
 string format = baseStarFormat + " StandardStar 1";
 //sprintf(format,"%s StandardStar %d",baseStarFormat, 1);
return format;
}

double StandardStar::Magnitude(StandardColor couleur)
{
  if (couleur == VBAND) return Vmag();
  if (couleur == BBAND) return Bmag();
  if (couleur == UBAND) return Umag();
  if (couleur == RBAND) return Rmag();
  if (couleur == IBAND) return Imag();
  cout << " " << (char *) couleur << " is not a referenced color. Use : V,B,U,R or I. "<< endl ;
  return 0.0;
}

double StandardStar::Magnitude()
{
  return this->Magnitude(this->color);
}

double StandardStar::eMagnitude(StandardColor couleur)
{
  if (couleur == VBAND) return dVmag();
  if (couleur == BBAND) return dBmag();
  if (couleur == UBAND) return dUmag();
  if (couleur == RBAND) return dRmag();
  if (couleur == IBAND) return dImag();
  cout << " " << (char *) couleur << " is not a referenced color. Use : V,B,U,R or I. "<< endl ;
  return 0.0;
}

double StandardStar::eMagnitude()
{
  return this->eMagnitude(this->color);
}

double StandardStar::DeltaColor()
{
  if (color == IBAND) return VImag();
  //  if (color == IBAND) return RImag();
  if (color == RBAND) return VRmag();
  if (color == VBAND) return BVmag();
  if (color == BBAND) return UBmag();
  if (color == UBAND) return UBmag();
  cout <<"You asked for the color of the star and there is no band specified !!!!!!!!!!!"<< endl;
  return 0.0;
}

StandardColor GetColor (const FitsHeader &header)
{
  string scolor = header.KeyVal("TOADBAND");
  if (scolor=="V") return VBAND;
  if (scolor=="B") return BBAND;
  if (scolor=="U") return UBAND;
  if (scolor=="R") return RBAND;
  if (scolor=="I") return IBAND;

  cout << scolor << " doesn't have any correspondance in the standard stars filters. Only : V,B,U,R or I." << endl;
  return NONE;
}

static double AirmassTerm (const FitsHeader &header)
{

  string color = header.KeyVal("TOADBAND");
  string telName = header.KeyVal("TOADINST");
  if( (telName=="CFH12K") || (telName=="12k") || (telName=="Cfht12K"))
    {
      if (color=="B") return 0.17;
      if (color=="V") return 0.10;
      if (color=="R") return 0.06;
      if (color=="I") return 0.05;
    }
  return 0.0;
}


static double ColorTerm (const FitsHeader &header)
{

  string color = header.KeyVal("TOADBAND");
  string telName = header.KeyVal("TOADINST");
  if( (telName=="CFH12K") || (telName=="12k") || (telName=="Cfht12K"))
    {
      if (color=="B") return 0.057;
      if (color=="V") return 0.005;
      if (color=="R") return 0.04;
      if (color=="I") return -0.01;
    }
  return 0.0;
}

StandardStarList* GetSelectedStandardStarList(const FitsHeader &header)
{
  Frame W(header);
  GtransfoLin Pix2RaDec;
  WCSLinTransfoFromHeader(header, Pix2RaDec);
  Frame Wradec = W.ApplyTransfo(Pix2RaDec);

  StandardColor couleur = GetColor(header);

  double RaMin = Wradec.xMin; 
  double RaMax = Wradec.xMax;
  double DecMin= Wradec.yMin;
  double DecMax= Wradec.yMax;
  if (Wradec.xMin > Wradec.xMax)
    {
      RaMin = Wradec.xMax; 
      RaMax = Wradec.xMin;
    }
  if (Wradec.yMin > Wradec.yMax)
    {
      DecMin = Wradec.yMax; 
      DecMax = Wradec.yMin;
    }

  char *env_var = getenv("STANDARDFILE");
  if (!env_var)
    {
      cerr << " you should define STANDARDFILE env var to run this code " << endl;
      return NULL;
    }
  string standardfile(env_var);

  StandardStarList *stdstarlist = new StandardStarList(standardfile);

  if (!stdstarlist || stdstarlist->size()==0)
    {
      cout << "Bad standard file !!!" << endl;
      return NULL;
    }

  for (StandardStarIterator si = stdstarlist->begin(); si != stdstarlist->end(); )
    {
      StandardStar * pstar = (StandardStar *) *si;
      double RaDeg = pstar->Ra();
      double DecDeg= pstar->Dec();
      if (RaDeg > RaMin && RaDeg < RaMax && DecDeg > DecMin && DecDeg < DecMax) 
	{
	  cout << " Name : " << pstar->Name() << " Landolt Magnitude : " << pstar->Magnitude(couleur) << endl;
	  ++si;
	}
      else
	{
	  si = stdstarlist->erase(si);
	}
    }

  stdstarlist->ApplyTransfo(Pix2RaDec.invert());
  return stdstarlist;
}

int GetStandardZeroPoint(StandardStarList *standardList, SEStarList &sestarlist, const FitsHeader &header, double *zp, double *err_zp, int &nzero)
{
  GtransfoIdentity ident;
  StandardColor couleur = GetColor(header);
  double exposureTime = header.KeyVal("TOADEXPO");

  StarMatchList *MatchingList = ListMatchCollect(*Standard2Base(standardList), *SE2Base(&sestarlist), &ident, 20);

   int n = 0;
   cout << "entering GetStandardZeroPoint with " << MatchingList->size() << " matches" << endl;
   for(StarMatchIterator smi = MatchingList->begin(); smi != MatchingList->end() ; smi++)
     {
       StandardStar * stdstar = smi->s1;
       SEStar * catstar = smi->s2;
       //get rid of saturated and negative flux
       double saturation = header.KeyVal("SATURLEV");
       stdstar->fluxpersec = catstar->flux/exposureTime;
       //	   stdstar->efluxpersec = catstar->EFlux()/exposureTime;
       stdstar->efluxpersec = sqrt(catstar->flux)/exposureTime;
       stdstar->airmass = header.KeyVal("AIRMASS");
       stdstar->color = couleur;
       stdstar->flag = catstar->Flag();
       zp[n] = stdstar->Magnitude() + 2.5*log10(stdstar->fluxpersec) + (AirmassTerm(header)*stdstar->airmass) - (ColorTerm(header)*stdstar->DeltaColor());
       err_zp[n] = 1.0857*(stdstar->efluxpersec/stdstar->fluxpersec);
       cout << " Name : " << stdstar->Name() << " Flux/TOADEXPO : " << stdstar->fluxpersec << " +/ "<< stdstar->efluxpersec << " Landolt Magnitude : " << stdstar->Magnitude(couleur) <<" Zero point : "<< zp[n] << " +/- " << err_zp[n] << " x = " << catstar->x << " y = " << catstar->y << " flag =" << catstar->Flag() <<endl;
       if ( (catstar->flux > 0)
	    //	 && (catstar->flux +catstar->Fond() < saturation) )
	    //	    && (!(catstar->IsSaturated())))
	    && (catstar->Flag()==0))
	    //	    && catstar->IsOK(saturation) )
	 {
	   n++;
	 }
       zp[n]=0.0;
       err_zp[n]=0.0;
     }
   nzero = n;

  for (StandardStarIterator si = standardList->begin(); si != standardList->end(); )
    {
      StandardStar * pstar = (StandardStar *) *si;
      if (pstar->fluxpersec > 0.0)
	{
	  ++si;
	}
      else
	{
	  si = standardList->erase(si);
	}
    }

   return 1;
}

double NormFactor(StandardStar *pstar)
{
  double err_flux = pstar->efluxpersec;
  double err_flat = 0.01*pstar->fluxpersec;
  //double err_flat = 0.0 ;
  double norm = 1.0857*1.0857*(err_flux*err_flux+err_flat*err_flat)/(pstar->fluxpersec*pstar->fluxpersec) ;
  //  double norm = 1.0857*(pstar->efluxpersec)/pstar->fluxpersec ;
  if (norm > 0.0) norm = 1/norm;
  return norm;
}

double FitZeroPoint(StandardStarList &stdstarlist,double &a1, double &a2, double mean, double sig, double &err_zp, double &err_a1, double &err_a2, string outfilename)
{

  cout << "========================================"<<endl;
  cout << "   Fit of Zp, airmass and color term    "<< endl;
  cout << "(check that you have differents airmass)"<<endl;
  cout << "========================================"<<endl;

  double zp = 0.0;

  double A[3][3], B[3];
  for (int i=0;i<3;i++){B[i]=0.0;for(int j=0;j<3;j++) A[i][j]=0.0;}
  // fill matrices

  for (StandardStarIterator si = stdstarlist.begin(); si != stdstarlist.end(); si++)
    {
      StandardStar * pstar = (StandardStar *) *si;
      double norm = NormFactor(pstar);
      double Dm = pstar->Magnitude() + 2.5*log10(pstar->fluxpersec);
      if (fabs(Dm-mean)<3*sig)
	{
	  double couleur = pstar->DeltaColor();
	  double airmass = pstar->airmass;
	  
	  cout <<"Name : "<<pstar->Name()<<" Dm = "<<Dm<<" couleur = "<<couleur<<" airmass = "<<airmass<<" norm = "<<norm<<endl;

	  B[0] += Dm*norm;
	  B[1] += Dm*couleur*norm;
	  B[2] += Dm*airmass*norm;
	  
	  A[0][0] += norm;
	  A[0][1] += couleur*norm;
	  A[0][2] += airmass*norm;

	  A[1][1] += couleur*couleur*norm;
	  A[1][2] += couleur*airmass*norm;
	  A[2][2] += airmass*airmass*norm;
	}
    }

  A[1][0] = A[0][1];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];

  int ierr;
  int n = 3;
  ierr = SymMatInv(&A[0][0],n);
  if (ierr == 0) cout<<"WARNING !!! matrix A not inversed !!!"<< endl;

  zp = A[0][0]*B[0] + A[0][1]*B[1] + A[0][2]*B[2];
  a1 = A[1][0]*B[0] + A[1][1]*B[1] + A[1][2]*B[2];
  a2 = A[2][0]*B[0] + A[2][1]*B[1] + A[2][2]*B[2];

  err_zp = sqrt(A[0][0]);
  err_a1 = sqrt(A[1][1]);
  err_a2 = sqrt(A[2][2]);

  cout << endl;
  cout <<"Matrice de Variance-Covariance : "<<endl;
  cout <<"================================ "<<endl;
  cout <<" "<< A[0][0] <<" "<< A[0][1] <<" "<< A[0][2] << endl;
  cout <<" "<< A[1][0] <<" "<< A[1][1] <<" "<< A[1][2] << endl;
  cout <<" "<< A[2][0] <<" "<< A[2][1] <<" "<< A[2][2] << endl;

  for (int i=0; i<3; ++i)
    {
      cout << " ";
      for (int j=0; j<=i; ++j)
	cout << A[i][j]/sqrt(A[i][i]*A[j][j]) << " ";
      cout <<endl;
    }
      
  float chi2 = 0.0;
  float dof = 0;
  for (StandardStarIterator si = stdstarlist.begin(); si != stdstarlist.end(); si++)
    {
      StandardStar * pstar = (StandardStar *) *si;
      double norm = NormFactor(pstar);
      double Dm = pstar->Magnitude() + 2.5*log10(pstar->fluxpersec);
      if (fabs(Dm-mean)<3*sig)
	{
	  double couleur = pstar->DeltaColor();
	  double airmass = pstar->airmass;
	  double delta = Dm - (zp + a1*couleur + a2*airmass);
	  chi2 += delta*delta*norm;
	  dof += 1.0;
	}
    }

  cout << endl;
  cout << "Chi2/ndof(="<< (int)(dof -3.0)<<") = " << chi2/(dof-3.0) << endl;

  if (outfilename.length()!=0)
    {
      ofstream out(outfilename.c_str());
      out <<"#M = -2.5log10(flux/time) + zp + a1*Color + a2*Airmass"<<endl;
      out <<"#zp dzp a1 da1 a2 da2 zpa1 zpa2 a1a2 chi2/ndof"<<endl;

      out << zp <<" ";
      out << err_zp <<" ";
      out << a1 <<" ";
      out << err_a1 <<" ";
      out << a2 <<" ";
      out << err_a2 <<" ";
      out << A[0][1]/sqrt(A[0][0]*A[1][1]) <<" ";
      out << A[0][2]/sqrt(A[0][0]*A[2][2]) <<" ";
      out << A[1][2]/sqrt(A[1][1]*A[2][2]) <<" ";
      out << chi2/(dof-3.0) << endl;
    }
  return zp;
}



double ReadTheoZeroPoint(const FitsHeader &header)
{
  double expo = header.KeyVal("TOADEXPO");
  double airmass = header.KeyVal("AIRMASS");
  double a = AirmassTerm(header);

  if (header.HasKey("TOADPZPT"))
    {
      double zerop = header.KeyVal("TOADPZPT");
      zerop -= 2.5*log10(expo);
      zerop += a*airmass;
      return zerop;
    }
  cout <<" NO ZERO POINT FOUND IN DATABASE !!!!!!!!!!!!! "<<endl;
  return 0.0;
}



/************************** FINDEFINITION StandardStar ************************/

#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

template class StarList<StandardStar>;

#ifdef USE_ROOT
template class StarListWithRoot<StandardStar>;
ClassImpT(StarListWithRoot,StandardStar);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT

LINKDEF_CONTENT : #pragma link C++ class StandardStar;
LINKDEF_CONTENT : #pragma link C++ class list<StandardStar*>;
LINKDEF_CONTENT : #pragma link C++ class StarList<StandardStar>-;
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<StandardStar>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<StandardStar>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<StandardStar>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef StandardStarIterator;
*/
#include "root_dict/standardstardict.cc"
#endif /* USE_ROOT */
