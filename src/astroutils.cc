#include <cmath>
#include <cstdio>

#include "fitsimage.h"
#include "astroutils.h"


#ifndef M_PI
#define M_PI    3.14159265358979323846  
#endif

bool IsLeapYear (const int year)
{
   return year % 400 == 0 || (year % 4 == 0 && year % 100 != 0);
}

long JulianDay(const int day, const int month, const int year) 
{
  long greg = 15 + 31 * (10 + 12 * 1582);
  long jy, jm, ja, jul;
  if ( month > 2 )
    {
      jy = year;
      jm = month + 1;
    }
  else
    {
      jy = year - 1;
      jm = month + 13;
    }
  
  jul = long((365.25 * jy) + (30.6001 * jm) + day + 1720995);
  if ( (day + 31 * (month + 12 * year)) >= greg )
    {
      ja = (long int)(0.01 * jy);
      jul = jul + long(2 - ja + (0.25 * ja));
    }  

  return jul;
}

double JulianDay(const int day, const int month, const int year, 
		 const int hour, const int min, const double second)
{
  return JulianDay(day,month,year) + ( hour + min/60.0 + second/3600.0) / 24.0;
}

double JulianDay(const FitsHeader& Header)
{
  string toaddate = Header.KeyVal("TOADDATE");

  // compute julian day at noon UT
  int day,year,month,hour,minute;
  double second;
  if (sscanf(toaddate.c_str(),"%d/%d/%d",&day, &month, &year) != 3) 
    {
      cerr << " Can't compute Julian date from date " << toaddate << endl; 
      return 0.;
    }

  // add UT if any
  string toadutim = Header.KeyVal("TOADUTIM");
  if (sscanf(toadutim.c_str(),"%d:%d:%lf",&hour,&minute,&second) != 3) 
    {
      return double(JulianDay(day,month,year));
    }

  return JulianDay(day,month,year,hour,minute,second);
}

double RedJulianDay(const FitsHeader& Header)
{
  return JulianDay(Header) - 2400000.;
}


double UtStringToDeci(const string UtString)
{
int hours, minutes; double seconds;
if (sscanf(UtString.c_str(),"%d:%d:%lf", &hours, &minutes, &seconds) != 3)
  {
    cerr << " cannot decode a universal time in " << UtString << endl;
    return 100.0; /* invalid value ! */
  }
return ( double(hours) + double(minutes)/60. + seconds/3600.);
}

string UtDeciToString(double UtDeci)
{
int hours = int(UtDeci);
UtDeci = (UtDeci-hours)*60.;
int minutes = int(UtDeci);
double seconds = (UtDeci-minutes)*60.;
char a_string[128];
sprintf(a_string,"%2d:%02d:%06.3f",hours,minutes,seconds);
return string(a_string);
}

double RaStringToDeg(const string RaString)
{
int hours, minutes; double seconds;
if (sscanf(RaString.c_str(),"%d:%d:%lf", &hours, &minutes, &seconds) == 3)
  {
    return 15.*( double(hours) + double(minutes)/60. + seconds/3600.);
  }
// try to decode a decimal number
 char *endptr;
 const char* startptr = RaString.c_str();
 double ra = strtod(startptr, &endptr);
 if (endptr != startptr) // means successful conversion by strtod
   return ra;
 cerr << " cannot decode a right ascencion in " << RaString << endl;
 return 100.0; /* invalid value ! */
}


double DecStringToDeg(const string DecString)
{
  /* separators may be either : or ' */
char dec[64];
strcpy(dec, DecString.c_str());
double minus_char = 1; /* no minus sign */
for (char *p = dec; *p; p++) 
  {
  if (*p == ':' || *p == '\'' ) *p = ' ';
  if (*p == '-') minus_char = -1.;
  }
int deg, minutes; double seconds;
if (sscanf(dec,"%d %d %lf",&deg, &minutes, &seconds) == 3)
  {
    double sign = 1;
    if (deg < 0) sign = -1.;
    else if (deg == 0) sign = minus_char;
    return sign*(fabs(double(deg))+fabs(double(minutes))/60. + fabs(double(seconds))/3600.);
  }
 char *endptr;
 const char *startptr = DecString.c_str();
 double decVal = strtod(startptr, &endptr);
 if (startptr != endptr) // successful conversion bu strtod
   return decVal;
 cerr << " cannot decode a declination in : " << DecString << endl;
 return 200;
}

string RaDegToString(double RaDeg)
{
RaDeg /=15;
int hours = int(RaDeg);
RaDeg = (RaDeg-hours)*60.;
int minutes = int(RaDeg);
double seconds = (RaDeg-minutes)*60.;
char a_string[128];
sprintf(a_string,"%2d:%02d:%06.3f",hours,minutes,seconds);
return string(a_string);
}

string DecDegToString(double DecDeg)
{
char sign = (DecDeg < 0)? '-' : ' ';
DecDeg = fabs(DecDeg);
int degrees = int(DecDeg);
DecDeg = (DecDeg-degrees)*60.;
int minutes = int(DecDeg);
double seconds = (DecDeg-minutes)*60.;
char a_string[128];
sprintf(a_string,"%c%02d:%02d:%05.2f",sign,degrees,minutes,seconds);
return string(a_string);
}

void GetRaDecInDeg(const FitsHeader &Header, double &ra, double &dec)
{
  string raString =  Header.KeyVal("TOADRASC");
  string decString =  Header.KeyVal("TOADDECL");
  ra = RaStringToDeg(raString);
  dec = DecStringToDeg(decString);
}

#include "wcscon.h"

void RaDec2000(const FitsHeader &Header, double &Ra, double &Dec)
{
  GetRaDecInDeg(Header,Ra,Dec);
  if (!Header.HasKey("TOADEQUI"))
    {
      cerr << "RaDec2000 : no equinox found in header. Cannot know if precession needed" << endl;
      return;
    }
  double equi = Header.KeyVal("TOADEQUI");
  if (equi == 2000.0)  return;
  if (equi == 0.0)
    {
      cerr << " RaDec2000 no equinox : assuming J2000 in file " 
	   << Header.FileName() << endl;
      return ;
    }
  if (equi == 1950.)
    {
      wcscon(WCS_B1950,WCS_J2000, 0.0, 0.0, &Ra, &Dec, 0.0);
      return;
    }
  cerr << " cannot precess convert safely ra,dec in " << Header.FileName() << " TOADEQUI : " << equi << endl;
}
 





#ifdef STORAGE
SkyWindow::SkyWindow(const FitsHeader& Head)
{
  double ra,dec;
  GetRaDecInDeg(Head,ra,dec);
  cosdec = cos(dec*M_PI/180.);
  epoch = Head.KeyVal("TOADEQUI");
  if (epoch != 2000) 
    {
      cerr << "WARNING: Epoch "<< epoch << " is not yet implemented" << endl;
      return SkyWindow();
    }
   int nx = Header.KeyVal("NAXIS1");
   int ny = Header.KeyVal("NAXIS2");
   Pix2Deg.apply(0,0,RaMin,DecMin);
   Pix2Deg.apply(nx-1,ny-1,RaMax,DecMax);
   RaMin /= cos(DecMin*M_PI/180.);
   RaMax /= cos(DecMax*M_PI/180.);
   
  Frame(ra + minra/cosdec, dec + mindec, ra + maxra/cosdec, dec + maxdec);
}

SkyWindow SkyWindow::operator * (const SkyWindow &Right) const
{
  SkyWindow result = *this;
  result *= Right;
  return result;
}

SkyWindow& SkyWindow::operator *= (const SkyWindow &Right)
{
  RaMin = max(RaMin,Right.RaMin);
  RaMax = min(RaMax,Right.RaMax);
  DecMin = max(DecMin,Right.DecMin);
  DecMax = min(DecMax,Right.DecMax);
  return *this;
}

double SkyWindow::Area() const
{
  return fabs((RaMax - RaMin)*(DecMax - DecMin));
}

bool SkyWindow::InSkyWindow(const double &ra,const double &dec) const
{
  return (( ra<= RaMax) && (dec<=DecMax) && (ra>=RaMin) && (dec>=DecMin));
}

void SkyWindow::dump(ostream & stream) const
{
  stream << "RA center" << RaCent << "(" << RaDegToString(RaCent) << ")  " 
	 << "DEC center" << DecCent << "(" << DecDegToString(DecCent) << ")" << endl
	 << "RA min DEC min " 
	 << RaMin << ' ' << DecMin 
	 << " RA max DEC max" 
	 << RaMax << ' ' << DecMax << endl;
}

#endif
