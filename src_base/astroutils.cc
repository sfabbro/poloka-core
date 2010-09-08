#include <cmath>
#include <cstdio>

#include "fitsimage.h"
#include "astroutils.h"

bool IsLeapYear (const int year)
{
   return year % 400 == 0 || (year % 4 == 0 && year % 100 != 0);
}

// based on IDL code juldate.pro in astrolib
// changed in april 06, because it was heavily bugged. now return a double
// (which is a half integer).
// Checked that it matches the idl output. Maybe the idl code was fixed
// since we first "borrowed" it.
// 07/2006: still buggy. seems to return a modjulianday.
// This routine should never have been used. JulDate was the correct one (borrowed from a C routine)

double JulianDay(const int day, const int month, const int year) 
{
  long jy, jm, ja;
  if ( month > 2 )
    {
      jy = year;
      jm = month;
    }
  else
    {
      jy = year - 1;
      jm = month + 12;
    }

  double a = long(jy/100);
  double ry = jy;
  
  double jul = floor(ry*0.25) + 365.*(ry-1860.) + int(30.6001*(jm+1)) + 
    day -105.5;
  if (jul > -100830.5) jul += 2 -a + floor(a/4);

  return jul; 
}

double JulianDay(const int day, const int month, const int year, 
		 const int hour, const int min, const double second)
{
  //return JulianDay(day,month,year) + ( hour + min/60.0 + second/3600.0) / 24.0;
  return JulDate(day, month, year, hour, min, second);
}

double ModifiedModifiedJulianDay(const FitsHeader& Header) {
  //return JulianDay(Header)-JulianDay(01,01,2003;)
  return JulianDay(Header) - JulianDay(01, 01, 2003, 0, 0, 0);
}


double JulianDay(const FitsHeader& Header)
{
  if (Header.HasKey("MJD-OBS")) {
    double mjd = Header.KeyVal("MJD-OBS");
    return mjd + 2400000.5; 
  }

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

double JulianDay(const string& DateString)
{
  // compute julian day at noon UT
  int day,year,month;
  if (sscanf(DateString.c_str(),"%d/%d/%d",&day, &month, &year) != 3) 
    {
      cerr << " JulianDay(" << DateString << ") : Error : can't compute Julian date from string " << endl; 
      return 0.;
    }

  return JulianDay(day,month,year);
}

// see def e.g. in http://tycho.usno.navy.mil/mjd.html
double ModJulianDay(const FitsHeader& Header)
{
  return JulianDay(Header) - 2400000.5;
}

string FitsDateFromModJulDay(const double& mjd) 
{
  int day, month, year, hour, min;
  double sec;
  DateFromJulianDay(mjd+2400000.5, day, month, year, hour, min, sec);
  char datechr[25];
  sprintf(datechr,"%04d-%02d-%02dT%02d:%02d:%06.3f", year, month, day, hour, min, sec);
  return string(datechr);
}

double ModJulDayFromFitsDate(const char *date)
{  
  int day, month, year, hour, min;
  double sec;
  if (sscanf(date, "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &min, &sec) != 6) {
    cerr << " Can't compute modified Julian date from date " << date << endl;
    return 0.;
  }
  return JulDate(day, month, year, hour, min, sec) - 2400000.5;
}

// directly converted from caldate in jday.sf.net
void DateFromJulianDay(const double &jd, int& day, int& month, int& year, 
		       int& hour, int& min, double& sec)
{
  long ljd = (long) (jd + 0.5);
  double frac = jd + 0.5 - (double) ljd + 1.0e-10;
  long ka = (long) ljd;
  long ialp;
  if ( ljd >= 2299161L )
    {
      ialp = long( (double(ljd) - 1867216.25 ) / ( 36524.25 ));
      ka = ljd + 1L + ialp - ( ialp >> 2 );
    }
  long kb = ka + 1524L;
  long kc = long((double(kb) - 122.1 ) / 365.25);
  long kd = long( double(kc) * 365.25);
  long ke = long(double( kb - kd ) / 30.6001);

  day = kb - kd - ((long) ( (double) ke * 30.6001 ));

  if ( ke > 13L ) month = ke - 13L;
  else	month = ke - 1L;

  if ( (month == 2) && (day > 28) )    day = 29;
  if ( (month == 2) && (day == 29) && (ke == 3L) ) year = kc - 4716L;
  else if ( month > 2 ) year = kc - 4716L;
  else year = kc - 4715L;

  double dhour= frac*24.0;
  hour = long(dhour);

  double dmin = ( dhour - double(hour) ) * 60.0;
  min = long(dmin);
  sec  = ( dmin -  double(min) ) * 60.0;
}

string DateFromJulianDay(const double &jd)
{
  double sec;
  int day,month,year,hour,min;
  DateFromJulianDay(jd,day,month,year,hour,min,sec);    
  char datechr[11];
  sprintf(datechr,"%02d/%02d/%04d",day,month,year);
  return string(datechr);
}


double JulDate(const int day, const int month, const int year, 
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



// CFHTLS deep fields
// not really checked the coordinates.... P.A. (2005)
static const double DeepRa[4] = {36.5, 150.12, 214.87, 333.88};
static const double DeepDec[4] = {-4.50, 2.20, 52.68, -17.73};

bool IdentifyDeepField(const FitsHeader &Head, string &Field)
{
  Field = "";
  string obj = string(Head.KeyVal("TOADOBJE"));
  if (obj == "D1" || obj == "D2" || obj == "D3" || obj == "D4")
    {
      Field = obj;
      return true;
    }
  string ras = string(Head.KeyVal("RA"));
  double ra = RaStringToDeg(ras);
  string decs = string(Head.KeyVal("DEC"));
  double dec = DecStringToDeg(decs);
  for (int i=0; i<4; ++i)
    if (fabs(ra-DeepRa[i])<0.2 && fabs(dec-DeepDec[i])<0.2)
      {
	char field[8];
	sprintf(field,"D%d",i+1);
	Field = field;
	return true;
      }
  return false;
}

