// This may look like C code, but it is really -*- C++ -*-
#ifndef ASTROUTILS__H
#define ASTROUTILS__H

class FitsHeader;
#include <string> /* for prototyping */


/*! \file */

//! return true if year is a leap year
bool IsLeapYear (const int year);

//! return the integer julian date given day, month and year
double JulianDay(const int day, const int month, const int year);

//! return the integer julian date given day, month, year, hour, minutes and seconds
double JulianDay(const int day, const int month, const int year, 
		 const int hour, const int min, const double sec);

//! compute the julian date given a FITS header
double JulianDay(const FitsHeader& Header);

//! compute the julian date given a string "DD/MM/YYYY"
double JulianDay(const string& DateString);

//! compute the modified julian date given a header (jd-2400000.5)
double ModJulianDay(const FitsHeader& Header);

//! 
double ModifiedModifiedJulianDay(const FitsHeader& Header);

//! compute date stuff given a julian date
void DateFromJulianDay(const double &jd, int& day, int& month, int& year, 
		       int& hour, int& min, double& sec);

//! convert modified julian date to a FITS date (format yyyy-mm-ddThh:mm:ss.sss)
string FitsDateFromModJulDay(const double& mjd);

//! convert a FITS date to modified julian date
double ModJulDayFromFitsDate(const char *date);

//! test routine to compare with JulianDay
double JulDate(const int day, const int month, const int year, 
	       const int hour, const int min, const double sec);

//! compute "DD/MM/YYYY" given a julian date
string DateFromJulianDay(const double &jd);

//! 
double UtStringToDeci(const string UtString);
//! -
double RaStringToDeg(const string RaString);
//! -
double DecStringToDeg(const string DecString);
//! -
string UtDeciToString(double UtDeci);
//!
string RaDegToString(double RaDeg);
//! -
string DecDegToString(double DecDeg);

//! -
void GetRaDecInDeg(FitsHeader const &Header, double &RaInDeg, double &DecInDeg);

//! returns alpha delta in 2000 whatever the equinox used in header.
void RaDec2000(FitsHeader const &Header, double &RaInDeg, double &DecInDeg);

//! Deep field name (D1, D2, D3, D4). specific to CFHTLS.
bool IdentifyDeepField(const FitsHeader &Head, string &Field);
#endif //ASTROUTILS__H
