// This may look like C code, but it is really -*- C++ -*-
#ifndef ASTROUTILS__H
#define ASTROUTILS__H

class FitsHeader;
#include <string> /* for prototyping */


/*! \file */

bool IsLeapYear (const int year);
long JulianDay(const int day, const int month, const int year);
double JulianDay(const int day, const int month, const int year, 
		 const int hour, const int min, const double sec);
//! computes the julian date given a header
double JulianDay(const FitsHeader& Header);
//! computes the reduced julian date given a header
double RedJulianDay(const FitsHeader& Header);

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
#endif //ASTROUTILS__H
