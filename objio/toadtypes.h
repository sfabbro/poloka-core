// -*- C++ -*-
// 
// $Id: toadtypes.h,v 1.1 2004/02/20 10:48:43 nrl Exp $
// 
// file toadtypes.h 
// 
/*!
  \file toadtypes.h
  \brief portables types used within the toads framework
  
  Last modified: $Date: 2004/02/20 10:48:43 $
  By:            $Author: nrl $
  
 */
#ifndef TOADTYPES_H
#define TOADTYPES_H

#include "config.h"

// 
// 1 byte int: int1 & uint1
// 
#if defined(SIZEOF_CHAR) && (SIZEOF_CHAR==1)
 typedef signed char int1;
 typedef unsigned char uint1;
#else
#error "Please define a one byte integer type"
#endif


// 
// 2 byte int: int2 & uint2
// 
#if defined(SIZEOF_SHORT) && (SIZEOF_SHORT==2)
 typedef signed short int2;
 typedef unsigned short uint2;
#else
#error "Please define a two byte integer type"
#endif


// 
// 4 byte int: int4 & uint4
// 
#if defined(SIZEOF_INT) && (SIZEOF_INT==4)
 typedef signed int int4;
 typedef unsigned int uint4;
#else
#error "Please define a four byte integer type"
#endif


// 
// 8 byte int: int8 & uint8
// 
#if defined(SIZEOF_INT) && (SIZEOF_INT==8)
 typedef signed long int8;
 typedef unsigned long uint8;
#elif defined(SIZEOF_LONG_LONG) && (SIZEOF_LONG_LONG==8)
 typedef signed long long int8;
 typedef unsigned long long uint8;
#else
#error "Please define an eight byte integer type"
#endif


// 
// 4 byte float: r4 
// 
#if defined(SIZEOF_FLOAT) && (SIZEOF_FLOAT==4)
 typedef float float4;
#else
#error "Please define a four byte float type"
#endif


// 
// 8 byte float: r8
// 
#if defined(SIZEOF_DOUBLE) && (SIZEOF_DOUBLE==8)
 typedef double float8;
#else
#error "Please define a eight byte float type"
#endif

#endif

