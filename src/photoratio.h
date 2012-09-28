// This may look like C code, but it is really -*- C++ -*-
#ifndef PHOTORATIO__H
#define PHOTORATIO__H

#include "basestar.h"
#include "reducedimage.h"
#include "gtransfo.h"


enum PhotoRatioMethod { ZeroPointDiff = 1,
			TotalLeastSquares = 2,
			MedianRatio = 3,
			AverageRatio = 4,
			PUnSet = 3};

string NamePhotoRatioMethod(const PhotoRatioMethod RatioMethod);

// photometric ratio = image / reference

//! quick and dirty photometric ratio (weighted robust mean of two catalogs)
double AveragePhotoRatio(const StarMatchList& MatchList, double& Error);

//! wrapper of the above with ReducedImage's 
double AveragePhotoRatio(const ReducedImage &Im,
			 const ReducedImage &Ref, 
			 double &Error,
			 const Gtransfo* Im2Ref=0);

//! quick and robust photometric ratio with StarMatchList (s1->flux/s2->flux)
double MedianPhotoRatio(const StarMatchList& MatchList, double& Error);

//! wrapper of the above with ReducedImage's 
double MedianPhotoRatio(const ReducedImage &Im,
			const ReducedImage &Ref, 
			double &Error,
			const Gtransfo* Im2Ref=0);

//! compute photometric ratio such that f1=R*f2 on average with total least squares
double TLSPhotoRatio(const StarMatchList& MatchList, double& Error, const double& NSigChi2Cut=5);

//! wrapper of the above with ReducedImage's 
double TLSPhotoRatio(const ReducedImage &Im,
		     const ReducedImage &Ref, 
		     double &Error,
		     const Gtransfo* Im2Ref=0);


//! simple difference of zero points scaled back to flux
double ZpPhotoRatio(const double& Zp, const double& SigZp, 
		    const double& ZpRef, const double& SigZpRef,
		    double& Error);

//! wrapper of the above with ReducedImage's and scanning zero points key names in header
double ZpPhotoRatio(const ReducedImage& Im, const ReducedImage& Ref, double& Error);


double PhotoRatio(const ReducedImage& Im, const ReducedImage& Ref,
		  double& Error,
		  const Gtransfo* Im2Ref=0,
		  const PhotoRatioMethod Method=TotalLeastSquares);

#endif // PHOTORATIO__H
