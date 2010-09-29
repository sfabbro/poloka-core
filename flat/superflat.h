#ifndef SUPERFLAT__H
#define SUPERFLAT__H

#include "image.h"
#include <vector>

class FitsSet;

double GainMultiply(FitsImage &Current);

Image *MakeSuperFlat(const FitsSet &FitsFileSet, const Image *Bias, const Image *SkyFLat);

Image *DeadPixelImageAndFlatSmoothing(Image &Flat, double FlatMin, double FlatMax);

double ImageMaxCut   (const Image &Current, const Image &Flat);

Image *MakeRawAverage(const FitsSet &FitsFileSet);
Image *MakeRawAverageAndSigma(const FitsSet &FitsFileSet, Image *Sigma, const int normalize = 0, const int normTime = 0);

Image *MakeRawMedian(const FitsSet &FitsFileSet, const int normTime = 0);
Image *MakeRawMaskedMedian(const FitsSet &FitsFileSet, const FitsSet &MaskFitsFileSet, const int normTime = 0);

void MakeRawAverageAndSigma(const FitsSet &FitsFileSet, const Image &Bias, const Image &SkyFlat, const vector<double> &norm, Image &Mean, Image &Sigma);


int FlatFieldImage(const string &InFileName, const string &FlatName, 
		   const string &BiasName, const string &FringeName, 
		   const string &OutFileName);

void  RawHighFreqFilter(Image &Img, const int GridSize);

void SubNormFringes(const string &InFileName, const Image &Fringe,
		    const bool norm);

void ImageAlreadyFlatFielded(const string &InFileName, const Image &Flat,
			     const string &OutFileName);
double ImageNormalisation(FitsImage &Current);

double BiasCorrect_and_trim(FitsImage &Current, const Image *Bias=NULL);

double SurfaceFit(const Image &Im, const int BackDegree);

#ifdef STORAGE
void FilterFringes(Image &FringePattern, const int GridSize, const int NInc);
  //Creates a inclination map, by sweeping a GridSize X GridSize/2 window on each pixels, catching the 
  // optimal inclination. To search for optimality we do for each angle:
  // - rotate and interpolate
  // - get the median over each column
  // - get the rms of medians
  // - picks up the angle where rms in minimal
  // InclinationMap(i,j) = bestangle
  //Then a single one pixel wide window (length is typically 10pixels) 
  //is properly orientated on the original given fringe frame. Finally the median of the 
  //intensity distribution over this window is adopted as the final intensity at the 
  //given position. (Ostrov PASP 109,p.338)



Image* MakeFringePattern(const Image &FringedImage,const Image &BlankImage, const double &Nsig, const bool &BackFlag, const bool &FilterFlag);
  //Make a fringe map by simply dividing an image (usually a superflat) by a another flat.
  // This other flat should be a blank, objectless, fringeless image, typically a dome or twillight flat.
  // It can also froce the min and max, remove a background and filter the fringes.
  // Note that the 2 input images should be normalized the same.

void RemoveFringes(Image &FringedImage, const Image &FringeMap,const bool &NormFlag);
  // Remove fringe pattern. Calculate new sigma and sky value by taking the normalisation factor 
  // of the fringe map as the one representing the best sky sigma.
  // This is done cause fringes are at the specific OH wavelength whereas sky noise is spread out
  // Assume also FringeMap is normalized to 1.
#endif

#endif /* SUPERFLAT__H */
