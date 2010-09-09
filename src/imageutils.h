#ifndef IMAGEUTILS__H
#define IMAGEUTILS__H



class Image;
class Gtransfo;
#include "frame.h" // for whichTransformed

//! No resampling, only shift an image within a larger frame
Image IntegerShiftImage(const Image& inputimage, const Gtransfo &g,
			int nx, int ny, float DefaultVal);

//! Image resampling. Can handle Variance maps
Image GtransfoImage(const Image& inputimage, const Gtransfo & g, int nx, int ny, 
		    float DefaultVal, const int interpLevel=3, 
		    const bool IsVarianceMap = false);


//! assumes that Transfo is a shift or involves a 'simple rotation'
Frame ApplyTransfo(const Frame& inputframe, const Gtransfo &T, const WhichTransformed W = SmallFrame);


//! still used as a check, as saturation value are still sometimes misset.
double ComputeSaturation(const Image& image); 

/*!  "Convolve" a mask image consisting of "patches" filled
with a uniform value separated by 0's. we enlarge the patches
with the same value, and put -1 in case of conflicting values */
bool ConvolveSegMask(const Image &In, Image &Out, const int ExtraSize);

#endif /* IMAGEUTILS__H */
