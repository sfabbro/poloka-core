#ifndef IMAGEUTILS__H
#define IMAGEUTILS__H



class Image;
class Gtransfo;

//! No resampling, only shift an image within a larger frame
Image IntegerShiftImage(const Image& inputimage, const Gtransfo &g,
			int nx, int ny, float DefaultVal);

//! Image resampling. Can handle Variance maps
Image GtransfoImage(const Image& inputimage, const Gtransfo & g, int nx, int ny, 
		    float DefaultVal, const int interpLevel=3, 
		    const bool IsVarianceMap = false);

//! still used as a check, as saturation value are still sometimes misset.
double ComputeSaturation(const Image& image); 

double ComputeDrift(const Image& Im);

/*!  "Convolve" a mask image consisting of "patches" filled
with a uniform value separated by 0's. we enlarge the patches
with the same value, and put -1 in case of conflicting values
Ny deals with the cases when we don't want to convolve up to the whole image height*/
bool ConvolveSegMask(const Image &In, Image &Out, const int ExtraSize, int Ny=-1);

#endif /* IMAGEUTILS__H */
