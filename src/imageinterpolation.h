#ifndef IMAGEINTERPOLATION__H
#define IMAGEINTERPOLATION__H

class Image;
class Frame;

Pixel Interpolate(const Image& inputimage, const double x, const double y, const int level=3, const bool IsVarianceMap=false);

#endif /* IMAGEUTILS__H */
