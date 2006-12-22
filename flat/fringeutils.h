/* 
 * $Source: /cvs/snovae/toads/poloka/flat/fringeutils.h,v $
 * $Revision: 1.2 $
 * $Author: guy $
 * $Date: 2006/12/22 13:35:40 $
 * $Name:  $
 */

#ifndef FRINGEUTILS__H
#define FRINGEUTILS__H

#include <vector>
#include <string>

class Image;
class Mat;
class vect;

//! Utility class for fringe analysis
/*! Used by fringefinder and defringe2
 */
class FringeUtils {

 public:
  //! Default Constructor
  FringeUtils(){};
  //! Default Destructor
  ~FringeUtils(){};
  
  //! Checks whether images have same size or not
  static bool IsSameSize(const Image &Img1,const Image &Img2);
  
  //! Computes Scalar Product
  /*! \return 1/p*sum(Img1(i)*Img2(i)) for( Img1(i)>nsigma*sigma1i && Img2(i)>nsigma*sigma2i )
   * if nsigma<0 (default), no cut is applied 
   */
  static double ScalarProduct(Image &Img1, Image &Img2,float nsigma, const Image* deadimage=0);
  //! Smooth an image
  /*! Convolution of Image Img with a smoothing filter 3x3
   *  0.025, 0.100, 0.025
   *  0.100, 0.500, 0.100
   *  0.025, 0.100, 0.025
   */
  static void SmoothFilter(Image &Img);
  
  //! Returns scalar product matrix
  /** ScalarProductMatrix is a symmetric matrix n*n 
   * computed with n images in a FitsFileSet.
   *  Each element is the ScalarProduct(Image_i,Image_j)  
   * \param  FitsFileSet set of images
   * \return double array[n*n] containing the matrix
   */
  static Mat ScalarProductMatrix(vector<string> &filelist, const Image* deadimage=0);
  
  //! Clear an image 
  /*! Set pixels values to 0 if fabs(value-mean)>NSigma*rms 
   */
  static void ClearImage(Image &Img, const double NSigma);
  
  //! Normalize an image -> (mean=0, rms=1) 
  /*! Only considers pixels for fabs(value-mean)>nsigma_cut*rms (default nsigma_cut=3)
   */
  static void Normalize(Image &img,const float nsigma_cut=3);

  //! Greatest Common Divider
  static int GreatestCommonDivider(int a, int b);
  
  //!Clipped average and rms
  static void ClippedAverageRms(double *pixelValues, const int count,
				double &average, double &rms);

  //! Returns the CVS version of this piece of code
  static string GetVersion();

  //! Remove fringes in an image 
  /*! \return 0 is everything is ok
   * \param image must be in RW mode
   * \param fringemapname name of the FITS file containing fringe patterns (one or several)
   * \param nvec max number of vectors to use (default=0 means all of them)
   * \param nsig cut on the image for the ScalarProduct (default=3), this cut is set to 5 for the fringes
   * \param substractbg if true, substract background using imageback before removing fringes 
   *        (and put it back afterwards) (default is false)
   * \param verbose if true prints a lot of stuff (default is true)
   */
  static int RemoveFringes(FitsImage &image, const string &fringefilename, int nvec=0, float nsig = 3, bool substractbg=false, bool verbose=true);
  
  //! Checks whether this image is defringed or not
  static bool IsDefringed(const FitsImage &image);
  
  //! Checks whether this image is defringed or not with the new method
  static bool IsDefringedWithNewMethod(const FitsImage &image);
  
   //! Checks whether this image is a fringe pattern
  static bool IsANewFringePattern(const string &filename);
  
};


#endif
