/* 
 * $Source: /cvs/snovae/toads/poloka/flat/fitsimagearray.h,v $
 * $Revision: 1.2 $
 * $Author: guy $
 * $Date: 2006/12/22 13:35:40 $
 * $Name:  $
 */

#ifndef FITSIMAGEARRAY__H
#define FITSIMAGEARRAY__H

#include <vector>
#include "fitsimage.h"

//! Array of Image in a FITS file 
/*! Inherites from FitsImage
 */ 
class FitsImageArray : public FitsImage {
  
  public :
    
    //!  Constructor for an existing fits image
     /*! \li Opens in (ReadOnly (RO) mode by default, use RW to modify or create)
     * \li In RW mode, use this creator for an image already containing several images, i.e. with EXTEND=true
     *  (if EXTEND=false, IsValid()=false)
     */
    FitsImageArray(const string &FileName, const FitsFileMode Mode = RO);
  
  //! Constructor for a new FitsImageArray with a minimal header from an existing FitsHeader. 
   /*! This header is used as a reference for image size, filter, chip.
   * A criterion for accepting images can be set with the function SetCriterion 
   * using the values of the enum FitsImageArray::Criterion. 
   * For instance SetCriterion(FitsImageArray::FILTER|FitsImageArray::CHIP) <=> ask for the same filter and chip.
   *
   * \param ImageInFirstHDU, if true (default) the first image that is added is saved in the first HDU,
   * otherwise it is saved in the second one.
   */
  FitsImageArray(const string &FileName, const FitsHeader& header, bool ImageInFirstHDU=true); 
  
  //! Constructor for a new FitsImageArray with a minimal header from an existing Image. 
   /*! This Image is used as a reference for image size (see previous constructor).
    * The user can define other keys for accepting other images.
    *
    * \param ImageInFirstHDU, if true (default) the first image that is added is saved in the first HDU,
    * otherwise it is saved in the second one. 
    */  
  FitsImageArray(const string &FileName, const Image& image, bool ImageInFirstHDU=true); 
  
  //! Standard destructor
  ~FitsImageArray();
  
  //! Status flag for all functions
  enum Status {OK = 0, OUTOFBOUNDS = 1, FAILURE = 2, ERRORWRITE = 3, EMPTY = 4};
  
  //! Criterion fo image selection
  enum Criterion {CHIP = 1, FILTER = 2, SIZE = 4, EXTENDED = 8};

  //! Write all images to the FITS file
  Status Write(bool force_bscale = false);

  //! Write all images to separated FITS files (same as split_fits)
  Status SplitAndWrite(const string &directory="",int HDUMax=0);
  
  //! Add a FitsImage after the latest HDU
   /*! \param filename, name of the fits file to append
   * Uses cfitsio fits_copy_hdu, i.e. toads is not used a lot
   * except for checking the compatibility of FitsHeader
   */
  Status Append(const string &filename);

  //! Add an Image after the latest HDU
   /*! \param Image that is appended.
   * \li A default FitsHeader is written, namely with NAXIS*, 
   * BSCALE, BZERO, BITPIX=16, and WRITEDAT
   * \li The user may add information to the header after
   * the call to this function 
   */
  Status Append(const Image &image);
  
  //! Point to the image at index HDU ( 1 <= HDU <= GetNImages())
  Status At(int HDU);
  
  //! Point to the next image
  Status Next();
  
  //! Set the criterion for merging of FitsImage
  /*! \param value : example FitsImageArray::CHIP|FitsImageArray::FILTER
   */
  void SetCriterion(int value) {fSelectedCriterion=value;};

  //! Check whether this FitsImageArray is valid
  bool IsValid();
  
  //! Returns the version of this piece of code
#ifndef SWIG
  string GetVersion();
#endif

  
  //! Check whether header is consistent with fMainHeader
  /*! \li Used by FitsImageArray::Append
   * \li Uses fSelectedCriterion (see SetCriterion)
   */
  bool CheckHeader(FitsHeader &header);
  
  private :
    //! Provides a name for the splitted file (used by FitsImageArray::SplitAndWrite)
    void SplitName(string &Dir, string &Base, string &Type);
  
  //! Read image at current HDU
  Status ReadImage();

  //! Fill the header for an image
  /*! \li Used in Append(const Image &image)
   * \li If EXTNAME does not exist, creates default im%03d (can be overwritten by user after Append)
   * \li Overwrites EXTVER
   * \li Adds date of creation
   * \li Computes BSCALE & BZERO (copied from fitsimage)
   */
  void FillHeader(const Image& image,bool force_bscale = false);
  
  //! used by creators
  Status Init(bool already_extended);
  
  int fCurrentHDU;
  int fSelectedCriterion;
  bool fImageInFirstHDU;
  FitsHeader *fMainHeader;
  bool fValid;
};


#endif

