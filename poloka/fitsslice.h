#ifndef FITSSLICE__H
#define FITSSLICE__H

#include <string>
#include <vector>

#include <poloka/fitsimage.h>

/*! \file 
    \brief  Load successively slices of fits images in memory.

    When one needs image data for a single traversal, FitsSlices
    and FitsParallelSlices should be considered. At the moment
    both classes only provide input services. See \ref example_slices
    for an example.
*/

//! Memory saving traversal of a fits image.
/*! This class is intended to access an image in a fits file, but
only a slice of the image will reside in memory at a time.
To be used when many images have to be accessed at the same time,
in parallel (see FitsParallelSlices), or maybe to traverse a very big image. */

class FitsSlice : public FitsHeader, public Image
{
  private :
    int sliceSize, overlap;
    int ySliceStart; /* in whole image coordinates */
    int nyTotal;
    int read_pixels(const int StartRow, int &NRows, Pixel *Where);
    bool lastSlice;

  public :
 //! constructor. Overlap is the number of rows in common to successive slices   (usually 0).
   FitsSlice(const string FileName, const int SliceYSize, const int Overlap);

 //! load next slice into memory.
   int LoadNextSlice();

 //! the current slice size (differs from the constructor value at last slice).
   int SliceSize() const { return sliceSize;};

 //! the coordinate in the full image of row j in present slice.
   int ImageJ(const int j) const { return (j+ySliceStart);};

 //! are we in the last slice?
   bool LastSlice() const {return lastSlice;};


};


class FitsOutSlice : public FitsHeader, public Image
{
  private :
    int sliceSize, overlap;
    int ySliceStart; /* in whole image coordinates */
    int nyTotal;
    int bitpix;
    double bscale, bzero;
     int write_pixels(const int StartRow, int &NRows, Pixel *Where);
    bool lastSlice;

  public :
 //! Constructor for a non existing image. Overlap is the number of rows in common to successive slices   (usually 0).
   FitsOutSlice(const string FileName, int Nx, int Ny, 
		const int SliceYSize, const int Overlap);

 //! constructor. Overlap is the number of rows in common to successive slices   (usually 0).
   FitsOutSlice(const string FileName, const FitsHeader &ModelHead,
	     const int SliceYSize, const int Overlap);

 //! write current slice to disk and clears the pixel buffer.
   int WriteCurrentSlice();

 //! the current slice size (differs from the constructor value at last slice).
   int SliceSize() const { return sliceSize;};

 //! the coordinate in the full image of row j in present slice.
   int ImageJ(const int j) const { return (j+ySliceStart);};

 //! are we in the last slice?
   bool LastSlice() const {return lastSlice;};

   ~FitsOutSlice();


};


//! Several fits images of the same size to be processed in parallel.
/*! This class is to be used as a handler for FitsSlice's from files
having the same sizes, that one wants to process in parallel (typically to
average or sum  images). Slices are cut along the second index of the whole image 
and have the same size in x as the image. 
See \ref example_slices for an example. */

class FitsParallelSlices : public vector<FitsSlice*>  /* should set it private I guess ?? */
{

private:
  int sliceSize;
  const int overlap;
  int imageSizeX, imageSizeY;

public :


  //! constructor.
  FitsParallelSlices(const int SliceSize, const int Overlap = 0);
  ~FitsParallelSlices();
  //! Adds a file in the list. Loads the first slice.
  bool AddFile(const string &FileName);


  //! load next slice for all involved files. Returns 0 if already on last slice.
  int LoadNextSlice();
  //! size of the current slice size. Different from the constructor value for last slice.
  int SliceSize() const {return sliceSize;};
  //! return the index in the whole image of the row SliceJ in the current slice.
  int ImageJ(const int SliceJ)  const {return (*begin())->ImageJ(SliceJ); };
  //! are we in the last slice
  bool LastSlice() const {return (*begin())->LastSlice();}
  //!
  int NFiles() const {return size();}
}; 

typedef vector<FitsSlice*>::iterator       FitsSliceIterator;
typedef vector<FitsSlice*>::const_iterator FitsSliceCIterator;


//! To traverse in parallel one "in" and one "out" image
class FitsInOutParallelSlices  
{
public : 
    // no point to set them private
  FitsSlice in;
  FitsOutSlice out;

 public:
  FitsInOutParallelSlices(const std::string &InName,
			  const std::string &OutName,
			  const int YSliceSize = 100,
			  const int Overlap = 0);

  //! load next slice of input and write current slice of output
  int LoadNextSlice();
  //! size of the current slice size. Different from the constructor value for last slice.
  int SliceSize() const {return in.SliceSize();};
  //! return the index in the whole image of the row SliceJ in the current slice.
  int ImageJ(const int SliceJ)  const {return in.ImageJ(SliceJ); };
  //! are we in the last slice
  bool LastSlice() const {return in.LastSlice();}
  
};


#endif
