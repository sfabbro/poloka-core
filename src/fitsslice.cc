
#include "fitsslice.h"
#include  <fitsio.h>

/************** info when we get cfitsio errors ************/

#define CHECK_STATUS(status, message, return_statement) \
if (status)\
  {\
  cfitsio_report_error(status, message);\
  cout << " for file " << this->fileName << endl;\
  return_statement;\
  }

static void cfitsio_report_error(int status, const char *Message)
{
cerr << Message ;
fits_report_error(stderr, status);
}


/************************** FitsSlice ******************************************/

//#define DEBUG_SLICES

FitsSlice::FitsSlice(const string FileName, const int SliceSize, const int Overlap) : 
            FitsHeader(FileName), Image(0,0)
{
  overlap = Overlap;
  
  int nx = KeyVal("NAXIS1");
  nyTotal = KeyVal("NAXIS2");
  sliceSize = min(SliceSize,nyTotal);
  Image::allocate(nx, sliceSize);
  read_pixels(0, sliceSize, data);
  ySliceStart = 0;
}


int FitsSlice::read_pixels(const int StartRow, int &NRows, Pixel* Where)
{
#ifdef DEBUG_SLICES
  cout << " read requested for file : " << FileName() << " from row " << StartRow << " for " << NRows << endl;
#endif
  int startPix = StartRow*Nx();
  lastSlice = (StartRow+NRows >= nyTotal);
  int endRow = min(StartRow+NRows, nyTotal);
  NRows = endRow - StartRow; /* to return the actual size of the read chunk */
  int nPix =  NRows* Nx();
  int status = 0;
  int anynull = 0;
  float nullval = 0;
  /* fits_read_img uses the fortran like numberring: fits pixel is pixel "1" */
  fits_read_img(fptr, TFLOAT, startPix+1, nPix, &nullval, Where,  &anynull, &status);
  CHECK_STATUS(status," FitsSlice::read_pixels", );
  return (!status);
}
  

int FitsSlice::LoadNextSlice()
{
  if (lastSlice) return 0;
  ySliceStart += sliceSize - overlap;
  int nx = Nx();
  int n_pix_overlap = overlap * nx;
  
  /* take the part at the end (corresponding to the overlap between current slice
     and next one and place it at the beginning) */
  int nRowRead = sliceSize-overlap;
  int nbyte_copy = n_pix_overlap*sizeof(Pixel);
  if (nbyte_copy) memcpy(data, data + nRowRead*nx, nbyte_copy);
  int status = read_pixels(ySliceStart+overlap, nRowRead, data + n_pix_overlap);
#ifdef DEBUG_SLICES
  cout << " actually read  " << nRowRead << " rows from file " << FitsHeader::FileName() << endl;
#endif
  // nRowRead differs on output for the last slice : assign the sliceSize: 
  sliceSize = overlap + nRowRead;
  return status;
}


/************************** FitsParallelSlices **************************/


FitsParallelSlices::FitsParallelSlices(const int SliceSize, const int Overlap) : 
  sliceSize(SliceSize), overlap(Overlap) 
{
  imageSizeX=0; imageSizeY = 0;
}

FitsParallelSlices::~FitsParallelSlices()
{
  for (FitsSliceIterator it = begin() ;  it != end(); ++it)
    {
      delete *it;
    }
}

bool FitsParallelSlices::AddFile(const string &FileName)
{
  FitsSlice *newFile = new FitsSlice(FileName, sliceSize, overlap);
  // check existence and 'fits'ness
  if (!newFile->IsValid())
    {
      delete newFile;
      return false;
    }
  // check sizes.
  int thisSizeX,thisSizeY;
  newFile->ImageSizes(thisSizeX, thisSizeY);
  if (!imageSizeX) imageSizeX =  thisSizeX;
  if (!imageSizeY) imageSizeY =  thisSizeY;
  if (imageSizeX !=  thisSizeX || imageSizeY !=  thisSizeY)
  {
    cerr << " FitsParallelSlices::AddFile " << FileName << "does not have the same size as file previously added " << endl;
    delete newFile;
    return false;
  }
// OK.
  push_back(newFile);
  return true;
}

int FitsParallelSlices::LoadNextSlice()
{
int global_status = 1;
 for (FitsSliceIterator it = begin() ;  it != end(); ++it)
   {
     FitsSlice *slice = *it;
     int status = slice->LoadNextSlice();
     if (!status && !slice->LastSlice()) cerr << " problems when going to next slice for file " << slice->FileName() << endl;
     global_status &= status;
     sliceSize = slice->SliceSize();
 }
 return global_status;
}

/*
  comments for the doc:
  - the name  "show" is supposed to be a file name. I could not get it
  to be found.
  - if you replace \code by \verbatim, the whole example disappears.
  - with \code, the comments are not displayed as comment in the html output.
*/
/* \file */

/*! \page example_slices Usage of FitsParallelSlices

  How to use  FitsParallelSlices
  to compute (e.g.) a clipped median of  a set of images, to build, e.g. a superflat:
  
  \code
  ...
  FitsParallelSlices slices(sliceSize,overlap);
  
  // // insert the files ...
  for (int i=3; i<nargs; ++i)
  {
  string fileName = args[i];
  slices.AddFile(fileName);
  }
  
//  // process 
Pixel *pixelvalues = new Pixel[slices.NFiles()];
do
  {
    for (int j=0; j<slices.SliceSize(); j++) for (int i=0; i<Nx_Image; i++) 
      {
	for (int k=0; k<nImages; k++) pixelValues[k] = (*slices[k])(i,j);
	float mean;
	mean = FArrayMedian(pixelValues, nImages);
	int j_true = slices.ImageJ(j);
	(*Flat)(i,j_true) = mean;
      }
  }  
  while (slices.LoadNextSlice());
delete [] pixelValues;
...
}
\endcode
*/

