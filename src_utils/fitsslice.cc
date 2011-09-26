
#include "fitsslice.h"
#include  <fitsio.h>


/************** info when we get cfitsio errors ************/

#define CHECK_STATUS(status, message, return_statement) \
if (status)\
  {\
  cfitsio_report_error(status, message);\
  std:: cout << " for file " << this->FileName() << std::endl;\
  return_statement;\
  }

static void cfitsio_report_error(int status, const char *Message)
{
  std::cerr << Message ;
  fits_report_error(stderr, status);
}


/************************** FitsSlice ******************************************/

//#define DEBUG_SLICES

FitsSlice::FitsSlice(const string FileName, const int SliceSize, const int Overlap) : 
            FitsHeader(FileName), Image(0,0)
{
  overlap = Overlap;
  if (overlap > SliceSize)
    {
      cout << " ERROR : FitsSlice : requested Overlap > SliceSize , this is impossible " << endl;
      abort();
    }
  
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
  lastSlice = (StartRow+NRows >= nyTotal);
  int endRow = min(StartRow+NRows, nyTotal);
  NRows = endRow - StartRow; /* to return the actual size of the read chunk */

  // /* fits_read_img uses the fortran like numberring: fits pixel is pixel "1" */
  // fits_read_img(fptr, TFLOAT, startPix+1, nPix, &nullval, Where,  &anynull, &status);
  /* to accomodate  compressed images (for which fits_read_img does not work properly,
     use our own reading routine */
  int status = read_image(0, StartRow, nx, endRow, Where);

  CHECK_STATUS(status," FitsSlice::read_image", );
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
  int newSliceSize = overlap + nRowRead;
  // it is in fact cleaner to assign the ny attribute of the image:
  ny = sliceSize = newSliceSize;
  return status;
}


/******************* FitsOutSlice ************************************/
FitsOutSlice::FitsOutSlice(const string FileName, 
			   const int Nx, const int Ny,
			   const int SliceSize, const int Overlap) 
  : FitsHeader(FileName,RW) , Image(0,0), bitpix(0)
{
  overlap = Overlap;  
  int nx = Nx;
  nyTotal = Ny;
  sliceSize = SliceSize;
  Image::allocate(nx, sliceSize);
  ySliceStart = 0;
}

FitsOutSlice::FitsOutSlice(const string FileName, const FitsHeader &ModelHead,
			   const int SliceYSize, const int Overlap) 
  : FitsHeader(ModelHead, FileName) , Image(0,0), bitpix(0)
{
  overlap = Overlap;  
  int nx = KeyVal("NAXIS1");
  nyTotal = KeyVal("NAXIS2");
  sliceSize = SliceYSize;
  Image::allocate(nx, sliceSize);
  ySliceStart = 0;
  if (sliceSize >= nyTotal)
    {
      sliceSize = nyTotal;
      lastSlice = true;
    }
}

int FitsOutSlice::WriteCurrentSlice()
{
  int nx = Nx();
  int n_pix_overlap = overlap * nx;
  
  int nRowWrite = sliceSize;
  int status = write_pixels(ySliceStart, nRowWrite, data);
  /* take the part at the end (corresponding to the overlap between
     current slice and next one and place it at the beginning) */
  //lastslice updated at the end of write_pixels
  // if lastslice is true, then we've just written the lastslice
  if (! LastSlice())
    {
  int nbyte_copy = n_pix_overlap*sizeof(Pixel);
  if (nbyte_copy) memcpy(data, data + (nRowWrite-overlap)*nx, nbyte_copy);
  // clear the remainder
  Pixel *zeroStart = data + n_pix_overlap;
  int nbyte_zero = (end() - zeroStart)*sizeof(Pixel);
  memset(zeroStart, 0, nbyte_zero);

  // where we stand now (what the pixel buffer now holds):
  ySliceStart += sliceSize - overlap;

#ifdef DEBUG_SLICES
  cout << " actually written " <<  nRowWrite << " rows " << endl;
  cout << " zeroed  " <<  nbyte_zero/(nx*sizeof(Pixel)) << " rows " << endl;
#endif
    }
  // nRowRead differs on output for the last slice : assign the sliceSize: 
  //  sliceSize = overlap + nRowRead;
  return status;
}

int FitsOutSlice::write_pixels(const int StartRow, int &NRows, Pixel *Where)
{
  if (bitpix == 0)
    {
      // to avoid reading them at each pixel chunk.
      bitpix = KeyVal("BITPIX");
      bscale = KeyVal("BSCALE");
      bzero =  KeyVal("BZERO");
      create_image(bitpix, nx, nyTotal);
    }
#ifdef DEBUG_SLICES
  cout << " write requested for file : " << FileName() << " from row " << StartRow << " for " << NRows << endl;
#endif
  int startPix = StartRow*Nx();
  int endRow = min(StartRow+NRows,nyTotal);
  lastSlice = (StartRow+NRows >= nyTotal);
  NRows = endRow-StartRow;
  int nPix =  NRows * Nx();
  int status = 0;

  if (nPix) status = write_image(startPix, nPix,  Where, 
				 nx, bitpix, bscale, bzero);
  CHECK_STATUS(status," FitsOutSlice::write_pixels", );
  return (!status);
}


FitsOutSlice::~FitsOutSlice()
{
  WriteCurrentSlice(); // can be done several times (only for the last slice!)
  std::cout << " FitsOutSlice : writen " 
	    << int(KeyVal("NAXIS1")) << 'x' << int(KeyVal("NAXIS2")) << ' '
	    << FileName() << endl
	    << " with BITPIX=" << int(KeyVal("BITPIX"))
    //	    << "BSCALE=" << double(KeyVal("BSCALE")) 
    //	    << "BZERO=" << double(KeyVal("BZERO"))
	    << std::endl;
  // don't close the output file here : ~Fitsheader does it.
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



/**************  class FitsInOutParallelSlices  ************/

/* The added value w.r.t FitsSlice and FitsOutSlice is just that
FitsInOutParallelSlices class forces SliceSize and Overlap to be the
same for input and ouput.  it also sequences read(in) and Write(out)
in a single call.  If you say it is not much, I won't rush into an
argument!
*/


FitsInOutParallelSlices::FitsInOutParallelSlices(const std::string &InName,
			  const std::string &OutName,
			  const int YSliceSize,
			  const int Overlap)
  : in(InName, YSliceSize, Overlap), out(OutName, in, YSliceSize, Overlap)
{
}

   

int FitsInOutParallelSlices::LoadNextSlice()
{
  return (in.LoadNextSlice() && out.WriteCurrentSlice());
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

