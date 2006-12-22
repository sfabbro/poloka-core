/* 
 * $Source: /cvs/snovae/toads/poloka/flat/fitsimagearray.cc,v $
 * $Revision: 1.2 $
 * $Author: guy $
 * $Date: 2006/12/22 13:35:40 $
 * $Name:  $
 */

#include <iostream>
#include <fitsio.h>
#include "fitsimagearray.h"
#include "fileutils.h"
#include <time.h>

// uncomment this to get a lot of output
#define DEBUG 

// get cvs version of the code
#define CVSVERSION "$Revision: 1.2 $"

// ============ STUFF copied in fitsimage =======================
static string  local_time()
{
time_t now = time(NULL);
//struct tm *date = localtime(&now);
char *ptime = ctime(&now);
// eventually remove \n at the end...
char *plast = ptime + strlen(ptime)-1;
if (*plast == '\n') *plast = '\0';
return string(ptime);
}


#define CHECK_STATUS(status, message, return_statement) \
if (status)\
  {\
  cfitsio_report_error(status, message);\
  cerr << " for file " << this->fileName << endl;\
  return_statement;\
  }
static void cfitsio_report_error(int status, const string &Message)
{
cerr << Message ;
fits_report_error(stderr, status);
}
// ============ end of STUFF copied in fitsimage =======================


FitsImageArray::FitsImageArray(const string &FileName, const FitsFileMode Mode) :
  FitsImage(FileName, Mode){
#ifdef DEBUG
  cout << "Creating a FitsImageArray with an existing FITS file : " << FileName.c_str() << endl;
#endif
  fMainHeader=0;
  fValid=true;
  if (!IsValid()) {
    cout << "ERROR in FitsImageArray::FitsImageArray : FitsImage is not valid" << endl;
    return; 
  }  
  
  bool extended = true;

  // Check whether this file already has several images, otherwise quit
  if (! HasKey("EXTEND") || ! KeyVal("EXTEND")) {
    if(FileMode()==RW) {
      cout << "ERROR in FitsImageArray::FitsImageArray :  this image is not extended, please use another creator to open it" << endl;
      fValid = false;
      return;
    }else{
      extended = false;
    }
  }
  
  // check if the first HDU contains an image or not
  if (! HasKey("NAXIS")) {
    cout << "ERROR in FitsImageArray::FitsImageArray : no key NAXIS" << endl;
    fValid = false;
    return;
  }
  int naxis = (int)KeyVal("NAXIS");  
  if(naxis!=0 && naxis!=2 ) {
    cout << "ERROR in FitsImageArray::FitsImageArray : NAXIS = " << naxis << endl;
    fValid = false;
    return;
  }
  fImageInFirstHDU = (naxis==2);
  
  Init(extended);
}

FitsImageArray::FitsImageArray(const string &FileName, const FitsHeader& imageheader, bool ImageInFirstHDU) :
  FitsImage(FileName,imageheader){
#ifdef DEBUG
  cout << "Creating a new FitsImageArray with an existing FitsHeader" << endl;
#endif
  fMainHeader=0;
  fValid=true;
  if(!IsValid())
    return;
  fImageInFirstHDU = ImageInFirstHDU;
  Init(false);
}

FitsImageArray::FitsImageArray(const string &FileName, const Image& image, bool ImageInFirstHDU) :
  FitsImage(FileName,image) {
#ifdef DEBUG
  cout << "Creating a new FitsImageArray with an existing Image" << endl;
#endif 
  fMainHeader=0;
  fValid=true;
  if(!FitsImage::IsValid())
    return;
  fImageInFirstHDU = ImageInFirstHDU;
  Init(false);
}

FitsImageArray::Status FitsImageArray::Init(bool already_extended) {
  fSelectedCriterion = FILTER|SIZE;
  
#ifdef DEBUG
  cout << "Init" << endl;
  cout << "fImageInFirstHDU = " << fImageInFirstHDU << endl;
#endif

  if(!already_extended && FileMode()==RW) {
    AddOrModKey("EXTEND",true);
    AddOrModKey("NEXTEND",1);
    AddCommentLine("Processed by FitsImageArray in TOADS");
    AddCommentLine("This is the main header");
    char comment[80];
    sprintf(comment,"CVS version of fitsimagearray.cc : %s",CVSVERSION);
    AddCommentLine(comment);
  }
  
  // save the main header, this creator copies the header of *this to a new file
  // warning: this may be quite a mess if several people us the same file 
  // so we add local time to the file name (reduce probability of problems) 
  char name[256];
  sprintf(name,"mainheader_%f.fits",(float)time(NULL));
  fMainHeader = new FitsHeader(*this,name);
  fMainHeader->EnableWrite(false); // do not save it when destroyed

#ifdef DEBUG
  cout << "main header name is = \"" << name << "\"" << endl; 
#endif
  
  if(!already_extended && !fImageInFirstHDU && FileMode()==RW) {
    cout << "Writing the first HDU containing the first header (and no image)" << endl;
    AddOrModKey("NAXIS",0); // no image here !!
    if(HasKey("NAXIS1"))
      RmKey("NAXIS1");
    if(HasKey("NAXIS2"))
      RmKey("NAXIS2"); 
    if(HasKey("EXTNAME"))
      RmKey("EXTNAME");
    AddOrModKey("EXTVER",0);
    // writing a NULL image, I don't know if this is standard
  }
  
  return OK;
}


FitsImageArray::~FitsImageArray() {
#ifdef DEBUG
  cout << "Destroying a FitsImageArray" << endl;
#endif
  
  // to prevent the FitsImage destructor to rewrite the file
  fileModeAtOpen = RO;
  
  if(fMainHeader!=0) {
    delete fMainHeader;
  }
}

FitsImageArray::Status FitsImageArray::Write(bool force_bscale) {
#ifdef DEBUG
  cout << "Writing a FitsImageArray <==> writing the FitsImage" << endl;
#endif

  if (!writeEnabled || FileMode() != RW)
    return ERRORWRITE;
  
  int status = FitsImage::Write(force_bscale);
  if(status==1)
    return OK; 
  else
    return FAILURE;
}

FitsImageArray::Status FitsImageArray::SplitAndWrite(const string &directory,int HDUMax) {
#ifdef DEBUG
  cout << "Writing several FitsImage" << endl;
#endif
  string dir,base,ext;
  SplitName(dir, base, ext);
  if(directory.size()>0)
    dir = directory;
  dir = dir+"/"+base;
  
  if(!IsDirectory(dir))
    MKDir(dir.c_str());
  
  base = "ccd_";
  ext = "fits";
  FitsHeader mainHeader(*((FitsHeader*)this), dir+"/mainheader.fits");
  
  char schip[256];
  char outFileName[512];
  const char *file_name_format = "%s/%s%02d.%s";

  int nhdu = NHDU();
  if(HDUMax==0)
    HDUMax=nhdu;
  
  cout << "FitsImageArray::SplitAndWrite nhdu = " << nhdu << endl;
  
  for (int hdu= 2; (hdu<= nhdu && hdu<= HDUMax); ++hdu) {
    if (At(hdu) != OK) {
      cout << "WARNING in FitsImageArray::SplitAndWrite : NDU limit reached" << endl;
      break;
    }
    int chip = KeyVal("TOADCHIP",true);
    if(chip==0)
      chip=hdu-2;
    
    cout << " considering image from chip " << chip << endl;
    sprintf(schip,"%d",chip);
    sprintf(outFileName,file_name_format,dir.c_str(),base.c_str(),chip, ext.c_str());
    cout << outFileName << endl;
    FitsHeader outFile(*((FitsHeader*)this), outFileName);
    outFile.Append_LowPriority(mainHeader);
    if (outFile.HasKey("EXTEND")) 
      outFile.ModKey("EXTEND",false);
    if (CopyDataTo(outFile)!=1) {
      cout << "ERROR in FitsImageArray::SplitAndWrite : CopyDataTo failed " << endl;
      return FAILURE;
    }
  }
  return OK; 
}

bool FitsImageArray::CheckHeader(FitsHeader &header) {
  if(!header.IsValid()) {
     cout << "ERROR in FitsImageArray::CheckHeader, FITS image is not valid" << endl;
     return false;
  }
  if(fMainHeader==0) {
    cout << "ERROR in FitsImageArray::CheckHeader, fMainHeader=0" << endl;
    return false;
  }
  bool result = true;
  if(fSelectedCriterion & SIZE)
    result &= fMainHeader->SameImageSizes(header);
  if(fSelectedCriterion & CHIP)
    result &= fMainHeader->SameChip(header);
  if(fSelectedCriterion & FILTER)
    result &= fMainHeader->SameFilter(header);
  
  if(fSelectedCriterion & EXTENDED)
    if( header.HasKey("EXTEND") && (bool)(header.KeyVal("EXTEND")) ) {
      cout << "WARNING in FitsImageArray::CheckHeader : input file is already extended, I don't add it" << endl;
      return false;
    }
  
  return result;
}

FitsImageArray::Status FitsImageArray::Append(const Image& img) {  
  
  if (!writeEnabled || FileMode() != RW) {
    cout << "Cannot Append an Image to the file (ReadOnly Mode)" << endl;
    return ERRORWRITE;
  }
  
  
#ifdef DEBUG
  cout << "Append an Image" << endl;
  cout << "initial NHDU = " << NHDU() << endl;
#endif
  
  bool writefirsthdu = ( NHDU()==0 && fImageInFirstHDU );
  
  // see cfitio user's guide version 2.4  
  // p34 for fits_copy_hdu
  // and p94 for fits_write_img
  // and p39 for fits_create_img
  // and p40 for fits_write_pix
  
  
  int status;
  if(writefirsthdu) { // writes image in the first HDU
    
    // Fill the header according to img content
    FillHeader(img);
    
#ifdef DEBUG  
    cout << "fits_write_img ..."<< endl;
#endif   
    status=0;
    fits_write_img(fptr, TFLOAT, 1, img.Nx()*img.Ny(), img.data, &status);
    
  } else {
    int naxis = 2;
    int bitpix = 16;
    long naxes[2];
    naxes[0]=img.Nx(); // NAXIS1 // maybe something stupid there
    naxes[1]=img.Ny(); // NAXIS2
    
#ifdef DEBUG  
    cout << "fits_create_img ..."<< endl;
#endif
    status=0;
    fits_create_img(fptr,bitpix,naxis,naxes,&status);
    CHECK_STATUS(status," Append ", );
    
#ifdef DEBUG  
    cout << "intermediate NHDU = " << NHDU() << endl;
#endif
    
    // Fill the header according to img content
    FillHeader(img);


    long fpixel[2]; // maybe something stupid there again
    fpixel[0]=1;
    fpixel[1]=1;
    
#ifdef DEBUG  
    cout << "fits_write_pix ..."<< endl;
#endif
    status = 0; 
    fits_write_pix(fptr,TFLOAT,fpixel,img.Nx()*img.Ny(),img.data,&status);
  }
  
  CHECK_STATUS(status," FitsImageArray::Append ", );
#ifdef DEBUG  
  cout << "final NHDU = " << NHDU() << endl;
#endif
  

if(status)
  return FAILURE;
else
return OK;
}



FitsImageArray::Status FitsImageArray::Append(const string& inputfilename) {
  
#ifdef DEBUG
  cout << "Append an Image" << endl;
  cout << "initial NHDU = " << NHDU() << endl;
#endif
  
  if (!writeEnabled || FileMode() != RW) {
    cout << "Cannot Append an Image to the file (ReadOnly Mode)" << endl;
    return ERRORWRITE;
  }

  FitsHeader temp_header(inputfilename,RO);
  if(!CheckHeader(temp_header))
    return FAILURE;
  
  // see cfitio user's guide version 2.4  p34
  // and p94 for fits_write_img
  int status = 0;
  fits_copy_hdu(temp_header.fptr,fptr,0,&status);
  CHECK_STATUS(status," Append ", );
  
#ifdef DEBUG  
  cout << "final NHDU = " << NHDU() << endl;
#endif 

  if(status)
    return FAILURE;
  else
    return OK; 
}

void FitsImageArray::FillHeader(const Image& image, bool force_bscale) {

#ifdef DEBUG  
  cout << "FitsImageArray::FillHeader ... " << endl;
#endif
  // hard input for the moment
  int bitpix = 16;
  
  // already saved by create image: SIMPLE,BITPIX,NAXIS,NAXIS1,NAXIS2
  AddOrModKey("BSCALE", 1.0);
  AddOrModKey("NAXIS", 2);
  AddOrModKey("NAXIS1", image.Nx());
  AddOrModKey("NAXIS2", image.Ny());
  AddOrModKey("BITPIX", bitpix);
  AddOrModKey("FIACOUNT", NHDU());
  
  // time
  string stime = local_time();
  char * time =  (char *) stime.c_str();
  AddOrModKey("WRITEDAT", time, " when this file was written"); 
  char name[80];
  
  // default extension EXTNAME & EXTVER
  int extension = 0;
  int nhdu = NHDU();
  if(fImageInFirstHDU)
    extension = nhdu;
  else
    extension = nhdu-1;
  sprintf(name,"im%03d",extension);
  if( ! HasKey("EXTNAME") ) // keeps existing EXTNAME, overwrite EXTVER
    AddKey("EXTNAME",name);

  // DO NOT MODIFY EXTVER which is used by KeyVal("TOADCHIP") !!!!
  //  AddOrModKey("EXTVER",extension);
  
  


  if (bitpix == 16)
    {
      double BScale = 1. ;
      double BZero = 0. ;
      if (HasKey("KEEPZERO") && bool(KeyVal("KEEPZERO")))
	{// compute Bscale and BZero so that 0 remain 0 (used for weight maps)
	  Pixel min,max;
	  image.MinMaxValue(&min,&max);
	  if (max != 0)
	    {
	      double short_max = SHRT_MAX;
	      double span = short_max;
	      BScale = max/span;
	      // this in fact assumes that the content is >=0
	      BZero = 0;
	      //	    cout << BZero/BScale << endl;
	    }
	  cout << nx << "X"<< ny << " 16 bits" 
	       << "   Min=" << min << "  Max=" << max
	       << " BSCALE=" << BScale << "   BZERO=" << BZero << "(forced)" << endl;
	}	
      else if (!force_bscale)
	{
	  /* CFITSIO DOES NOT COMPUTE automatically BSCALE and BZERO 
	     according to BITPIX in most of the cases */
	  Pixel min,max;
	  image.MinMaxValue(&min,&max);
	  double short_min =  SHRT_MIN +2 ;
	  double short_max =  SHRT_MAX -2 ;
	  double span = short_max -  short_min ;
	  BScale = (max - min)/span;
	  // BScale is used as the denominator of a division...
	  if (BScale == 0) BScale = 1; // constant image;
	  BZero = min - BScale * (short_min);
	  cout << nx << "X"<< ny << " 16 bits" 
	       << "   Min=" << min << "  Max=" << max
	       << " BSCALE=" << BScale << "   BZERO=" << BZero << endl;
	}
      else
	cout << "saved with BSCALE = 1. and BZERO = 0. " << endl ;
      ModKey("BSCALE",BScale);
      ModKey("BZERO", BZero);
    }  
  
  Flush();
}

FitsImageArray::Status FitsImageArray::Next() {
#ifdef DEBUG
  cout << "FitsImageArray::Next() for NHDU=" << NHDU() << endl;
#endif
  if(fCurrentHDU+1>NHDU()) { // a priori out of range
    cout << "ERROR in FitsImageArray::Next() out of range" << endl;
    return OUTOFBOUNDS;
  }
  int status = 0;
  fits_movrel_hdu(fptr, 1, NULL, &status);
  CHECK_STATUS(status,"FitsImageArray::Next",);
  fCurrentHDU++;
  if(status==0)
    return OK; 
  else
    return OUTOFBOUNDS;
}

FitsImageArray::Status FitsImageArray::At(int HDU){
#ifdef DEBUG
  cout << "FitsImageArray::At(" << HDU << ")" << endl;
#endif
  
  if(HDU<1 || HDU>NHDU()) { // a priori out of range
    cout << "ERROR in FitsImageArray::At(" << HDU << ") out of range (NHDU=" << NHDU() << ")"<< endl;
    return OUTOFBOUNDS;
  }
  if(HDU==fCurrentHDU) { // nothing to do
#ifdef DEBUG
    cout << "Nothing to do cause (HDU==fCurrentHDU" << endl;
#endif
    return OK;
  }
  
#ifdef DEBUG
  cout << "fits_movabs_hdu ..." << endl;
#endif
  int status = 0;
  fits_movabs_hdu(fptr, HDU, NULL, &status);
  //fits_movrel_hdu(fptr, HDU-fCurrentHDU, NULL, &status);
  CHECK_STATUS(status,"FitsImageArray::At",);
  if(status!=0)
    return OUTOFBOUNDS;
  fCurrentHDU = HDU;
  
  return ReadImage();
}
  

FitsImageArray::Status FitsImageArray::ReadImage() {

#ifdef DEBUG
  cout << "FitsImageArray::ReadImage()" << endl;
#endif
  if(!HasKey("NAXIS") || ((int)KeyVal("NAXIS"))!=2 ) {
    return FitsImageArray::EMPTY;
  }
  
  int nx_in_fits = KeyVal("NAXIS1",true);
  int ny_in_fits = KeyVal("NAXIS2",true);
  
  if( nx != nx_in_fits || ny != ny_in_fits ) {
    cout << "WARNING in FitsImageArray::ReadImage nx!=NAXIS1 or ny!=NAXIS2" << endl;
    nx = nx_in_fits;
    ny = ny_in_fits;
    allocate(nx,ny);
  }
  
  float nullval = 0;
    int anynull;
#ifdef DEBUG
  cout << "fits_read_img ..." << endl;
#endif
  int status = 0;
  fits_read_img(fptr, TFLOAT, 1, nx*ny, &nullval,  data, &anynull, &status);
  CHECK_STATUS(status," FitsImageArray::ReadImage ", );
  if(status!=0)
    return FAILURE;
  else
    return OK;
}
  

void FitsImageArray::SplitName(string &Dir, string &Base, string &Type) {
  Dir = DirName(fileName);
  Base = BaseName(fileName);
  char name[256];
  sprintf(name,"%s",Base.c_str());
  char *gz = strstr(name,".gz");
  if (gz) *gz = '\0';
  char *dot = strrchr(name,'.');
  if (dot)
    {
      Type = dot+1;
      *dot='\0';
      Base = name;
    }
  else
    {
      Type = "";
    }
}

bool FitsImageArray::IsValid() {
  bool headerisvalid = FitsHeader::IsValid();
  if(!headerisvalid) {
#ifdef DEBUG
    cout << "Header is ot valid" << endl;
#endif
    return false;
  }
  return fValid;
}


string GetVersion() {
  string version = CVSVERSION;
  return version;
}
