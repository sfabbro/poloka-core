#include <iostream>
#include <string>
#include <vector>

#include <fitsio.h>

#include <poloka/fileutils.h>
#include <poloka/fitsimage.h>

static bool has_key(fitsfile *fptr, char *key) {
  int status=0;
  char a_C_string[100];
  fits_read_key(fptr, TSTRING, key, a_C_string, NULL, &status);
  return (status != KEY_NO_EXIST);
}

static void usage(const char * prog) {
  cerr << "Usage: " << prog << " [OPTION]... FITS FITS...\n"
       << "Merge FITS files\n\n"
       << "    -c          : use rice compression\n"
       << "    -d DIRECTORY: merge files into DIRECTORY\n\n";
  exit(EXIT_FAILURE);
}

static int convert_fits_file(const string &inFile, const string &outFile) {
  fitsfile *infptr, *outfptr;
  int status=0;

   /* Open the input file and create output file */
  fits_open_file(&infptr, inFile.c_str(), READONLY, &status);
  fits_create_file(&outfptr, outFile.c_str(), &status);

  if (status != 0) {
    fits_report_error(stderr, status);
    return status;
  }


  /* Get the current HDU position */
  int hdupos;
  fits_get_hdu_num(infptr, &hdupos);
  int first_hdu = hdupos;

  /* Main loop through each extension */
  for (; !status; hdupos++) {
    int hdutype;
    long int naxes[9];
    int bitpix;
    long int totpix;
    int naxis=0;
    
    fits_get_hdu_type(infptr, &hdutype, &status);
    
    if (hdutype == IMAGE_HDU) {
      
      /* get image dimensions and total number of pixels in image */
      for (int ii = 0; ii < 9; ii++)
	naxes[ii] = 1;

      fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);

      totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4]
	     * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    }

      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0) { 

	/* just copy tables and null images */
	if (! (hdutype == IMAGE_HDU && naxis == 0)) // no .. not NULL Images (not in original imcopy)
          fits_copy_hdu(infptr, outfptr, 0, &status);

      } else {

	/* Explicitly create new image, to support compression */
	fits_create_img(outfptr, bitpix, naxis, naxes, &status);

          /* copy all the user keywords (not the structural keywords) */
	  int nkeys;
          fits_get_hdrspace(infptr, &nkeys, NULL, &status); 

          for (int ii = 1; ii <= nkeys; ii++) {
	    char card[100];
	    fits_read_record(infptr, ii, card, &status);
	    if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	      fits_write_record(outfptr, card, &status);
          }

	  // DEBUG
	  {
	    int out_hdu;
	    fits_get_hdu_num(outfptr, &out_hdu);
	    cout << " out_hdu " << out_hdu << endl;
	  }

	  /* MODIFICATION of imcopy : merge the main header into the extension header */
	  int current_hdu;
	  fits_get_hdu_num(infptr, &current_hdu);  /* Get the current HDU position */
	  fits_movabs_hdu(infptr, first_hdu, NULL, &status);
          /* copy all the user keywords (not the structural keywords) */
          fits_get_hdrspace(infptr, &nkeys, NULL, &status); 

          for (int ii = 1; ii <= nkeys; ii++) {
	    char card[100];
	    fits_read_record(infptr, ii, card, &status);
	    if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
	      char key[80];
	      strncpy(key, card, 30);
	      char *equal=strchr(key,'=');
	      if (equal) { // if not, it is a comment/history
		*equal='\0';
		if (has_key(outfptr, key)) continue;
	      }
	      fits_write_record(outfptr, card, &status);
	    }
	  }

	  // DEBUG
	  {
	    int out_hdu;
	    fits_get_hdu_num(outfptr, &out_hdu);
	    cout << " out_hdu " << out_hdu << endl;
	  }

	  fits_movabs_hdu(infptr, current_hdu, NULL, &status);

	  int datatype = 0;    
          switch (bitpix) {
	  case BYTE_IMG:
	    datatype = TBYTE;
	    break;
	  case SHORT_IMG:
	    datatype = TSHORT;
	    break;
	  case LONG_IMG:
	    datatype = TINT;
	    break;
	  case FLOAT_IMG:
	    datatype = TFLOAT;
	    break;
	  case DOUBLE_IMG:
	    datatype = TDOUBLE;
	    break;
          }

          int bytepix = abs(bitpix) / 8;

          long int npix = totpix;
          int iteration = 0;

          /* try to allocate memory for the entire image */
          /* use double type to force memory alignment */
          double *array = (double *) calloc(npix, bytepix);

          /* if allocation failed, divide size by 2 and try again */
          while (!array && iteration < 10)  {
	    iteration++;
	    npix = npix / 2;
	    array = (double *) calloc(npix, bytepix);
          }
	  
          if (!array)  {
	    cerr << "convert_fits_file: memory allocation error\n";
	    return EXIT_FAILURE;
          }

          /* turn off any scaling so that we copy the raw pixel values */
	  double bscale=1;
	  double bzero=0;
          fits_set_bscale(infptr,  bscale, bzero, &status);
          fits_set_bscale(outfptr, bscale, bzero, &status);

          int first = 1;
          while (totpix > 0 && !status) {
	    double nulval = 0.;
	    int anynul;
             /* read all or part of image then write it back to the output file */
             fits_read_img(infptr, datatype, first, npix, 
			   &nulval, array, &anynul, &status);
	     
             fits_write_img(outfptr, datatype, first, npix, array, &status);
             totpix = totpix - npix;
             first  = first  + npix;
          }
          free(array);
	  break;
      }

      /* try to move to next HDU */
      fits_movrel_hdu(infptr, 1, NULL, &status);
  }
  
  /* Reset after normal error */
  if (status == END_OF_FILE) status = 0;

  fits_close_file(outfptr,  &status);
  fits_close_file(infptr, &status);

  /* if error occurred, print out error message */
  if (status)
    fits_report_error(stderr, status);
  return status;
}

int main(int nargs, char **args) {

  string outDir;
  vector<string> fileNames;
  bool riceCompression=false;

  for (int i=1; i< nargs; i++) {
    char *arg = args[i];
    if (arg[0] == '-')
      switch (arg[1]) {
      case 'd' : outDir = AddSlash(args[++i]); break;
      case 'c' : riceCompression=true; break;
      default: usage(args[0]);
      }
    else
      fileNames.push_back(args[i]);
  }

  if (outDir.empty()) {
    if (fileNames.size() != 2) usage(args[0]);
    if (riceCompression)
      return convert_fits_file(fileNames[0], fileNames[1]+"[COMPRESS]");
    else 
      return convert_fits_file(fileNames[0], fileNames[1]);
  }
  
  MKDir(outDir.c_str());

  int status = 0;
  for (size_t k=0; k< fileNames.size(); ++k) {
    string outputName;
    if (riceCompression)
      outputName= outDir+CutExtension(BaseName(fileNames[k]))+".fz[COMPRESS]";
    else
      outputName= outDir+BaseName(fileNames[k]);	  
    status += convert_fits_file(fileNames[k], outputName);
  }

  return status? EXIT_SUCCESS : EXIT_FAILURE;
}
