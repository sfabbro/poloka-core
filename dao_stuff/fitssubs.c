#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>
#include "cdaophot.h"

/******************************************************************************/

static void full_path_image_name(char *new_name, const char *name)
{
  char *c;
  strcpy(new_name, getenv("PWD"));
  strcat(new_name,"/");
  strcat(new_name, name);
  if (strrchr(name,'.')) {c = strrchr(new_name,'.'); *c ='\0';}
  strcat(new_name,".fits");
}

static void switch_file_extension(char *new_name, const char *name, 
				  const char *ext)
{
  char *dot;
  const char *c = ext;
  char cname[NCHARFILE];
  size_t slen, i;

  /* Remove any whitespace around the filename. */
  sscanf(name, "%s", cname);
  
  /* Replace extension. If not found, append it. */
  dot = strrchr(cname, '.');
  if (!dot) dot = strrchr(cname, '\0');
  *dot = '.'; ++dot;
  while (*c != '\0') 
    {
      *dot = *c;
      ++dot; ++c;
    }
  *dot = '\0';

  /* Convert into a fortran string */
  strncpy(new_name, cname, NCHARFILE);
  slen = strlen(new_name);
  for(i=slen;i<NCHARFILE;i++) new_name[i]=' ';  
}

extern struct Common_Filename FILNAM;
extern struct Common_Size SIZE;

/* pointers to fits file structure used by cfitsio */
fitsfile *imdata = 0;
fitsfile *imcopy = 0;
int bitpix = 0;


/******************************************************************************/
void ATTACH(const char* image, int* open)
{
  char full_image_name[NCHARFILE], short_name[NCHARFILE];
  int iostatus = 0;
  long naxes[2];

  /* Remove any whitespace around the filename. */
  sscanf(image, "%s", short_name);

  if (*open) CLPIC("DATA");
  full_path_image_name(full_image_name, short_name);

  /* Open file */
  fits_open_file(&imdata, full_image_name, READONLY, &iostatus);
  if (iostatus != 0) 
    {
      fits_report_error(stderr, iostatus);
      return;
    }

  /*
  printf(" fname ='%s'\n",fname);
  printf(" full_image_name ='%s'\n",full_image_name);
  */

  *open = 1;

  /* Get size and image type */
  fits_get_img_param(imdata, 2, &bitpix, NULL, naxes, &iostatus);
  if ((iostatus != 0) 
      && (bitpix != SHORT_IMG)
      && (bitpix != LONG_IMG)
      && (bitpix != FLOAT_IMG)
      && (naxes[0] < 4 )
      && (naxes[1] < 4 ))
    {
      fits_report_error(stderr, iostatus);
      fprintf(stderr, " ATTACH : Unable to get proper size or type \n");
      CLPIC("DATA");
      return;
    }
  
  /* Assign common block sizes and type */
  SIZE.ncol = (int) naxes[0];
  SIZE.nrow = (int) naxes[1];

  printf("\n%38s%5i%5i\n", "Picture size:", SIZE.ncol, SIZE.nrow);      
  printf("%38s%5i\n\n", "Picture bit/pixel:", bitpix);      

  switch_file_extension(FILNAM.COOFIL, short_name, "coo");
  switch_file_extension(FILNAM.MAGFIL, short_name, "ap");
  switch_file_extension(FILNAM.PSFFIL, short_name, "psf");
  switch_file_extension(FILNAM.PROFIL, short_name, "als");
  switch_file_extension(FILNAM.GRPFIL, short_name, "grp");
  /*
  printf(" coofil = '%s'\n", FILNAM.COOFIL);
  printf(" magfil = '%s'\n", FILNAM.MAGFIL);
  printf(" psffil = '%s'\n", FILNAM.PSFFIL);
  printf(" profil = '%s'\n", FILNAM.PROFIL);
  printf(" grpfil = '%s'\n", FILNAM.GRPFIL);
  */
}
/******************************************************************************/

/******************************************************************************/
void RDARAY(const char* text, int *lx, int *ly, int *mx, 
	    int* my, const int* maxx, float* func, int* ier)
{
  int j, nx, anynul = 0;
  long jy;
  fitsfile *fits_file;
  
  /* Get fits file to extract data */
  if (strncmp(text,"DATA",4) == 0)  fits_file = imdata;
  else fits_file = imcopy;
  
  /* Convert upper right corner */
  *mx = *lx + *mx - 1;
  *my = *ly + *my - 1;

  /* Check bounds with image borders */
  *lx = (1 < *lx) ? *lx : 1;
  *ly = (1 < *ly) ? *ly : 1;
  *mx = (SIZE.ncol < *mx) ? SIZE.ncol : *mx;
  *my = (SIZE.nrow < *my) ? SIZE.nrow : *my;

  *my = *my - *ly + 1;
  nx = *mx-*lx+1;
  *ier = 0;
  for (j=0; j<*my; ++j)
    {
      jy = *lx + SIZE.ncol * (*ly-1 + j);
      fits_read_img_flt(fits_file, 1, jy, nx, 1.e38, 
			func + j*SIZE.ncol, &anynul, ier);
      /* printf(" reading elts %d to %d\n", jy, jy+nx-1); */
    }  
  *mx = nx;
  if (*ier != 0) fits_report_error(stderr, *ier);
}
/******************************************************************************/


/******************************************************************************/
void WRARAY(const char* text, int *lx, int *ly, int *mx, 
	    int* my, const int* maxx, float* func, int* ier)
{
  int j,jy,nx,ny;
  fitsfile *fits_file;
  
  /* Get fits file to extract data */
  if (strncmp(text,"DATA",4) == 0)  fits_file = imdata;
  else fits_file = imcopy;
  
  /* Convert upper right corner */
  *mx = *lx + *mx - 1;
  *my = *ly + *my - 1;

  /* Check bounds with image borders */
  *lx = (1 < *lx) ? *lx : 1;
  *ly = (1 < *ly) ? *ly : 1;
  *mx = (SIZE.ncol < *mx) ? SIZE.ncol : *mx;
  *my = (SIZE.nrow < *my) ? SIZE.nrow : *my;
  
  nx = *mx - *lx + 1;
  ny = *my - *ly + 1;
  *ier = 0;
  for (j=0; j<ny; ++j)
    {
      jy = *lx + SIZE.ncol * (*ly-1 + j);
      fits_write_img_flt(fits_file, 1, jy, nx, func + j*SIZE.ncol, ier);
      /* printf(" writing elts %d to %d\n", jy,jy+nx-1); */
    }
  *mx = nx;
  *my = ny;
  
  if (*ier != 0) fits_report_error(stderr, *ier);
}
/******************************************************************************/


/******************************************************************************/
void COPPIC(const char* picture, int* ier)
{
  char fitspic[512], fname[256];
  long naxes[2];

  *ier = 0;

 /* Remove any whitespace around the filename. */
  sscanf(picture, "%s", fname);

  full_path_image_name(fitspic, fname);

  /* Delete the picture if it already exists and if not, open it*/
  DELPIC(fitspic, ier);
  if (imcopy) fits_delete_file(imcopy, ier);
  else fits_create_file(&imcopy, fitspic, ier);

  if (*ier != 0) fits_report_error(stderr, *ier);

  /* Copy header, previous, current, and following hdu from imdata to imcopy */
  fits_copy_file(imdata, imcopy, 1, 1, 1, ier);
  naxes[0] = SIZE.ncol;
  naxes[1] = SIZE.nrow;

  /* Force bitpix to be float
     (because too lasy to adapt scale and zero in case of 16bits)*/

  fits_resize_img(imcopy, FLOAT_IMG, 2, naxes, ier);

   if (*ier != 0) fits_report_error(stderr, *ier);
}
/******************************************************************************/


/******************************************************************************/
void CLPIC(const char* text)
{
  int iostatus = 0;
  fitsfile *fits_file;
  
  /* Get fits file to close */
  if (strncmp(text,"DATA",4) == 0)  fits_file = imdata;
  else fits_file = imcopy;

  /* Close file */
  fits_close_file(fits_file, &iostatus);
  fits_report_error(stderr, iostatus);

  /* Cancel pointer */
  if (fits_file == imdata) imdata = 0;
  else imcopy = 0;
}
/******************************************************************************/


/******************************************************************************/
void DELPIC(const char* image, int* ier)
{
  /* Don't use fits_delete_file because it assumes the file is still open
   and what we want is just a plain remove */
  if (access(image, F_OK) == 0)
    {
      *ier = remove(image);
      if (*ier != 0) 
	fprintf(stderr,"  DELPIC: Unable to remove file %s \n", image); 
    }
}
/******************************************************************************/

/******************************************************************************/
void LIST(const char* file)
{
  printf("\n  Image file = %s\n\n", file);
}
/******************************************************************************/

/******************************************************************************/
void SPLASH()
{
  int iostatus = 0;
  fits_flush_file(imcopy, &iostatus);
  if (iostatus != 0)  fits_report_error(stderr, iostatus);
}
/******************************************************************************/
