#ifndef FITSSUBS_H
#define FITSSUBS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Define the instances of the FORTRAN COMMON blocks that DAOphot II and 
   Allstar use. */
struct Common_Size
{
  int ncol;  /* The number of columns (ie. x axis range) in the image file. */
  int nrow;  /* The number of rows (ie. y axis type) in the image file. */
};

#define SIZE size_
extern struct Common_Size SIZE;
  
/*! Opens the image file without reading it.*/
#define ATTACH attach_
void ATTACH(const char* image, int* open);

/*! Reads pixels from the attached image (or the working copy of the image).
  If text is "DATA" then the pixels will be copied from the image, if
  text is "COPY" then the pixels will be copied from the working image copy. */
#define RDARAY rdaray_
void RDARAY(const char* text, int *lx, int *ly, int *mx, 
	    int* my, const int* maxx, float* func, int* ier);

/*! Writes pixels to the attached image (or the working copy of the image).
  If text is "DATA" then the pixels will be copied to the image, or if
  text is "COPY" then the pixels will be copied to the working image copy. */
#define WRARAY wraray_
void WRARAY(const char* text, int *lx, int *ly, int *mx, 
            int* my, const int* maxx, float* func, int* ier);

/*! Copy a fits image */
#define COPPIC coppic_
void COPPIC(const char* picture, int* ier);

/*! Closes an image file. */
#define CLPIC clpic_
void CLPIC(const char* text);

/*! Deletes the image whose name is specified as image. */
#define DELPIC delpic_
void DELPIC(const char* image, int* ier);

/*! Give the name of the current file */
#define LIST list_
void LIST(const char* file);

/*! Flush the working copy image */
#define SPLASH splash_
void SPLASH();

#ifdef __cplusplus
}
#endif


#endif /* FITSSUBS__H */
