#include <hbook.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NWPAWC 100000
float pawc_[NWPAWC];
#define TOPDIR "//TOPDIR"


/* routines to handle toads ascii lists : */
/* shared with l2ncdf */
#include "list_utils.c"


static char* unslashed_topdir(char *Topdir)
{
if (strstr(Topdir,"//") == Topdir) return Topdir+2; else return Topdir;
}


int open_hbook_file(char *name)
{
int istat;
 int lrec = 1024; 
HROPEN(1, unslashed_topdir(TOPDIR),name, "N",lrec,istat);
if (istat == 0) return 1;
return 0;
}

int tuple_book(int Id, char *TupleName, char **Tags, int Dim)
{
#define MAX_LENGTH 30
char tags[MAXVAR][MAX_LENGTH]; /* emulate the implementation of fortran arrays of strings */
int i;
char title[128];
if (Dim > MAXVAR) { printf(" tuple %d truncated\n",Id); Dim = MAXVAR; }
for (i=0; i<Dim; i++) 
  {strncpy(tags[i],Tags[i],MAX_LENGTH-1);}
strcpy(title,TupleName); /* may be useless, but not sure */
HBOOKN(Id,title,Dim,TOPDIR,60000,tags);
return 1;
}

void tuple_end(int Id)
{int icycle = 0;
HROUT(Id,icycle," ");
HREND(unslashed_topdir(TOPDIR));
}


static void my_hfn(int *Id,float *X)
{
  HFN(*Id,X);
}


int make_tuple(const char *AsciiName, char* HbkFileName, int Id)
{
FILE *in;
char **tags;
int dim;
in = fopen(AsciiName,"r");
if (!in) 
  {
  printf(" cannot open %s\n",AsciiName);
  return 0;
  }
tags = decode_tags(in, &dim);
if (!tags) return 0;
if (!open_hbook_file(HbkFileName))
  {
    printf( " could not open %s\n",HbkFileName); fclose(in); return 0;
  }
tuple_book(Id, "toto", tags, dim);
read_data(in, dim, my_hfn, &Id);
tuple_end(Id);
return 1;
}


int main(int argc, char **argv)
{
char *hbkName=NULL;
 char string[256];
 if (argc == 2)
   {
     sprintf(string,"%s.hbk",CutExtension(argv[1]));
     hbkName = string;
   }
 if (argc == 3) hbkName = argv[2];
 if (!hbkName)
  {
  printf(" syntax : l2tup <ascii_list> <hbbokfile> \n");
  exit(1);
  } 
HLIMIT(NWPAWC);
make_tuple(argv[1], hbkName, 1);
return 0;
}
