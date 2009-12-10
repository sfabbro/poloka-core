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


int open_hbook_file(char *name, const int Lrec)
{
int istat;
 int lrec = Lrec;
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


int make_tuple(const char *AsciiName, char* HbkFileName, int Id, int Lrec)
{
FILE *in;
char **tags;
int dim;
int nrec = 0;
in = fopen(AsciiName,"r");
if (!in) 
  {
  printf(" cannot open %s\n",AsciiName);
  return 0;
  }
tags = decode_tags(in, &dim);
if (!tags) return 0;
 if (!open_hbook_file(HbkFileName, Lrec))
  {
    printf( " could not open %s\n",HbkFileName); fclose(in); return 0;
  }
tuple_book(Id, "toto", tags, dim);
 nrec = read_data(in, dim, my_hfn, &Id, /* print_bad_lines =  */ false);
 printf(" %d entries written to %s\n", nrec, HbkFileName);
tuple_end(Id);
return 1;
}

static void usage(const char* prog)
{
  printf(" syntax : \n %s <ascii_list> <hbbokfile> [-l lrec] \n",prog);
  exit(1);
} 



int main(int argc, char **argv)
{
  char *fileNames[2] = {NULL,NULL};
  int lrec = 1024;
  int count = 0;
  for (int i=1; i< argc; ++i)
    {
      char *arg= argv[i];
      if (arg[0] == '-' && arg[1] != '\0')
	{
	  switch (arg[1])
	    {
	    case 'l' : ++i; lrec = atoi(argv[i]); break;
	    case 'h' :
	    default:  usage(argv[0]);
	    }
	}
      else
	{
	  fileNames[count++] = arg;
	}
    }
  if (count == 0 || count > 2) usage(argv[0]);
  char *hbkName=NULL;
  char string[256];
  if (count == 1)
    {
      sprintf(string,"%s.hbk",CutExtension(fileNames[0]));
      hbkName = string;
    }
  else hbkName = fileNames[1];
  HLIMIT(NWPAWC);
  if (!make_tuple(fileNames[0], hbkName, 1, lrec)) return 1;
  return 0;
}
