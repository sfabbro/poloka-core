#include <iostream>
#include <stream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <string>

#include <time.h>

//#include "stlist.h"
//#include "starse.h"
//#include "poly.h"
//#include "clumproc.h"
#include "hbook.h"

#include "fileutils.h"
#include "fitsimage.h"
#include "dbimage.h"

#define NWPAWC 100000
float pawc_[NWPAWC];

static void HistoFileOpen(const string &name)
{
int istat=0;
HLIMIT(NWPAWC);
#define HBK_FILE_NAME "histos.hbk"

char dirname1[256];
DirName(name.c_str(),dirname1);
string histodir = AddSlash(string(dirname1));

string link = "/tmp/"+BaseName(dirname1);
string command = "ln -sf " + string(dirname1) + " " + link;
system(command.c_str());
string histoname = link +"/"+"histo.hbk";

char toto1[50]="TOPDIR";
char toto2[256];
sprintf(toto2,histoname.c_str());
char toto3[50]="N";

HROPEN(1,toto1,toto2,toto3,1024,istat);
if (istat !=0)
  {
    cout << "can not open file " << toto2 << endl;
    exit(2);
  }
command = "rm " + link;
system(command.c_str());

}


static void HistoFileClose()
{
int icycle=0;
char toto1[50]="TOPDIR";
char toto4[50]=" ";
HROUT(0,icycle,toto4);
HREND(toto1);
}


int image_histos(char *FileName)
{
if (!FileName) return 0;

DbImage dbimage(FileName);

if (!dbimage.IsValid())
  {
    cerr << " mkskyhist : can not find "<< FileName << endl;
    return 0;
  }

string name = dbimage.FitsImageName(Calibrated);

FitsImage img(name);
if (!img.IsValid()) return 0;

int nx = img.Nx();
int ny = img.Ny();
int ntot = nx*ny;

double average = 0;
double variance = 0;
double value;

for (int j=0; j<ny; j++)
for (int i=0; i<nx; i++) 
  {
    value = img(i,j);
    average += value;
    variance += value*value;
  }
average /= ntot;
variance = (variance/(ntot-1)) - (average*average);
double sigma = sqrt(variance);
cout << "image "<< FileName  << " ave " << average 
     << " sqrt(ave) " << sqrt(average) << " sigma " 
     << sigma << " count " << ntot << endl;
if (sigma == 0) return 0;

for(int loop=0; loop<3; loop++)
  {
    int count=0;
    double mean = average;
    average = 0.0;
    variance = 0.0;
    for (int j=0; j<ny; j++) for (int i=0; i<nx; i++) 
      {
	value = img(i,j);
	if (fabs(value-mean)<3*sigma)
	  {
	    average += value;
	    variance += value*value;
	    count++;
	  }
      }
    average /= count;
    variance = (variance/(count-1)) - (average*average);
    sigma = sqrt(variance);
    cout << "image "<< FileName  << " ave " << average << " sqrt(ave) " << sqrt(average) << " sigma " << sigma << " count " << count << endl;
    if (sigma == 0) return 0;
  }


HistoFileOpen(name.c_str());
static int id =10;
HBOOK1(id,FileName,50, average-5.*sigma, average+5.*sigma, 0);
for (int i=0; i<nx; i++) for (int j=0; j<ny; j++)
  {
  HFILL(id,img(i,j),0.,1.);
  }
HistoFileClose();

return 1;
}

int main(int nargs,char **args)
{

for (int i=1 ; i<nargs; i++)
  {
    cout << " processing " << args[i] << endl;
    image_histos(args[i]);
  }

return 1;
}
