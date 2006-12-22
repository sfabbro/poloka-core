#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <getopt.h>
#include <string.h>
#include <list>

#include "fitsimage.h"
#include "stringlist.h"
#include "gtransfo.h"
#include "wcsutils.h"
#include "reducedimage.h"
static void usage(const char * prog)

{
  cerr << " usage : " << endl;
  cerr << prog << ' ' << "[-R Ra Dec][-P x y]  [-W][-C] <dbimage name list>)" << endl
        << "        -R : Ra Dec coordinate" << endl
        << "        -P : Pixel coordinate on the first dbimage " << endl    ;
 
  exit(1);
}


static string first_word(const char *line, const char sep = ' ')
{
  int start;
  int end_line = strlen(line);

  while ( line[start] == sep && start < end_line) start++;
  int end = start;
  while ( line[end] != sep &&  end < end_line) end++;
  return string(line).substr(start,end);
}




int
main(int argc, char **argv)
{
string FileName;
double x,y;
bool pix = true;
bool coor = false;
StringList  image_list;
if (argc != 5) usage(argv[0]);
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] == '-')  
      switch (arg[1])
        {
	case 'R' :
	pix =false;
	case 'P' : 
	i =i + 2;
	if (i > argc) usage(argv[0]);
	x = atof(argv[i-1]);
	y = atof(argv[i]);
	coor= true;	
	break;
        default : cerr << " do not understand " << arg << endl; usage(argv[0]);
        }
      else FileName=argv[i];
    }
    
if (!coor)   usage(argv[0]);  



	FILE *file = fopen(FileName.c_str(),"r");
	if (!file)
	{
		cerr << " cannot open \"" << FileName << "\"" << endl;
		exit(1);
	}
	char line[512];

while (fgets(line,512,file))
{
	if (line[0] == '#') continue;
	if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
		


 	char *start_line = line+strspn(line," "); // skip spaces
 	if (start_line[strlen(start_line)-1] == '\n' )
		{
		 start_line[strlen(start_line)-1] = '\0' ;
		}
 	if (strlen(start_line) == 0) continue; // skip blank lines
 
 	string currentName = first_word(start_line);
 	RemovePattern(currentName," ");
 image_list.push_back(currentName);

}	


 fclose(file);


double i_count = 0.0;
double r_count = 0.0;
double g_count = 0.0;
double z_count = 0.0;

double Ra,Dec;
if(!pix) { Ra = x ; Dec = y; cout << "Ra = "<< Ra <<"   Dec = "<< Dec << endl;}






for (StringIterator i = image_list.begin() ; i != image_list.end(); i++ )
 {

ReducedImage current(*i);
	FitsImage current_fits (current.FitsWeightName());
	string band=current_fits.KeyVal("FILTER");

        
        if(band[0]> 64 && band[0]< 91) band[0]=band[0]+32;
        band = (char) band[0];

 	Frame current_Frame(current_fits, WholeSizeFrame);
 	Gtransfo *Pix2RaDec;
 	if (!WCSFromHeader(/*dynamic_cast <FitsHeader> */current_fits, Pix2RaDec))
 	{
 		cout << "cannot handle "<< current.FitsWeightName() <<" without a WCS " << endl;
 		continue;
 	}
 	Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, current_Frame);

if(pix){ Pix2RaDec->apply(x,y,Ra,Dec);cout << "Ra = "<< Ra <<"   Dec = "<< Dec << endl;}
else RaDec2Pix->apply(Ra,Dec,x,y);
pix = false;
if (current_Frame.InFrame(x,y))
{
double pix_val = current_fits(x,y);
cout << current.FitsWeightName() <<"    \tx = "<<x <<"    \ty = " << y <<  "   \tvalue = "  << pix_val << endl;
switch(band[0])
{
case 'i' : i_count+=pix_val ; break;
case 'r' : r_count+=pix_val ; break;
case 'g' : g_count+=pix_val ; break;
case 'z' : z_count+=pix_val ; break;
default : cout << "unbleme" << endl;
}

}
else cout << "point ouside " << current.FitsWeightName() << endl; 


	 

}

cout << " ival :" << i_count << endl;
cout << " rval :" << r_count << endl;
cout << " gval :" << g_count << endl;
cout << " zval :" << z_count << endl;

}
