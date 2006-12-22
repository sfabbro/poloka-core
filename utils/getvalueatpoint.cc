#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <fstream.h>
#include <getopt.h>
#include <string.h>
#include <list>

#include "fitsimage.h"
#include "reducedimage.h"
#include "stringlist.h"
#include "gtransfo.h"
#include "wcsutils.h"

static void usage(const char * prog)

{
  cerr << " usage : " << endl;
  cerr << prog << ' ' << "[-R Ra Dec][-P x y]  (-D [-W][-C] <dbimage name list>) or (fits image list)" << endl
        << "        -R : Ra Dec coordinate" << endl
        << "        -P : Pixel coordinate on the first dbimage " << endl
        << "        -D : to use DBImage (check .dbconfig)" << endl
        << "        -W : Show value on Weight image for DBimage" << endl     
        << "        -C : Show value on Calibrate image for DBimage" << endl     ;
 
  exit(1);
}




int
main(int argc, char **argv)
{
double x,y;
bool pix = true;
bool cal = false;
bool weight = false;
bool coor = false;
bool isdbim = false;
StringList  image_list;
if (argc < 5) usage(argv[0]);
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] == '-')  
      switch (arg[1])
        {
        case 'W': weight = true; break;
	case 'C': cal= true  ; break;
	case 'D': isdbim = true ; break;
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
      else image_list.push_back(argv[i]);
    }
if (!coor)   usage(argv[0]);  
if ( cal == false && weight == false) cal = true;
double Ra,Dec;
if(!pix) { Ra = x ; Dec = y; cout << "Ra = "<< Ra <<"   Dec = "<< Dec << endl;}

for (StringIterator i = image_list.begin() ; i != image_list.end(); )
 {
if (isdbim)
{
ReducedImage current(*i);
if(cal)
 {

	FitsImage current_fits (current.FitsName());
 	Frame current_Frame(current_fits, WholeSizeFrame);
 	Gtransfo *Pix2RaDec;
 	if (!WCSFromHeader(/*dynamic_cast <FitsHeader> */current_fits, Pix2RaDec))
 	{
 		cout << "cannot handle "<<current.FitsName() <<" without a WCS " << endl;
 		continue;
 	}
 	Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, current_Frame);

if(pix) { Pix2RaDec->apply(x,y,Ra,Dec);cout << "Ra = "<< Ra <<"   Dec = "<< Dec << endl;}
else RaDec2Pix->apply(Ra,Dec,x,y);
pix = false;
if (current_Frame.InFrame(x,y)) cout << current.FitsName() <<"    \tx = "<<x <<"    \ty = " << y <<  "   \tvalue = " << current_fits(x,y)<< endl;
else cout << "point ouside " <<  current.FitsName() << endl; 
}
if(weight)
 {

	FitsImage current_fits (current.FitsWeightName());
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
if (current_Frame.InFrame(x,y))cout << current.FitsWeightName() <<"    \tx = "<<x <<"    \ty = " << y <<  "   \tvalue = "  << current_fits(x,y)<< endl;
else cout << "point ouside " << current.FitsWeightName() << endl; 
}
}
else 
{
	FitsImage current_fits (*i);
 	Frame current_Frame(current_fits, WholeSizeFrame);
 	Gtransfo *Pix2RaDec;
 	if (!WCSFromHeader(/*dynamic_cast <FitsHeader> */current_fits, Pix2RaDec))
 	{
 		cout << "cannot handle "<< *i <<" without a WCS " << endl;
 		continue;
 	}
 	Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, current_Frame);

if(pix){ Pix2RaDec->apply(x,y,Ra,Dec);cout << "Ra = "<< Ra <<"   Dec = "<< Dec << endl;}
else RaDec2Pix->apply(Ra,Dec,x,y);
pix = false;
if (current_Frame.InFrame(x,y))cout << *i <<"    \tx = "<<x <<"    \ty = " << y <<  "   \tvalue = " << current_fits(x,y)<< endl;
else cout << *i  << " point ouside" << endl; 
}	 
 ++i;
}
}
