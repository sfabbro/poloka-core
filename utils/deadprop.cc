#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <fstream.h>
#include <getopt.h>
#include <string.h>

#include "fitsimage.h"




/*
//NAXIS1  =                 2279 / length of data axis 1                          
//NAXIS2  =                 5385

        FitsHeader header(current.FitsName());
        if (header.HasActualKey("MJDATE"))
        {
                date = header.KeyVal("MJDATE");*/

int
main(int argc, char **argv)
{

if (argc < 2) exit(1);
for ( int i = 1; i < argc; i++)
{
string name = argv[i];
 cout << name << "\t\t";
FitsImage current_fits (name);
int x_axe = current_fits.KeyVal("NAXIS1");
int y_axe = current_fits.KeyVal("NAXIS2");
int count=0;
for(int x = 0 ; x <=x_axe;x++)
for(int y = 0 ; y <=y_axe;y++)

 if (current_fits(x,y) <= 0.5 ) count++;

 cout << count << " / " <<x_axe*y_axe << "\t\t";
 cout << double(count)/double(x_axe*y_axe)*100.0 << "%" << endl;

}
}
