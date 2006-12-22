#include "lambert_dust.h"


static void usage(const char * prog)
{
  cerr << " usage : " << endl;
  cerr << prog << ' ' << " [-R Ra Dec] [-G Long Lat] " <<  endl;
 
  exit(1);
}


int main(int argc, char **argv)
{
bool Ra=false;
bool Gal=false;
double x,y;
int i;
  for ( i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-')  usage(argv[0]);
      switch (arg[1])
        {

        case 'G' :	Gal =true;
			i++; if (i >= argc) usage(argv[0]);x = atof(argv[i]);
             		i++; if (i >= argc) usage(argv[0]);y = atof(argv[i]);
			break;
        case 'R' :	Ra =true;
			i++; if (i >= argc) usage(argv[0]);x = atof(argv[i]);
             		i++; if (i >= argc) usage(argv[0]);y = atof(argv[i]);
			break;
        default : cerr << " do not understand " << arg << endl; usage(argv[0]);
        }
    }
    
if( (Ra&&Gal)|| (!Ra && !Gal)) {cerr << " do not understand "  << endl; usage(argv[0]);}
Lambert_Dust my_dust;


if (Ra) cout  <<my_dust.Value_From_RaDec(x,y)<< endl;
else 	cout  <<my_dust.Value_From_Gal(x,y)  << endl;
return 0;
};
