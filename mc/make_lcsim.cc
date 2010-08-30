#include<iostream>

#include "lcsim.h"
#include "fileutils.h"
#include "myrdm.h"


static void usage(const char * prog)
{
  cerr << " usage : " << endl;
  cerr << prog << ' ' << " [-o] [-S simfile] [-L snialist] [-C randomize_code]" << endl
	<< "        -S simfile : to provide a simfile name (default = \"./simfile\")" << endl
	<< "        -L snialist  : to use an allready computed snialist " << endl
	<< "        -C code  : to randomize a part of parameters of an already computed snialist ; add each integer : " << endl
	<< "                  1  : H_Extinction" << endl
	<< "                  2  : Stretch" << endl
	<< "                  4  : Redshift" << endl
	<< "                  8  : MaxDay" << endl
	<< "                 16  : Dispersion" << endl
	<< "                 32  : Color" << endl
	<< "                 64  : Alpha" << endl
	<< "                128  : Beta" << endl
	<< "                256  : H_Rv" << endl
	<< "                512  : MW_Rv" << endl
	<< "                1024 : Positions" << endl				
	<< "            default  : 0 " << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  RandomSeed();
  
  string simfile = "simfile";
  string listname="";
  int code=2047;
  int new_code = 0;
  bool overwrite = false;
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-')  usage(argv[0]);
      switch (arg[1])
	{
	case 'h' : usage(argv[0]);
	case 'S' : i++; if (i >= argc) usage(argv[0]);
		simfile = argv[i]; break;
	case 'o' : overwrite = true;break;
	case 'L' :i++; if (i >= argc) usage(argv[0]);
		listname = argv[i]; break;
	case 'C' :i++; if (i >= argc) usage(argv[0]);
		new_code = atoi(argv[i]); break;
	default : cerr << " do not understand " << arg << endl; usage(argv[0]);
	}
    }
 if ( FileExists("snia.list") && !overwrite && listname!="snia.list" ) { cerr << "sim allready done. use -o option to overwrite " << endl; exit(1);}
  Sim simu(simfile);
  if(simu.Master() !="")simu.MakeGeoRef("georef");
  if (listname!="") {
  if(FileExists(listname)) {simu.ReadSnList(listname); code = new_code;}
  else {cerr << listname << "  not present" << endl; exit(1);}
  } 
  simu.MakeSnList(code);
  simu.SaveList("snia.list");
  simu.MakeSubfiles();
  simu.PrintSim();
  ofstream pr("done.txt");
  simu.PrintSim(pr);
  pr.close();
return (0);
 } 
