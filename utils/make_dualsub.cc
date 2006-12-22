#include<iostream>

#include "sub.h"
#include "dual_sub.h"
 
 
static void usage(const char * prog)
{
  cerr << " usage : " << endl;
  cerr << prog << " [-o] [-S subfile] [-p <dual_path>] [-l <list>]" << endl
       << "        -o : overwrite existing images (default = false)" << endl
       << "        -S subfile : to provide a subfile name (default = \"subfile\")" << endl
       << "        -p dual_path : to provide path to dual images (default = \"../../dbim/\")" << endl       
       << "        -l list : to provide list to match detection (default = \"../../dbim/fakesn.list\")" << endl;       
       
  exit(1);
}

int main(int argc, char **argv)
{
  string subfile = "subfile";
  string dual_path = "../../dbim/MC";
  string listname = "../../dbim/fakesn.list";
  bool overwrite = false;

  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-')  usage(argv[0]);
      switch (arg[1])
	{
	case 'h' : usage(argv[0]);
	case 'o' : overwrite = true; break;
	case 'S' : i++; if (i >= argc) usage(argv[0]); subfile = argv[i]; break;
	case 'p' : i++; if (i >= argc) usage(argv[0]); dual_path = argv[i]; break;		
	case 'l' : i++; if (i >= argc) usage(argv[0]); listname = argv[i]; break;	
	
	default : cerr << " do not understand " << arg << endl; usage(argv[0]);
	}
    }
  // if sub/matcheddet.list exists the subtraction is normally well finished
  // You don't want to rerun the whole process
  // except when explicitly required with overwrite
  if (!FileExists("sub/matcheddet.list") || overwrite)
    {
      Sub sub(subfile, overwrite);
      sub.DoIt();
 
      ProcessDualSub(sub,dual_path,listname);
    }
} 
