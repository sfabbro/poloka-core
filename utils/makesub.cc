#include<iostream>

#include "sub.h"
#ifdef ENABLE_MC
#include "mcsub.h"
#endif 
 
static void usage(const char * prog)
{
  cerr << " usage : " << endl;
  cerr << prog << " [-o] [-S subfile] [-M]" << endl
       << "        -o : overwrite existing images (default = false)" << endl
       << "        -S subfile : to provide a subfile name (default = \"subfile\")" 
#ifdef ENABLE_MC
       << "        -M : MC mode for efficiency computation" << endl;
#else
       << "        -M : MC mode for efficiency computation (disabled at compilation)" << endl;
#endif
  exit(1);
}

int main(int argc, char **argv)
{
  string subfile = "subfile";
  bool overwrite = false;
#ifdef ENABLE_MC
  bool MCmode = false ;
#endif
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-')  usage(argv[0]);
      switch (arg[1])
	{
	case 'h' : usage(argv[0]);
	case 'o' : overwrite = true;break;
#ifdef ENABLE_MC
	case 'M' : MCmode = true;break;
#endif
	case 'S' : i++; if (i >= argc) usage(argv[0]);
	  subfile = argv[i]; break;
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
#ifdef ENABLE_MC
      if (MCmode)
	MCProcess(sub);
#endif
    }
} 
