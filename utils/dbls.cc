#include <iostream>

#include "dbimage.h"
#include "stringlist.h"
#include "fileutils.h"


void  usage()
{
  cerr <<"dbls [<options>] [<symbolic_path> ... ] [<image> ... ] \n" 
       <<" options among -raw -cal -flat -dead -sat -satgz -cat -cos -weight"
       << endl
       <<"  -a also outputs non existent file names(as well as dangling links)"
       << endl
       << "'dbls -dump' reads and dumps the config file contents" 
       << endl
       << "'dbls -example' provides a dbconfigfile example" 
       << endl;
}

static int PrintIfAbsent = 0;

int PrintOut(const string &FileName)
{
  if (PrintIfAbsent || FileExists(FileName.c_str()))
    {
      cout << ' ' << FileName; return 1;
    }
  return 0;
}


void dump_info(const DbImageList &list, char **which_info, const int ninfo)
{
  for (DbImageCIterator dbi = list.begin(); dbi != list.end(); ++dbi)
    {
      const DbImage& image = *dbi;
      if (ninfo == 0) image.dump();
      else
	{   
	  int count = 0;
	  for (int i=0; i<ninfo; ++i)
	    {
	      string fileName = image.GetFileName(which_info[i]);
	      if (fileName.length() == 0)
		{
		  cerr << " unrecognised option -"<< which_info[i] << endl;
		  usage();
		  exit(-1);
		}
	      count += PrintOut(fileName);
	    }
	  if (count) cout << endl;
	}
    }
}

int main(int argc,char **argv)
{
  if (argc <=1)
    {
      usage(); exit(1);
    }
  StringList path_list;
  DbImageList image_list;
  DbConfigSetDumpLevel(1); /* short print out */
  char **which_info = (char **) calloc(argc,sizeof(char*));
  int n_info = 0;

  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] == '-')  
	{ 
	  if (arg[2] == '\0')
	    {
	      switch (arg[1])
		{
		case 'a' : PrintIfAbsent = 1; break;
		default : usage(); exit(1);
		}
	    }
	  else if (string(arg) == "-dump")
	    {
	      DbConfigDump();
	      exit(0);
	    }
	  else if (string(arg) == "-example")
	    {
	      DbConfigExample();
	      exit(0);
	    }
	  else which_info[n_info++] = arg+1; 
	  continue;
	}
      DbImage dbimage(arg);
      if (!dbimage.IsValid()) path_list.push_back(arg);
      else image_list.push_back(dbimage);
    }
  DbInit();
  if (path_list.size() == 0 && image_list.size() == 0) 
    path_list.push_back(getenv("PWD"));
  int failures = 0;
  for (StringCIterator pi = path_list.begin(); pi != path_list.end(); ++pi)
    {
      DbImageList a_list((*pi).c_str());
      if (a_list.size() == 0) failures++;
      dump_info(a_list, which_info, n_info);
    }
  dump_info(image_list, which_info, n_info);
  if (failures != 0) return -1;
  return 0;
}
