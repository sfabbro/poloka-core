#include <iostream>

#include <poloka/dbimage.h>
#include <poloka/stringlist.h>
#include <poloka/fileutils.h>

static void usage(const char* progname)
{
  cerr << "Usage: " << progname << " [OPTION]... DBIMAGE|DBTAG...\n"
       << "Print the path a DBIMAGE file or all paths refered by a DBTAG\n\n"
       << "    -a    : also outputs non existent dbimages (as well as dangling links)\n"
       << "    -what : full path of what=raw|cal|flat|dead|sat|satgz|cat|cos|weight\n"
       << "    -dump : dumps the current $POLOKA_DB_CONFIG file content\n"
       << "    -example : provides a dbconfig file example\n";
  exit(EXIT_FAILURE);
}

static int PrintIfAbsent = 0;

int PrintOut(const string &FileName)
{
  if (PrintIfAbsent || FileExists(FileName))
    {
      cout << FileName; 
      return 1;
    }
  return 0;
}


int dump_info(const DbImageList &list, char **which_info, const int ninfo)
{

  int status = EXIT_SUCCESS;
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
	      if (fileName.empty())
		{
		  cerr << "unrecognised option file info: " << which_info[i] << endl;
		  status = EXIT_FAILURE;
		}
	      if (i>0) cout << ' ';
	      count += PrintOut(fileName);
	    }
	  if (count) cout << endl;
	}
    }
  return status;
}

int main(int argc,char **argv)
{
  if (argc <=1) usage(argv[0]);

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
		default : usage(argv[0]);
		}
	    }
	  else if (string(arg) == "-dump")
	    {
	      DbConfigDump();
	      exit(EXIT_SUCCESS);
	    }
	  else if (string(arg) == "-example")
	    {
	      DbConfigExample();
	      exit(EXIT_SUCCESS);
	    }
	  else which_info[n_info++] = arg+1; 
	  continue;
	}
      DbImage dbimage(arg);
      if (!dbimage.IsValid()) path_list.push_back(arg);
      else image_list.push_back(dbimage);
    }
  DbInit();
  if (path_list.empty() && image_list.empty()) 
    path_list.push_back(getenv("PWD"));
  int failures = 0;
  for (StringCIterator pi = path_list.begin(); pi != path_list.end(); ++pi)
    {
      DbImageList a_list((*pi).c_str());
      if (a_list.empty()) failures++;
      dump_info(a_list, which_info, n_info);
    }
  dump_info(image_list, which_info, n_info);
  if (failures != 0) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
