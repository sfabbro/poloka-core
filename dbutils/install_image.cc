#include <iostream>

#include "fileutils.h"
#include "dbimage.h"
#include "fitsimage.h"

int process_one_file(const char *path, const char *file, const DbImageKind Kind)
{
FitsHeader header(file);
if (!header.IsValid())
  {
    cerr << " file " << file << " does not seem to be a fits file " << endl;
    return 0;
  }
if (string(header.KeyVal("TOADTYPE")) == "ZERO")
  {
    cerr << " the file : " << file << " seems to be a bias " << endl;
    return 0;
  }
// cout << " installing in " << path << " from " << file << endl;
return InstallImage(path, file, Kind);
}



int main(int argc, char **argv)
{
if (argc < 2)
  {
    cerr << "install_image [-cal] [-elixir] [-r] <directory> <fits image(s)> " << endl
	 << "  Creates a DbImage in <directory> from <fits_image(s)> " << endl
	 << "   -r  : creates <directory> path relatively to <fits_image(s)> " << endl
         << "   -cal: installs the <fits_image(s)> as calibrated (default is raw)" << endl
         << "   -elixir: installs the <fits_image(s)> as elixir (default is raw)" << endl;
                             
    exit(1);
  }

bool relative_path = false;
const char* path = NULL;
DbImageKind kind = Raw;
for (int i=1; i<argc ; ++i)
  {
    const char *arg = argv[i];
    
    if (strstr(arg,"-r") == arg)
      {
      relative_path = true;
      continue;
      }
    if (strstr(arg,"-cal") == arg)
      {
	kind = Calibrated; continue;
      }
    if (strstr(arg,"-elixir") == arg)
      {
	kind = Elixir; continue;
      }

    if (!path) 
      {
        path = arg;
        continue;
      }
  const char *file = arg;

  if (!FileExists(file))
    {
      cerr << " the file : " << file << " does not exists " << endl;
      continue;
    }

  // transforms the path to a standard one (i.e. explicit)
  // (cannot be done once for all when relative)

  string this_path = path;
  string p;
  if (relative_path)
    {
      p = StandardPath((DirName(file)+ "/" + this_path));
    }
  else p = StandardPath(this_path);
  if (p == "")
    {
      cerr << " problems to process " << this_path << endl;
      exit(-1);
    }
  this_path = p;
  if (!FileExists(this_path))
    {
      cerr << " The found path " << this_path << " does not exist (for image :" << file << ")" << endl;
      exit (-1);
    }
  process_one_file(this_path.c_str(), file, kind);
  }
return 0;
}
