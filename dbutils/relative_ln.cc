#include <iostream>

#include "fileutils.h"

/* This program sets symbolic links just like ln does.
The only difference is that the link value is as
relative as possible to the link itself. So that
if the link and target common root is moved, the
link still points correctly.
*/

void usage(const char *pg_name, char *message=NULL)
{
  if (message) cerr << message << endl;
  cerr << pg_name << " <actual file name> <link name> " << endl;
  cerr << "\t (sets a symbolic link, and only overwrites existing symbolic links)" << endl;
}

int main(int nargs, char **args)
{

if (nargs < 3)
  usage(args[0]);
string file = args[1];
string link = args[2];
#ifdef ISLINK_OK
if (FileExists(link) && !IsLink(link))
#endif
if (FileExists(link))
  {
    cerr << args[0] << " only overwrites links, not files " << endl;
    return -1;
  }
string link_value = RelativeLink(file.c_str(), link.c_str());
string command = "ln -fs " + link_value + " " + link;
 cout << " debug " << command << endl;
if (system(command.c_str()) != 0) return -1;
return EXIT_SUCCESS;
}
