#include <cstdlib>
#include <string>
#include <iostream>

#include <poloka/fileutils.h>
#include <poloka/polokaconf.h>
#include <poloka/polokaexception.h>

static string ConfFileName;

void SetDatacardsFileName(const string &NewFileName)
{
  if (!FileExists(NewFileName))
    throw(PolokaException("SetDatacardsFileName : cannot find " + NewFileName));
  cout << " SetDatacardsFileName: setting the default '"
       << NewFileName << "'" << endl;
  ConfFileName = NewFileName;
}


string DefaultDatacards(const string& DefaultName)
{

  if (!ConfFileName.empty()) return ConfFileName;
  if (FileExists(DefaultName)) return DefaultName;
  char *where = getenv("POLOKA_CONF_DIR");
  if (where)
    return (AddSlash(where) + DefaultName);
  else // expect troubles
    return "";
}

    
