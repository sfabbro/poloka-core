#include <string.h>
#include <iostream>
#include "fitsimage.h"


int main(int argc, char**argv)
{
  if (argc >=4 )
    {
      char *in_name = argv[1];
      char *key = argv[2];
      const char *comment = "Changed key";
      FitsHeader toChange(in_name,RW);
      if (argc==5) comment = argv[4];
      cout << "Changing key " << key << " from " << toChange.KeyVal(key) << " to " << argv[3] << endl;
      toChange.AddOrModKey(key,argv[3]);
    }
  else 
    {
      cout << "usage : changekey <FITS file> <Key> <New value> [comment]" << endl;
    }    
}


