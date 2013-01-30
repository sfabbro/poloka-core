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
      if (argc==5) comment = argv[4]; // ?
      string type;
      char * endptr;
      double dval = strtod(argv[3],&endptr);
      string oldval(toChange.HasKey(key) ? string(toChange.KeyVal(key)) :  "<was absent>");
      if (endptr != argv[3]) // means successful conversion to double
	{
	  type = "double";
	  toChange.AddOrModKey(key,dval);
	}
      else if (strcmp(argv[3],"T") == 0)
	{
	  type = "bool ";
          toChange.AddOrModKey(key,true);
	}
      else if (strcmp(argv[3],"F") == 0)
	{
	  type = "bool";
          toChange.AddOrModKey(key,false);
	}
      else 
	{
	  type = "string";
          toChange.AddOrModKey(key,string(argv[3]));
	}
      cout << "Changing key " << key << " from " << oldval
	   << " to " << argv[3] << " (type=" << type << ')' <<  endl;

    }
  else 
    {
      cout << "usage : changekey <FITS file> <Key> <New value> [comment]" << endl;
    }    
}


