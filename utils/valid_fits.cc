#include <iostream>
#include "fileutils.h"


int main(int argc, char **args)
{
for (int i=1; i< argc; ++i)
  {
    if (IsFits(args[i])) cout << args[i] << endl;
  }
return 1;
}
