#include <iostream>
#include <fstream>
#include <matvect.h>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " " << endl ;
}

int main(int argc, char **argv)
{
  
  Mat m1(4,4);
  m1.Identity();

  Mat m2=m1;
  m2.Zero();
  
  
  for(int i=0;i<4;i++)
    m2(0,i)=1;
  
  m2 -= (2.*m1);
  m1 *= 2;
  Mat m3 = m1*m2;
  cout << m1 << endl;
  cout << m2 << endl;
  cout << m3 << endl;
  
  return 0;
}

