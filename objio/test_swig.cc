// -*- C++ -*-
// 
// 
// 
// 
#include <list>

#include "test_swig.h"
#include "xmlstream.h"

#include "test_swig_dict.h"
#include "persister.h"

#include "dictio.h"


struct Toto {
  double x;
  double y;
  double z;
  int    id;
};


int main()
{
  obj_output<xmlostream> oo("ttt.xml");
  
  AA a;
  oo << a;
  
  Point p;
  oo << p;
  
  Star s;
  oo << s;
  
  std::list<Star> ls;
  int i;
  for(i=0;i<2;i++) {
    Star s_tmp;
    ls.push_back(s_tmp);
  }
  oo.write(ls);
  
  Toto toto;
  for(i=0;i<10;i++)
    oo << toto;
  
  dict_output<xmlostream> dop("ddd.xml");
  persister<Star> ps(s);
  dop.write(ps);
}
