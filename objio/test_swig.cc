// -*- C++ -*-
// 
// 
// 
// 
#include <math.h>

#include <list>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "persistence.h"
#include "typemgr.h"
#include "test_swig.h"
#include "test_swig_dict.h"

#include "B__persister.h"
#include "BB__persister.h"


struct Toto {
  double x;
  double y;
  double z;
  int    id;
};


int main()
{
  cout << typemgr::size() << endl;
  
  
  int i,j,k;
  AA a;
  Point p;
  Star s;
  std::list<Star> ls;
  std::vector<double> double_vector;
  std::vector<std::string> string_vector;
  map<std::string,short> string_short_map;
  std::list< std::vector< std::map<string,unsigned int> > > stupidly_complex_example;
  B<int> b;
  BB<string,double> bb;
  
  p.x()=2;
  p.y()=5.243;
  
  s.id()=7;
  s.x()=55.555;
  s.y()=55.555;
  s.flux()=1234567;
  
  for(i=0;i<20;i++) {
    Star s_tmp;
    s_tmp.id()=i;
    s_tmp.x() = i;
    s_tmp.y() = i;
    s_tmp.flux()=i*i;
    ls.push_back(s_tmp);
  }
  
  for(i=0;i<10;i++) {
    std::stringstream sstrm;
    sstrm << "toto = " << i;
    string_vector.push_back(sstrm.str());
  }
  
  for(i=0;i<10;i++) 
    double_vector.push_back(sin(i));
  
  string_short_map["toto"]=2;
  string_short_map["titi"]=-2383;
  string_short_map["tutu"]=3435;
  string_short_map["tata"]=5;  
  
  bb.lt_.push_back("glop");
  bb.mtu_["tutu"]=3.14;
  
  obj_output<xmlostream> oo("ttt.xml");
  oo << a;
  oo << p;
  oo << s;
  oo << ls;
  oo << string_vector;
  oo << double_vector;
  oo << string_short_map;
  oo << b;
  oo << bb;
  
  //  persister<Star>* tps = (persister<Star>*)typemgr::getPersister(typeid(&s).name());
  persister_base* tps = typemgr::getPersister(typeid(&s).name());
  std::cout << " name=" << tps->name() << " version=" << tps->version() << endl;
  oo.write(tps);
  
  oo.close();

  exit(0);
  
  AA a2;
  Point p2;
  Star s2;
  std::list<Star> ls2;
  std::vector<std::string> string_vector2;
  std::vector<double> double_vector2;
  map<std::string,short> string_short_map2;
  
  
  std::cout << "about to read object..." << std::endl;
  obj_input<xmlistream> oi("ttt.xml");
  
  oi >> a2;
  oi >> p2;
  p2.print();
  oi >> s2;
  s2.print();
  oi >> ls2;
  std::cout << " read ls2: size=" << ls2.size() << std::endl;
  std::list<Star>::iterator it;
  for(it=ls2.begin();it!=ls2.end();it++)
    it->print();
  
  oi >> string_vector2;
  std::cout << " string_vector.size()=" << string_vector2.size() << std::endl;
  std::vector<std::string>::iterator it2;
  for(it2=string_vector2.begin();it2!=string_vector2.end();it2++)
    std::cout << *it2 << std::endl;
  
  oi >> double_vector2;
  std::cout << " double_vector.size()=" << double_vector2.size() << std::endl;
  std::vector<double>::iterator dit2;
  for(dit2=double_vector2.begin();dit2!=double_vector2.end();dit2++)
    std::cout << *dit2 << std::endl;
  
  oi >> string_short_map2;
  std::cout << "string_short_map2.size()=" << string_short_map2.size() << std::endl;
  std::map<std::string,short>::iterator ssm2it;
  for(ssm2it=string_short_map2.begin();ssm2it!=string_short_map2.end();ssm2it++)
    std::cout << ssm2it->first << " " << ssm2it->second << std::endl;
}




  
//  Toto toto;
//  for(i=0;i<10;i++)
//    oo << toto;
//  
//  dict_output<xmlostream> dop("ddd.xml");
//  persister<Star> ps(s);
//  dop.write(ps);
