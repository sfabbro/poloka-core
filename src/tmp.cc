#include <iostream>
#include <string>


template<class T>
class toto {
public:
  toto() {}
  ~toto() {}
  
  void write_members() {
    std::cout << "base class!" << std::endl;
  }
  
private:
  T i;
};



template<>
class toto<int> {
public:
  toto<int>() {}
  ~toto<int>() {}
  
  void write_members() {
    std::cout << "base class toto<int>!" << std::endl;
  }
  
private:
  int glop;
};



int main()
{
  toto<double> td;
  toto<int>    ti;
  float        a[100];
  
  
  td.write_members();
  ti.write_members();
  
  std::string str = "Hello!";
  std::cout << " str[0]=" << str[0] << std::endl;
}
