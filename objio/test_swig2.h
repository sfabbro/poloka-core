// -*- C++ -*-
// 
// 
// 
// 
#ifndef TEST_SWIG2_H
#define TEST_SWIG2_H


#define TOTO int


class test_1 {
public:
  test_1();
  ~test_1();

private:
  mutable int2 a1;
  int2* const a2;
  int4**  a3;
  float8**const* a4;
  const float8* a5;
  float4 const& a6;
  float4 toto[100][22];
  
  list<double> toto;
  list<double>& rtoto;
};


#endif

