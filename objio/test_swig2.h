// -*- C++ -*-
// 
// 
// 
// 
#ifndef TEST_SWIG2_H
#define TEST_SWIG2_H


#include <stdint.h>

#define TOTO int


class test_1 {
public:
  test_1();
  ~test_1();

private:
  mutable int16_t a1;
  int16_t* const a2;
  int32_t**  a3;
  double**const* a4;
  const double* a5;
  float const& a6;
  float toto[100][22];
  
  list<double> toto;
  list<double>& rtoto;
};


#endif

