#include <iostream>
#include <lightcurvepoint.h>

LightCurvePoint::LightCurvePoint() {
  julianday=0;
  flux=0;
  eflux=0;
  mag=99;
  emag_minus=0;
  emag_plus=0;
}
