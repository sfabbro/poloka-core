// -*- C++ -*-
#include "test_swig_dict.h"
#include "typemgr.h"
#include "xmlstream.h"

const char* persister<AA ,xmlstream>::className_ = "AA";
const char* persister<AA ,xmlstream>::memberNames_[] = {};
const char* persister<AA ,xmlstream>::memberTypes_[] = {};
type_registrar<AA,xmlstream> __0__registration__;

const char* persister<Point ,xmlstream>::className_ = "Point";
const char* persister<Point ,xmlstream>::memberNames_[] = {"x_", "y_"};
const char* persister<Point ,xmlstream>::memberTypes_[] = {"float8", "float8"};
type_registrar<Point,xmlstream> __1__registration__;

const char* persister<Star ,xmlstream>::className_ = "Star";
const char* persister<Star ,xmlstream>::memberNames_[] = {"id_", "flux_"};
const char* persister<Star ,xmlstream>::memberTypes_[] = {"uint4", "float8"};
type_registrar<Star,xmlstream> __2__registration__;

const char* persister<B<int> ,xmlstream>::className_ = "B<int>";
const char* persister<B<int> ,xmlstream>::memberNames_[] = {"t_"};
const char* persister<B<int> ,xmlstream>::memberTypes_[] = {"int"};
type_registrar<B<int>,xmlstream> __3__registration__;

const char* persister<B<double> ,xmlstream>::className_ = "B<double>";
const char* persister<B<double> ,xmlstream>::memberNames_[] = {"t_"};
const char* persister<B<double> ,xmlstream>::memberTypes_[] = {"double"};
type_registrar<B<double>,xmlstream> __4__registration__;

const char* persister<BB<double,double> ,xmlstream>::className_ = "BB<double,double>";
const char* persister<BB<double,double> ,xmlstream>::memberNames_[] = {"lt_", "mtu_", "t_", "u_"};
const char* persister<BB<double,double> ,xmlstream>::memberTypes_[] = {"std::list<double>", "std::map<double,double>", "double", "double"};
type_registrar<BB<double,double>,xmlstream> __5__registration__;

const char* persister<BB<int,double> ,xmlstream>::className_ = "BB<int,double>";
const char* persister<BB<int,double> ,xmlstream>::memberNames_[] = {"lt_", "mtu_", "t_", "u_"};
const char* persister<BB<int,double> ,xmlstream>::memberTypes_[] = {"std::list<int>", "std::map<int,double>", "int", "double"};
type_registrar<BB<int,double>,xmlstream> __6__registration__;

const char* persister<BB<string,string> ,xmlstream>::className_ = "BB<string,string>";
const char* persister<BB<string,string> ,xmlstream>::memberNames_[] = {"lt_", "mtu_", "t_", "u_"};
const char* persister<BB<string,string> ,xmlstream>::memberTypes_[] = {"std::list<string>", "std::map<string,string>", "string", "string"};
type_registrar<BB<string,string>,xmlstream> __7__registration__;

