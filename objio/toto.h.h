// -*- C++ -*-
#ifndef __toto__
#define __toto__

#include <string>
#include "persister.h"
#include "objio.h"
#include "test_swig2.h"

template<>
class persister<test_1> : public handle<test_1> {
public:
persister<test_1>() : handle<test_1>() {}
persister<test_1>(test_1& obj) : handle<test_1>(&obj) {}
persister<test_1>(test_1 const & obj) : handle<test_1>(const_cast<test_1*>(&obj)) {}// I know, that's UGLY.
~persister<test_1>() {}

unsigned int version() const { return 0; }
std::string  name()    const { return (std::string)"test_1"; }
unsigned int size()    const { return 9;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
oo.write(obj_->a1, "a1");
oo.write(obj_->toto, "toto");
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
oi.read(obj_->a1);
oi.read(obj_->toto);
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, test_1 const& p)
{
persister<test_1> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, test_1& p)
{
persister<test_1> pp(p);
oi.read(pp);
return oi;
};
