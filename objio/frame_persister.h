// -*- C++ -*-
#ifndef __frame_persister__
#define __frame_persister__

#include <string>
#include "persister.h"
#include "objio.h"
#include "../src/frame.h"

template<>
class persister<Frame> : public handle<Frame> {
public:
persister<Frame>() : handle<Frame>() {}
persister<Frame>(Frame& obj) : handle<Frame>(&obj) {}
persister<Frame>(Frame const & obj) : handle<Frame>(const_cast<Frame*>(&obj)) {}// I know, that's UGLY.
~persister<Frame>() {}

unsigned int version() const { return 2; }
std::string  name()    const { return (std::string)"Frame"; }
unsigned int size()    const { return 4;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
oo.write(obj_->xMin, "xMin");
oo.write(obj_->xMax, "xMax");
oo.write(obj_->yMin, "yMin");
oo.write(obj_->yMax, "yMax");
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
oi.read(obj_->xMin);
oi.read(obj_->xMax);
oi.read(obj_->yMin);
oi.read(obj_->yMax);
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, Frame const& p)
{
persister<Frame> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, Frame& p)
{
persister<Frame> pp(p);
oi.read(pp);
return oi;
};
