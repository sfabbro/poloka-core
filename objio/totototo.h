// -*- C++ -*-
#ifndef __totototo__
#define __totototo__

#include <string>
#include "persister.h"
#include "objio.h"
#include "test_swig.h"

template<>
class persister<A> : public handle<A> {
public:
persister<A>() : handle<A>() {}
persister<A>(A& obj) : handle<A>(&obj) {}
persister<A>(A const & obj) : handle<A>(const_cast<A*>(&obj)) {}// I know, that's UGLY.
~persister<A>() {}

unsigned int version() const { return 0; }
std::string  name()    const { return (std::string)"A"; }
unsigned int size()    const { return 0;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, A const& p)
{
persister<A> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, A& p)
{
persister<A> pp(p);
oi.read(pp);
return oi;
};
template<>
class persister<C> : public handle<C> {
public:
persister<C>() : handle<C>() {}
persister<C>(C& obj) : handle<C>(&obj) {}
persister<C>(C const & obj) : handle<C>(const_cast<C*>(&obj)) {}// I know, that's UGLY.
~persister<C>() {}

unsigned int version() const { return 0; }
std::string  name()    const { return (std::string)"C"; }
unsigned int size()    const { return 0;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, C const& p)
{
persister<C> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, C& p)
{
persister<C> pp(p);
oi.read(pp);
return oi;
};
template<>
class persister<AA> : public handle<AA> {
public:
persister<AA>() : handle<AA>() {}
persister<AA>(AA& obj) : handle<AA>(&obj) {}
persister<AA>(AA const & obj) : handle<AA>(const_cast<AA*>(&obj)) {}// I know, that's UGLY.
~persister<AA>() {}

unsigned int version() const { return 1; }
std::string  name()    const { return (std::string)"AA"; }
unsigned int size()    const { return 0;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, AA const& p)
{
persister<AA> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, AA& p)
{
persister<AA> pp(p);
oi.read(pp);
return oi;
};
template<>
class persister<Point> : public handle<Point> {
public:
persister<Point>() : handle<Point>() {}
persister<Point>(Point& obj) : handle<Point>(&obj) {}
persister<Point>(Point const & obj) : handle<Point>(const_cast<Point*>(&obj)) {}// I know, that's UGLY.
~persister<Point>() {}

unsigned int version() const { return 1; }
std::string  name()    const { return (std::string)"Point"; }
unsigned int size()    const { return 2;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
oo.write(obj_->x_, "x_");
oo.write(obj_->y_, "y_");
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
oi.read(obj_->x_);
oi.read(obj_->y_);
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, Point const& p)
{
persister<Point> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, Point& p)
{
persister<Point> pp(p);
oi.read(pp);
return oi;
};
template<>
class persister<Star> : public handle<Star> {
public:
persister<Star>() : handle<Star>() {}
persister<Star>(Star& obj) : handle<Star>(&obj) {}
persister<Star>(Star const & obj) : handle<Star>(const_cast<Star*>(&obj)) {}// I know, that's UGLY.
~persister<Star>() {}

unsigned int version() const { return 1; }
std::string  name()    const { return (std::string)"Star"; }
unsigned int size()    const { return 2;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
oo.write(obj_->id_, "id_");
oo.write(obj_->flux_, "flux_");
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
oi.read(obj_->id_);
oi.read(obj_->flux_);
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, Star const& p)
{
persister<Star> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, Star& p)
{
persister<Star> pp(p);
oi.read(pp);
return oi;
};
template<>
class persister<B> : public handle<B> {
public:
persister<B>() : handle<B>() {}
persister<B>(B& obj) : handle<B>(&obj) {}
persister<B>(B const & obj) : handle<B>(const_cast<B*>(&obj)) {}// I know, that's UGLY.
~persister<B>() {}

unsigned int version() const { return 0; }
std::string  name()    const { return (std::string)"B"; }
unsigned int size()    const { return 1;}
std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }
std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }

private:
template<class IOS>
void write_members(obj_output<IOS>& oo) const {
oo.write(obj_->t_, "t_");
}

template<class IOS>void read_members(obj_input<IOS> const& oi) {
oi.read(obj_->t_);
}

template<class T> friend class obj_input;
template<class T> friend class obj_output;

static const char* memberNames_[];
static const char* memberTypes_[];
};

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, B const& p)
{
persister<B> pp(p);
oo.write(pp);
return oo;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, B& p)
{
persister<B> pp(p);
oi.read(pp);
return oi;
};
