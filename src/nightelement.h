// This may look like C code, but it is really -*- C++ -*-
#ifndef NIGHTELEMENT__H
#define NIGHTELEMENT__H

class Night;

//! a template to use when an element belongs to a night
template <class Element> class NightElement : public Element {
public:
  //! empty constructor
  NightElement(): night(NULL) {};
  NightElement(const Element &AnElement): Element(AnElement), night(NULL) {};
  NightElement(const Night *ANight): Element(), night(ANight) {};
  NightElement(const Element &AnElement, const Night *ANight) : Element(AnElement), night(ANight) {};
  ~NightElement(){};
  const Night* night;
  NightElement* Clone() const {return new NightElement(*this);}
  friend bool ByIncreasingSeeing(const NightElement<Element> *one, const NightElement<Element> *two)
  {
    return (one->night->seeing < two->night->seeing);
  }
  
  friend bool ByIncreasingTime(const NightElement<Element> *one, const NightElement<Element> *two)
  {
    return (one->night->julianDate < two->night->julianDate);
  }
};



#endif // NIGHTELEMENT__H
