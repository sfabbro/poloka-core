#ifndef SPARSEVECT__H
#define SPARSEVECT__H



struct SparseVectElement
{
  public:
  int index;
  double value;

  SparseVectElement(const int &I, const double &V) : index(I), value(V) {};
  
  bool operator < (const SparseVectElement &R) const
  { return (index < R.index);}
};



#include <vector>

class SparseVect : public std::vector<SparseVectElement>
{
 private:

#ifdef STORAGE
  struct ElementRef{
    SparseVect &v;
    int i;
    ElementRef(SparseVect &V, int I) : v(V), i(I) {};
    void operator =(const double &Val)
    { v.set(i,Val);}
  };
#endif

 public:

  SparseVect() {sorted = true;} // not really necessary.
  void zero() {clear();}

#ifdef STORAGE  
  double operator() (const int i) const;
  ElementRef operator()(const int i)
    { return ElementRef(*this,i);}
#endif

  double set(const int I, const double Val)
    { push_back(SparseVectElement(I,Val)); sorted = false; return Val;}

  void sort()
    {
    if (!sorted) std::sort(begin(), end());
    sorted = true;
    }
  
  iterator get(const int index)
  {
    iterator I;
    for(I=begin();I!=end();I++)
      if( I->index == index)
	return I;
    return end();
  }

  const_iterator get(const int index) const
  {
    const_iterator I;
    for(I=begin();I!=end();I++)
      if( I->index == index)
	return I;
    return end();
  }


  double operator()(const int ii) const
  {
    for (const_iterator i = begin(); i != end(); ++i)
      if (i->index == ii) return i->value;
    return 0;
  }
  
 protected:
  bool sorted;


};

typedef SparseVect::const_iterator SVCIterator;
typedef SparseVect::iterator SVIterator;
    

#endif /* SPARSEVECT__H */
