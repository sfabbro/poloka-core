// -*- C++ -*-
// 
// streamer.h
// 
// a streamer makes the interface between old versions
// of an object written on file, and its current implementation.
// Usually, the end-user just wants to use obj_streamer, a generic 
// streamer, able to read the object description from the file
// and interface it to the current object implementation.
// However, if no dictionary is available, and if the object 
// was written in a pretty obfuscated way, one can 
// write its own streamer.
// 
#ifndef STREAMER_H
#define STREAMER_H



template<class T>
class generic_streamer {
public:
  generic_streamer();
  ~generic_streamer();
  
  virtual unsigned int version() const { version_; }
  
  template<class IOS>
  void                 read(obj_input<IOS> const oi&, persister<T>& p) {
    unsigned int i, sz = obj_map_.size();
    setter_ s;
    for(i=0;i<sz_;i++) {
      s = persister<T>::set_[obj_map_[i]];
      p.(*s)(oi);
    }
  }
  
  void                 copy(generic_streamer<T> const& gs) {
    clear_();
    std::copy(gs.obj_map_.begin(), gs.obj_map_.end(), 
	      std::back_inserter(obj_map_));
    std::copy(gs.names_.begin(), gs.names_.end(), 
	      std::back_inserter(names_))
  }
  
private:
  typedef persister<T>::setter_ setter_;
  
  void  clear_() {
    obj_map_.clear();
    names_.clear();
  }

  unsigned int version_;
  std::vector<unsigned int> obj_map_;
  std::vector<std::string>  names_;
};


#endif

