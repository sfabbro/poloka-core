// -*- C++ -*-
// 
// \file persistence.h 
// 
// 
#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "toadexceptions.h"

#include "objio_defs.h"
#include "toadtypes.h"
#include "dict.h"
#include "persister.h"
#include "objio.h"

#include "xmlexceptions.h"
#include "xmlstream.h"

#define CLASS_VERSION(className,id) \
static const unsigned short __version__=id;\
template<class U,class V>friend class persister;\


#endif

