//==================================================================================
// Module       : gkeygen.hpp
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining a key-generator, as used in GeoFLOW . This class is
//                intended to be a abstract base class for defined GKeyGen
//                objects
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(GKEYGEN_HPP)
#define GKEYGEN_HPP

#include <stdlib.h>
#include "gtypes.h"
#include "gtpoint.hpp"
#include "gtvector.hpp"

template<typename KT, typename FT> class GKeyGen
{
  static_assert( std::is_integral<KT>::value && std::is_floating_point<FT>::value, "Invalid template type" );

public:

                           GKeyGen() {};
                          ~GKeyGen() {};

virtual  void              key(KT id[], GTPoint<FT> point[], GINT  n)=0;             // Use GTPoint
virtual  void              key(GTVector<KT> &id, GTVector<GTPoint<FT>> &point)=0;    // Use GTPoint
virtual void               key(GTVector<KT> &id, GTVector<GTVector<FT>> &x)=0;        // Use GTVector for points
virtual void               key(GTVector<KT> &id, GTVector<GTVector<FT>> &x,
                               GTVector<GSIZET> &ix)=0;                               // Use ix members of GTVector x for points


private:

         // Member data:
         //   ...
};

#endif
