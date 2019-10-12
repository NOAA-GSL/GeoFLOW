//==================================================================================
// Module       : gtstat.hpp
// Date         : 10/1/19 (DLR)
// Description  : Encapsulates the methods and data associated with
//                GeoFLOW statistics methods.
// Copyright    : Copyright 2019, Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#if !defined(_GTSTAT_DECL_HPP)
#define _GTSTAT_DECL_HPP

#include <cstdlib>
#include <limits>
#include <iostream>
#include <vector>
#include <cassert>
#include "gtypes.h"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gcomm.hpp"


using namespace std;


template <typename T> class GTStat
{
  static_assert(std::is_floating_point<T>::value, "Requires floating point template parameter");

  public:

           GTStat<T>() = delete;
           GTStat<T>(GSIZET nbins, GC_COMM icomm=GC_COMM_WORLD);
           GTStat<T>(const GTStat<T> &obj);
          ~GTStat<T>() = default;
    
    void   set_nbins(GSIZET nbins)
             {nbins_ = nbins; 
              lpdf_.resize(nbins);}        // Set no. bins
    GSIZET get_nbins() const
             {return nbins_; }             // Get no. bins
    void   dopdf1d(GTVector<T> &u, GBOOL ifixdr, T &fmin, T &fmax, GINT iside, GBOOL dolog, GTVector<T> &utmp, GTVector<T> &pdf);
    void   dopdf1d(GTVector<T> &u, GBOOL ifixdr, T &fmin, T &fmax, GINT iside, GBOOL dolog, GTVector<T> &utmp, const GString &fname); 
    

  private:

    GINT             myrank_;   // task's rank
    GC_COMM          comm_;     // GC_COMM handle
    GSIZET           nbins_;    // no. bins
    GSIZET           nkeep_;    // no. stoch var indices kept
    T                gavg_;     // avg of PDF
    T                sig_;      // std deviation of PDF
    GTVector<GSIZET> ikeep_;    // indirection indices to valid values
    GTVector     <T> lpdf_ ;    // local pdf
    GTVector     <T> gpdf_ ;    // globally-reduced pdf

};


#include "gtstat.ipp"

#endif

