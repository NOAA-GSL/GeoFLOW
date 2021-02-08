//==================================================================================
// Module       : gdd_base.hpp
// Date         : 8/28/18 (DLR)
// Description  : Forms virtual abstract base class for domain 
//                decomposition objects. Default behaviour is to
//                do a simple partition that divides the number of 
//                elememt representations by the number of MPI tasks.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(_GDD_BASE_HPP)
#define _GDD_BASE_HPP
#include "gtypes.h"
#include "gtpoint.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"

#include "tbox/tracer.hpp"

template<typename T>
class GDD_base
{

public:

                          GDD_base(GINT nprocs);
                          GDD_base(const GDD_base &);
                         ~GDD_base();

virtual  GDD_base        &operator=(const GDD_base &g) { nprocs_ = g.nprocs_; return *this; }

                   
virtual  GSIZET           doDD(const GTVector<GTVector<T>> &x, GINT irank, GTVector<GINT> &iret );
virtual  GSIZET           doDD(const GTVector<GTPoint<T>>  &x, GINT irank, GTVector<GINT> &iret);

private:

GINT               nprocs_ ;  // number of tasks/ranks to partition into

};

#include "gdd_base.ipp"

#endif
