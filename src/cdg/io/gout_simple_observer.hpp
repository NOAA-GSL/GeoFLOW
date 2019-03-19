//==================================================================================
// Module       : ggio_simple.hpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object fr carrying out simple non-MPI-based 
//                I/O
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGIO_SIMPLE_HPP)
#define _GGIO_SIMPLE_HPP

#include "gtvector.hpp"
#include "observer_base.hpp"



template <typename EquationType>
class GGIOSimple : public ObserverBase<EquationType>
{

public:
                           GGIOSimple() = default;
                           GGIOSimple(GSIZET nstage);
                          ~GGIOSimple();
                           GGIOSimple(const GGIOSimple &a) = default;
                           GGIOSimple &operator=(const GGIOSimple &bu) = default;

        void               observe_impl(const Traits &traits, const Time &t, const State &u);

private:
// Private methods:

// Private data:
        GBOOL              bgrid_printed_;

};

#include "ggio_simple.ipp"

#endif

