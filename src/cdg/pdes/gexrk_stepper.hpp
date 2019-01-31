//==================================================================================
// Module       : gexrk_stepper.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGExRKSTEPPER_HPP)
#define _GGExRKSTEPPER_HPP

#include "gtvector.hpp"
#include "gtmatrix.hpp"


template class<typename T>
class GExRKstepper
{
public:
                           GExRKstepper();
                           GExRKstepper(GSIZET iorder=4);
                          ~GExRKstepper() = default;
                           GExRKstepper(const GExRKstepper &a) = default;
                           GExRKstepper &operator=(const GExRKstepper &bu) = default;
                           setOrder(INT iorder);
                           setRHSfunction(std::function<void(const GFTYPE &t, 
                                          GTVector<GTVector<GFTYPE>> &uin,
                                          GFTYPE &dt, GTVector<GTVector<GFTYPE>> &dudt)> &callback);
                           step(const GFTYPE &t, GTVector<GTVector<GFTYPE>> &uin,
                                GFTYPE &dt, GTVector<GTVector<GFTYPE>> &uout);


private:
// Private methods:

         GButcherRK         butcher_;
         std::function<void(const GFTYPE &t, 
                            GTVector<GTVector<GFTYPE>> &uin,
                            GFTYPE &dt, GTVector<GTVector<GFTYPE>> &dudt)>
                            *rhd_callback_;

};
#endif

