//==================================================================================
// Module       : gext.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                an extrapolation scheme with variable timestep.
//                PDE. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_G_EXT_HPP)
#define _G_EXT_HPP

#include "gtvector.hpp"
#include "gmultilev_coeffs_base.hpp"


template<typename T>
class G_EXT : GMultilevel_coeffs_base<T>
{
public:
                           G_EXT(GSIZET iorder=2);
                          ~G_EXT();
                           G_EXT(const G_EXT &a);



private:
// Private methods:
         void               computeCoeffs();

};
#endif

