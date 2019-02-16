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
class G_EXT : public GMultilevel_coeffs_base<T>
{
         using GMultilevel_coeffs_base<T>::iorder_;
         using GMultilevel_coeffs_base<T>::maxorder_;
         using GMultilevel_coeffs_base<T>::coeffs_;
         using GMultilevel_coeffs_base<T>::dthist_;

public:
                           G_EXT(GINT iorder, GTVector<T> &dthist);
                          ~G_EXT();
                           G_EXT(const G_EXT &a);



private:
// Private methods:
         void               computeCoeffs();

};

#include "gext.ipp"

#endif

