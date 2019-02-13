//==================================================================================
// Module       : gbdf.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                a Backwards Differencing Formula (BDF) scheme with variable 
//                timestep.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_G_BDF_HPP)
#define _G_BDF_HPP

#include "gtmatrix.hpp"
#include "gmultilev_coeffs_base.hpp"


template<typename T>
class G_BDF: public GMultilevel_coeffs_base<T>
{
public:
                           G_BDF(GSIZET iorder=2);
                          ~G_BDF();
                           G_BDF(const G_BDF &a);

private:
// Private methods:
         void               computeCoeffs();
         GTMatrix<T>        c_bdf_;              // database of BDF coeffs

};
#endif

