//==================================================================================
// Module       : gbdf.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                a Backwards Differencing Formula (BDF) scheme with variable 
//                timestep.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(G_BDF_HPP)
#define G_BDF_HPP

#include "gtvector.hpp"
#include "multilev_coeffs.hpp"


class G_BDF: public G_Multilevel_coeffs_base 
{
public:
                           G_BDF(GSIZET iorder=2);
                          ~G_BDF();
                           G_BDF(const G_BDF &a);

private:
// Private methods:
         void               computeCoeffs();
         GTMatrix<GFTYPE>   c_bdf_;              // database of BDF coeffs

};
#endif

