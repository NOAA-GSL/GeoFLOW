//==================================================================================
// Module       : gab.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                and Adams-Bashforth scheme with variable timestep.
//                PDE. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(G_AB_HPP)
#define G_AB_HPP

#include "gtvector.hpp"
#include "multilev_coeffs.hpp"


class G_AB: public G_Multilevel_coeffs_base 
{
public:
                           G_AB(GSIZET iorder=2);
                          ~G_AB();
                           G_AB(const G_AB &a);

private:
// Private methods:
         void               computeCoeffs();

};
#endif

