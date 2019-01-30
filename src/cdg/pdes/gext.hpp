//==================================================================================
// Module       : gext.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                an extrapolation scheme with variable timestep.
//                PDE. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(G_EXT_HPP)
#define G_EXT_HPP

#include "gtvector.hpp"
#include "multilev_coeffs.hpp"



class G_EXT 
{
public:
                           G_EXT(GSIZET iorder=2);
                          ~G_EXT();
                           G_EXT(const G_EXT &a);

         friend ostream&    operator<<(ostream&, G_EXT &);                                  // Output stream operator



private:
// Private methods:
         void               computeCoeffs();

};
#endif

