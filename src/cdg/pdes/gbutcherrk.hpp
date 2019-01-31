//==================================================================================
// Module       : gbutcherrk.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multistage Butcher tableau for
//                explicit RK methods. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(G_BUTCHERRK_HPP)
#define G_BUTCHERRK_HPP

#include "gtvector.hpp"
#include "gtmatrix.hpp"


template class<typename T>
class GButcherRK
{
public:
                           GButcherRK();
                           GButcherRK(GSIZET iorder=4);
                          ~GButcherRK() = default;
                           GButcherRK(const GButcherRK &a) = default;
                           GButcherRK &operator=(const GButcherRK &bu) = default;
                           setOrder(GINT iorder);

         GTVector<T>      &alpha() { return alpha_ };
         GTVector<T>      &beta () { return beta_ };
         GTVector<T>      &c()     { return c_ };

private:
// Private methods:
         void               computeCoeffs();

         GTVector<T> alpha_; // time coeffs
         GTMatrix<T> beta_;  // stage coeffs
         GTVector<T> c_;     // stage weights
};
#endif

