//==================================================================================
// Module       : multilev_coeffs_base.hpp
// Date         : 1/28/19 (DLR)
// Description  : Base class for multilevel time stepping coefficients.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(G_MULTLEVEL_COEFFS_HPP)
#define G_MULTLEVEL_COEFFS_HPP

#include "gtvector.hpp"


template<typename T>
class GMultilevel_coeffs_base 
{
public:
                           GMultilevel_coeffs_base() = delete;
                           GMultilevel_coeffs_base(GINT iorder, GTVector<T> &dthist)
                           { iorder_ = iorder; dthist_ = &dthist; }
                           GMultilevel_coeffs_base(const GMultilevel_coeffs_base &a) = default;
                           virtual ~GMultilevel_coeffs_base() = default;

         GTVector<T>       &getCoeffs() { return coeffs_ ; }
         GTVector<T>       &getTimestepHistory() { return *dthist_; }

         T                 &operator()(const GSIZET  i) { return coeffs_[i]; }
         T                 &operator[](const GSIZET  i) { return coeffs_[i]; }


protected:

// Private data:
         GINT               iorder_;              // specified GMultilevel_coeffs_base order
         GINT               maxorder_;            // max GMultilevel_coeffs_base order allowed by class
         GTVector<T>        coeffs_;              // buffer of GMultilevel_coeffs_base coeffs
         GTVector<T>       *dthist_;              // timestep history pointer, managed outside of class
};
#endif

