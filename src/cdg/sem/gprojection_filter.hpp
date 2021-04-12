//==================================================================================
// Module       : gprojection_filter.hpp
// Date         : 4/10/21 (DLR)
// Description  : Computes a projectioon filter for stabilization.
//                Taken from Deville, Fischer & Mund "High-Order
//                Methods for Incompressible Flow"
//                 
//                Define interpolation matrices,
//                    I_N^M(i,j) = h_N,j(xi_M,i)
//                where h_N is the Lagrange interpolating polynomial
//                of order N, evaluated at nodes, x_M,i from the Mth 
//                polynomial. Then define
//                    P_N^M = I_M^N I_N^M.
//                The filter is then defined in 1d as
//                    F = alpha Pi_N^M + 1-alpha) I_N^N
//                where
//                    I_N^N 
//                is the Nth-order identify matrix. Filter, F, is then 
//                applied in tensor product form.
// Copyright    : Copyright 2021. Colorado State University. All rights reserved.
// Derived From : FilterBase
//==================================================================================

#if !defined(_GPROJECTIONFILTER_HPP)
#define _GPROJECTIONFILTER_HPP

#include "gtypes.h"
#include "gtvector.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/filter_base.hpp"


template<typename TypePack>
class GProjectionFilter : public FilterBase<TypePack>
{
public:
        using Types      = EquationBase<TypePack>;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using Grid       = typename Types::Grid;
        using Mass       = typename Types::Mass;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;
        using Size       = typename Types::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
               "StateComp is of incorrect type");

        // GProjectionFilter traits:
        struct Traits {
          std::vector<int>    istate;   // state ids to filter
          std::vector<int>    pdelta;   // amount to reduce orig order
          std::vector<double>  alpha;   // filter 'strength'
        };


                          GProjectionFilter() = delete;
                          GProjectionFilter(Traits &traits, Grid &grid);
                          GProjectionFilter(const GProjectionFilter &);
                         ~GProjectionFilter();

protected:

        void              apply_impl(const Time &t, State &u, State  &utmp, 
                                State &uo);
        void              apply_impl(const Time &t, State &u, State  &utmp); 

private:
        void              init();

        GBOOL                        bInit_;   // is filter initialized?
        Traits                       traits_;
        GTVector<Ftype>              tmp_;     // 1-element tmp space
        GTVector<GTMatrix<Ftype>>    F_;       // interp matrices
        GTVector<GTMatrix<Ftype>>    FT_;      // interp matrix transposes
        Grid                        *grid_;    // grid set on construction

};


#include "gprojection_filter.ipp"


#endif
