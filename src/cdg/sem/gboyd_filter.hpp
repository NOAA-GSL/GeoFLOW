//==================================================================================
// Module       : gboyd_filter.hpp
// Date         : 9/14/20 (DLR)
// Description  : Computes the Boyd filter to diminish aliasing errors.
//                Taken from Giraldo & Rosemont 2004, MWR:132 133:
//                    u <-- F u
//                where
//                    F = L Lambda L^-1; s.t.
//                and 
//                    Lambda = 1 if i< ifilter
//                             mu [(i-ifilter)/(N - ifilter)]^2, i>= ifilter.
//                L is the Legendre transform matrix:
//                    L = | P_0(xi0), P_1(xi0) ... P_i(xi0)-P_{i-2)(xi0) ... |
//                        | P_0(xi1), P_1(xi1) ... P_i(xi1)-P_{i-2)(xi1) ... |
//                        |  ...                                             |.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : FilterBase
//==================================================================================

#if !defined(_GBOYDFILTER_HPP)
#define _GBOYDFILTER_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/filter_base.hpp"


template<typename TypePack>
class GBoydFilter : public FilterBase<TypePack>
{
public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Mass       = typename Interface::Mass;
        using Ftype      = typename Interface::Ftype;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
               "StateComp is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Ftype>*>>::value,
               "Derivative is of incorrect type");

        // GBoydFilter traits:
        struct Traits {
          int    ifilter;    // filter starting mode
          double mufilter;   // filter trunc amount
        };


                          GBoydFilter() = delete;
                          GBoydFilter(Traits &traits, Grid &grid);
                          GBoydFilter(const GBoydFilter &);
                         ~GBoydFilter();

protected:

        void              apply_impl(const Time &t, StateComp &u, State  &utmp, 
                                StateComp &po);
        void              apply_impl(const Time &t, StateComp &u, State  &utmp); 

private:
        void              init();

        GBOOL                         bInit_;    // is filter initialized?
        GINT                          ifilter_;  // filter mode
        Ftype                         mufilter_; // truncation amount 
        GTMatrix<Ftype>               Lambda_;   // mode-weighting matrix
        Grid                         *grid_;     // grid set on construction


};


#include "gboyd_filter.ipp"


#endif
