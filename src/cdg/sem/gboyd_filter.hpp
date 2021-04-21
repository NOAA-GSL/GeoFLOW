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
        using Ftype      = typename Interface::Ftype;
        using Time       = typename Interface::Time;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
               "StateComp is of incorrect type");

        // GBoydFilter traits:
        struct Traits {
          std::vector<GINT>             istate ;    // state ids to filter
          std::vector<GINT>             pdelta;     // starting mode to filter
          std::vector<double>           strength;   // filter strength
        };


                          GBoydFilter() = delete;
                          GBoydFilter(Traits &traits, Grid &grid);
                          GBoydFilter(const GBoydFilter &);
                         ~GBoydFilter();

protected:

        void              apply_impl(const Time &t, State &u, State  &utmp, 
                                State &po);
        void              apply_impl(const Time &t, State &u, State  &utmp); 

private:
        void              init();

        GBOOL                         bInit_;    // is filter initialized?
        GTMatrix<Ftype>               Lambda_;   // mode-weighting matrix
        Traits                        traits_;
        Grid                         *grid_;     // grid set on construction


};


#include "gboyd_filter.ipp"


#endif
