//==================================================================================
// Module       : gout_simple_observer.hpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based 
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GOUT_SIMPLE_OBS_HPP)
#define _GOUT_SIMPLE_OBS_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gio.h"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename EquationType>
class GOutSimpleObserver : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using State       = typename Equation::State;
        using Grid        = typename Equation::Grid;
        using Value       = typename Equation::Value;
        using Derivative  = typename Equation::Derivative;
        using Time        = typename Equation::Time;
        using Jacobian    = typename Equation::Jacobian;
        using Size        = typename Equation::Size;
        using EquationPtr = std::shared_ptr<Equation>;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

                           GOutSimpleObserver(Traits &traits, Grid &grid);
                          ~GOutSimpleObserver();
                           GOutSimpleObserver(const GOutSimpleObserver &a) = default;
                           GOutSimpleObserver &operator=(const GOutSimpleObserver &bu) = default;

        void               observe_impl(const Time t, const State &u);

private:
// Private methods:
        void               init(State &u);
// Private data:
        GBOOL              bgrid_printed_;
        GSIZET             cycle_last_;
        GSIZET             cycle_;
        Time               time_last_;
        GTVector<GINT>     istate_;
        GTVector<GString>  state_names_;
        GGrid             *grid_;
    

};

#include "gout_simple_observer.ipp"

#endif

