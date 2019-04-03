//==================================================================================
// Module       : gposixio_observer.hpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based 
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GGPOSIXIO_OBSERVER_HPP)
#define _GPOSIXIO_OBSERVER_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;

typedef GTVector<GTVector<GFTYPE>*> State;
typedef GTVector<GFTYPE> StateElem;


extern void gio_write(const GGrid &grid, const State &u, const GTVector<GINT> &nu, 
                      const GSIZET tindex, const GFTYPE time, const GTVector<GString> &svars, 
                      GC_COMM comm, GBOOL &bprgrid);


template<typename EquationType>
class GPosixIOObserver : public ObserverBase<EquationType>
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

//      using ObserverBase<EquationType>::ObsType;
//      using OBS_CYCLE = typename ObserverBase<EquationType>::ObsType::OBS_CYCLE;
//      using OBS_TIME  = typename ObserverBase<EquationType>::OBS_TIME;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

                           GPosixIOObserver() = delete;
                           GPosixIOObserver(typename ObserverBase<EquationType>::Traits &traits, Grid &grid);
                          ~GPosixIOObserver() = default;
                           GPosixIOObserver(const GPosixIOObserver &a) = default;
                           GPosixIOObserver &operator=(const GPosixIOObserver &bu) = default;

        void               observe_impl(const Time &t, const State &u, const State &uf);

private:
// Private methods:
        void               init(const Time t, const State &u);
// Private data:
        GBOOL              bprgrid_;    // print grid flag
        GSIZET             cycle_last_; // most recent output cycle
        GSIZET             cycle_;      // continuously-running cycle
        GSIZET             ocycle_;     // output cycle number
        GFTYPE             time_last_;  // most recent output time
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of states to print
        GString            sdir_;      ;// directory in which to write
    

};

#include "gposixio_observer.ipp"

#endif

