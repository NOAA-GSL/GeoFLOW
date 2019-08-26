//==================================================================================
// Module       : gposixio_observer.hpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based 
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GPOSIXIO_OBSERVER_HPP)
#define _GPOSIXIO_OBSERVER_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "gio.h"

using namespace geoflow::pdeint;
using namespace std;


#if 0
extern GIOTraits;
void gio_write_state(GIOTraits &, GGrid &grid, 
                     const GTVector<GTVector<GFTYPE>*> &u,
                     const GTVector<GINT> &iu, 
                     const GTVector<GString> &svars,
                     GC_COMM comm);

void gio_write_grid (GIOTraits &, GGrid &grid,
                     const GTVector<GString> &svars,
                     GC_COMM comm);
#endif


template<typename EquationType>
class GPosixIOObserver : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename Equation::State;
        using Grid        = typename Equation::Grid;
        using Value       = typename Equation::Value;
        using Derivative  = typename Equation::Derivative;
        using Time        = typename Equation::Time;
        using Jacobian    = typename Equation::Jacobian;
        using Size        = typename Equation::Size;

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
                           GPosixIOObserver(const EqnBasePtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits);
                          ~GPosixIOObserver() = default;
                           GPosixIOObserver(const GPosixIOObserver &a) = default;
                           GPosixIOObserver &operator=(const GPosixIOObserver &bu) = default;

        void               observe_impl(const Time &t, const State &u, const State &uf);

private:
// Private methods:
        void               init(const Time t, const State &u);
// Private data:
        GBOOL              bprgrid_;    // print grid flag
        GBOOL              bInit_;      // is initialized?
        GINT               ivers_;      // output version number
        GINT               wtime_;      // width of time field
        GINT               wtask_;      // width of task field
        GINT               wfile_;      // filename max
        GSIZET             cycle_last_; // most recent output cycle
        GSIZET             cycle_;      // continuously-running cycle
        GSIZET             ocycle_;     // output cycle number
        GFTYPE             time_last_;  // most recent output time
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of state files
        GTVector<GString>  grid_names_ ;// list of names of grid comp files
        GString            sidir_;     ;// directory from which to read (e.g., for restart)
        GString            sodir_;     ;// directory in which to write
    

};

#include "gposixio_observer.ipp"

#endif

