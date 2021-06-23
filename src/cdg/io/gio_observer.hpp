//==================================================================================
// Module       : gio_observer.ipp
// Date         : 1/28/20 (DLR)
// Description  : Observer object for carrying out binary output
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GIO_OBSERVER_HPP)
#define _GIO_OBSERVER_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/io_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;


template<typename EquationType>
class GIOObserver : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using IOBaseType  = IOBase<EquationType>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using State       = typename Equation::State;
        using StateInfo   = typename Equation::StateInfo;
        using Grid        = typename Equation::Grid;
        using Ftype       = typename Equation::Ftype;
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

                           GIOObserver() = delete;
                           GIOObserver(EqnBasePtr &equation, Grid &grid, const IOBasePtr &io_ptr,
                                       typename ObserverBase<EquationType>::Traits &traits);

                          ~GIOObserver() = default;
                           GIOObserver(const GIOObserver &a) = default;
                           GIOObserver &operator=(const GIOObserver &bu) = default;

        void               observe_impl(const Time &t, const Time &dt, const State &u, const State &uf);
        void               init_impl(StateInfo &);
        void               setIO(IOBasePtr ioobj) { pIO_ = ioobj; }

private:
// Private methods:
        void               print_derived(const Time &t, const State &u);
        void               print_griddat(StateInfo &);
// Private data:
        GBOOL              bprgrid_    ;// print grid flag
        GBOOL              bInit_      ;// is initialized?
        GSIZET             cycle_last_ ;// most recent output cycle
        GSIZET             cycle_      ;// continuously-running cycle
        GSIZET             ocycle_     ;// output cycle number
        Ftype              time_last_  ;// most recent output time
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of state files
        GTVector<GString>  grid_names_ ;// list of names of grid comp files
        IOBasePtr          pIO_        ;// ptr to IO object
        State              up_         ;// current state array
        State              gp_         ;// current grid 'state' array
        StateInfo          stateinfo_  ;// info struct for state
        StateInfo          gridinfo_   ;// info struct for grid
        EqnBasePtr         pEqn_       ;// equation object pointer

};

#include "gio_observer.ipp"

#endif

