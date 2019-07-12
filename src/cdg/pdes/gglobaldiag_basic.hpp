//==================================================================================
// Module       : gglobaldiag_basic.hpp
// Date         : 3/28/19 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                Burgers equation.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GGLOBALDIAG_BURGERS_OBS_HPP)
#define _GGLOBALDIAG_BURGERS_OBS_HPP

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


template<typename EquationType>
class GGlobalDiag_basic : public ObserverBase<EquationType>
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
        using ObserverBase<EquationType>::utmp_;
        using ObserverBase<EquationType>::grid_;
        using ObserverBase<EquationType>::traits_;

//      using ObserverBase<EquationType>::ObsType;
//      using OBS_CYCLE = typename ObserverBase<EquationType>::ObsType::OBS_CYCLE;
//      using OBS_TIME  = typename ObserverBase<EquationType>::OBS_TIME;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

                           GGlobalDiag_basic() = delete;
                           GGlobalDiag_basic(EquationPtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits);
                          ~GGlobalDiag_basic() = default;
                           GGlobalDiag_basic(const GGlobalDiag_basic &a) = default;
                           GGlobalDiag_basic &operator=(const GGlobalDiag_basic &bu) = default;

        void               observe_impl(const Time &t, const State &u, const State &uf);

private:
// Private methods:
        void               init(const Time t, const State &u);
        void               do_kinetic_L2 (const Time t, const State &u, const State &uf, const GString file);
        void               do_kinetic_max(const Time t, const State &u, const State &uf, const GString file);
// Private data:
        GBOOL              bInit_;
        GINT               myrank_;     // MPI rank
        GSIZET             cycle_last_; // most recent output cycle
        GSIZET             cycle_;      // continuously-running cycle
        GSIZET             ocycle_;     // output cycle number
        GSIZET             ikinetic_;   // stores GSC_KINETIC component types
        GFTYPE             time_last_;  // most recent output time
        GFTYPE             ivol_;       // inverse of grid volume
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of states to print
        GString            sidir_;      // directory from which to read
        GString            sodir_;      // directory in which to write
        State              ku_;         // vector of pointers to kinetic components

};

#include "gglobaldiag_basic.ipp"

#endif

