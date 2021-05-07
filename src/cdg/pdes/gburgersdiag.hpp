//==================================================================================
// Module       : gburgersdiag.hpp
// Date         : 3/28/19 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                Burgers equation.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GBURGERSDIAG_OBS_HPP)
#define _GBURGERSDIAG_OBS_HPP

#include "gtvector.hpp"
#include "gutils.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;



template<typename EquationType>
class GBurgersDiag : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename Equation::State;
        using StateInfo   = typename Equation::StateInfo;
        using Grid        = typename Equation::Grid;
        using Ftype       = typename Equation::Ftype;
        using Time        = typename Equation::Time;
        using CompDesc    = typename Equation::CompDesc;
        using Size        = typename Equation::Size;
        using ObserverBase<EquationType>::utmp_;
        using ObserverBase<EquationType>::traits_;


//      using ObserverBase<EquationType>::ObsType;
//      using OBS_CYCLE = typename ObserverBase<EquationType>::ObsType::OBS_CYCLE;
//      using OBS_TIME  = typename ObserverBase<EquationType>::OBS_TIME;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

                           GBurgersDiag() = delete;
                           GBurgersDiag(EqnBasePtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits);
                          ~GBurgersDiag() = default;
                           GBurgersDiag(const GBurgersDiag &a) = default;
                           GBurgersDiag &operator=(const GBurgersDiag &bu) = default;

        void               observe_impl(const Time &t, const Time &dt, const State &u, const State &uf);

        void               init_impl(StateInfo &);
private:
// Private methods:
        void               do_kinetic_L2 (const Time &t, const Time &dt, const State &u, const State &uf, const GString file);
        void               do_kinetic_max(const Time &t, const Time &dt, const State &u, const State &uf, const GString file);
// Private data:
        GBOOL              bInit_;
        GINT               myrank_;     // MPI rank
        GSIZET             cycle_last_; // most recent output cycle
        GSIZET             cycle_;      // continuously-running cycle
        GSIZET             ocycle_;     // output cycle number
        GTVector<GINT>     ikinetic_;   // stores GSC_KINETIC component types
        Ftype              time_last_;  // most recent output time
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of states to print
        GString            sidir_;      // directory from which to read
        GString            sodir_;      // directory in which to write
        State              ku_;         // vector of pointers to kinetic components
        Grid              *grid_;       // grid object

};

#include "gburgersdiag.ipp"

#endif

