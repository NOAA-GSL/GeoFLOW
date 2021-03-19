//==================================================================================
// Module       : gmconvdiag.hpp
// Date         : 3/18/21 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                GMConv solver.
// Copyright    : Copyright 2021. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GMCONVDIAG_OBS_HPP)
#define _GMCONVDIAG_OBS_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gutils.hpp"
#include "gmconv.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;

typedef GTVector<GTVector<GFTYPE>*> State;
typedef GTVector<GFTYPE> StateElem;


template<typename EquationType>
class GMConvDiag : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename Equation::State;
        using StateInfo   = typename Equation::StateInfo;
        using Grid        = typename Equation::Grid;
        using Value       = typename Equation::Value;
        using Derivative  = typename Equation::Derivative;
        using Time        = typename Equation::Time;
        using CompDesc    = typename Equation::CompDesc;
        using Jacobian    = typename Equation::Jacobian;
        using Size        = typename Equation::Size;
        using ObserverBase<EquationType>::utmp_;
        using ObserverBase<EquationType>::traits_;


        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

                           GMConvDiag() = delete;
                           GMConvDiag(EqnBasePtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits);
                          ~GMConvDiag() = default ;
                           GMConvDiag(const GMConvDiag &a) = default;
                           GMConvDiag &operator=(const GMConvDiag &bu) = default;

        void               observe_impl(const Time &t, const State &u, const State &uf);

        void               init_impl(StateInfo &);
private:
// Private methods:
        void               do_L2 (const Time t, const State &u, const State &uf, const GString file);
        void               do_max(const Time t, const State &u, const State &uf, const GString file);
// Private data:
        GBOOL              bInit_;
        GINT               myrank_;     // MPI rank
        GSIZET             cycle_last_; // most recent output cycle
        GSIZET             cycle_;      // continuously-running cycle
        GSIZET             ocycle_;     // output cycle number
        GFTYPE             time_last_;  // most recent output time
        GString            sidir_;      // directory from which to read
        GString            sodir_;      // directory in which to write
        GGrid             *grid_;       // grid object
        GMConv<EquationType> 
                          *solver_;     // specific equation pointer

};

#include "gmconvdiag.ipp"

#endif

