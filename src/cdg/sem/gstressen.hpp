//==================================================================================
// Module       : gstress.hpp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the full viscous
//                stress-energy operator. The effect of the viscous stress in the 
//                momentum eqution is
//                    F_i = [2  mu s_{ij}],j + (zeta Div u delta_{ij}),j,
//                where
//                    s_{ij} = (u_j,i + u_i,j)/2 - 1/2d Div u delta_{ij}, and
//                d is the problem dimension. The viscous stress-energy for the 
//                energy equation is
//                    [2 kappa u_i s_{ij}],j - [lambda u_i Div u delta_{ij}],j
//                where u_i is the velocity, and mu, the (shear) viscosity, zeta is
//                the 'bulk' viscosity. Strictly speaking, kappa=mu, and lambda=zeta,
//                but we allow these to be set independently for now. Repeated
//                indices are summed here.  mu, zeta, kappa, lambda, may vary
//                in space or be constant. Currently, the so-called Stokes 
//                approximation is used by default s.t.
//                      (zeta - mu/d) = -2/3 mu, and
//                      (lambda - kappa/d ) = -2/3 kappa.
//               
//                For the energy, this operator is nonlinear, 
//                so it should not derive from GLinOp. 
//
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

#if !defined(_GSTRESSENERGYOP_HPP)
#define _GSTRESSENERGYOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"

#undef  DO_BDY
#undef  DO_COMPRESS_MODES_ONLY

template<typename TypePack>
class GStressEnOp
{
public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Mass       = typename Interface::Mass;
        using Ftype      = typename Interface::Value;
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
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // GStressEn solver traits:
        struct Traits {
          GBOOL       Stokes_hyp = TRUE; // use Stokes hypothesis?
          GBOOL       indep_diss = TRUE; // use indep. diss'n for mom & energy?
          StateComp   mu;                // dynamic/shear viscosity
          StateComp   zeta;              // bulk (dilitation) viscosity
          StateComp   kappa;             // shear visc for energy
          StateComp   lambda;            // bulk visc. for energy
        };

                          GStressEnOp() = delete;
                          GStressEnOp(Traits &traits, Grid &grid);
                          GStressEnOp(const GStressEnOp &);
                         ~GStressEnOp();

        void              apply(StateComp &d, State &u, GINT idir, State  &utmp, 
                                StateComp &si);                              // stress op evaluation in idir
        void              apply(StateComp &d, State &u, State  &utmp,  
                                StateComp &e);                               // stress-energy op evaluation


private:
        Mass             *massop_;     // mass matrix, required
        Grid             *grid_;       // grid set on construction
        StateComp        *mu_;         // dynamic/shear viscosity
        StateComp        *zeta_;       // bulk/dilitation viscosity
        StateComp        *kappa_;      // dyn.shear visc for energy
        StateComp        *lambda_;     // bulk/diliataion visc for energy
        StateComp         tfact_;      // diliataion truncation factor
        Traits            traits_;     // operator traits

};


#include "gstressen.ipp"


#endif
