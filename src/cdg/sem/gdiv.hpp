//==================================================================================
// Module       : gdiv.hpp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the divergence
//                operator,
//                      Div(rho \vec{v})
//                where rho is a scalar field, and  \vec{v} is
//                a vector field. This isn't the strictly conservative
//                form that uses Gauss' theorem to add surfaces fluxes; 
//                it is volume-integrated.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

#if !defined(_GDIVOP_HPP)
#define _GDIVOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"


template<typename TypePack>
class GDivOp
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

        // MConv solver traits:
        struct Traits {
          GBOOL           docollocation = FALSE;   // colocation or weak forms?
        };

                          GDivOp() = delete;
                          GDivOp(Traits &traits, Grid &grid);
                          GDivOp(const GDivOp &);
                         ~GDivOp();

        void              apply(StateComp &d, State &u, State  &utmp, 
                                StateComp &div);                             // stress op evaluation in idir
        void              apply(State &u, State  &utmp,  
                                StateComp &div);                             // stress-energy op evaluation


private:
        Traits                        traits_;    // traits structure
        Mass                         *massop_;    // mass matrix, required
        Grid                         *grid_;      // grid set on construction

};


#include "gdiv.ipp"


#endif
