//==================================================================================
// Module       : gpdv.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the 'pdV' operator:
//                p Div u. This is a nonlinear operator, and requires that 
//                grid consist of elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

#if !defined(_GPDVOP_HPP)
#define _GPDVOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"


template<typename TypePack>
class GpdV 
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
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

public:

                          GpdV() = delete;
                          GpdV(Grid &grid);
                          GpdV(const GpdV &);
                         ~GpdV();

        void              apply(StateComp &p, State &u, State  &utmp, 
                                StateComp &po);                              // Operator-field evaluation

private:
        Mass                         *massop_; // mass matrix, required
        Grid                         *grid_;   // grid set on construction


};


#include "gpdv.ipp"


#endif
