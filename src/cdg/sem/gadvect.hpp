//==================================================================================
// Module       : gadvect.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the advection operator:
//                u.Grad p  This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

#if !defined(_GADVECTOP_HPP)
#define _GADVECTOP_HPP
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "gmass.hpp"
#include "gmtk.hpp"


using namespace std;


template<typename TypePack>
class GAdvect
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


                          GAdvect(Grid &grid);
                          GAdvect(const GAdvect &);
                         ~GAdvect();

        void              apply(StateComp  &p, const State &u, 
                                State &utmp, StateComp &po);                       // Operator-field evaluation
        void              init();                                            // must call after all 'sets'

private:
        Grid                         *grid_;   // grid set on construction


};

#include "gadvect.ipp"


#endif
