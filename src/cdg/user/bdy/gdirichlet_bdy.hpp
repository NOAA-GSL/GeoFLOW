//==================================================================================
// Module       : gdirichlet_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                Dirichlet boundary conditions
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_GDIRICHLET_BDY_HPP)
#define _GDIRICHLET_BDY_HPP

#include "gtypes.h"
#include <functional>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ggfx.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GDirichletBdy : public UpdateBdyBase<TypePack>
{
public:
        using Interface  = UpdateBdyBaseBase<TypePack>;
        using Base       = Interface;
        using State      = typename Interface::State;
        using Grid       = typename Interface::Grid;
        using Ftype      = typename Interface::Value;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Ftype>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // GDirichletBdy solver traits:
        struct Traits {
          GBOOL    compute_once=TRUE;     // compute bdy condition once?
          GINT             bdyid;    // bdy id
          GTVector<FType>  value;    // Diriclet value for each istate 
          GTVector<GINT>  istate;    // state indices to operate on
        };

        GDirichletBdy() = delete; 
        GDirichletBdy(GDirichletBdy<TypePack>::Traits &traits);
       ~GDirichletBdy();
        GDirichletBdy(const GDirichletBdy &bu) = default;
        GDirichletBdy &operator=(const GDirichletBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
        
private:


        Traits              traits_;        // Traits structure
        GBOOL               bcomputed_;     // tell us that operation was called

};

#include "gdirichlet_bdy.ipp"

#endif
