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
#include "ggfx.hpp"
#include "ggrid.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GDirichletBdy : public UpdateBdyBase<TypePack>
{
public:
        using Types      = TypePack;
        using Base       = UpdateBdyBase<Types>;
        using EqnBase    = EquationBase<TypePack>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Value;
        using Time       = typename Types::Time;

//      static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
//             "State is of incorrect type");
//      static_assert(std::is_same<Grid,GGrid>::value,
//             "Grid is of incorrect type");

        // GDirichletBdy solver traits:
        struct Traits {
          GINT             bdyid;    // bdy id
          GTVector<Ftype>  value;    // Diriclet value for each istate 
          GTVector<GINT>   istate;   // state indices to operate on
          GTVector<GSIZET> ibdyvol;  // indir. inidices into comput volume

        };

        GDirichletBdy() = delete; 
        GDirichletBdy(typename GDirichletBdy<Types>::Traits &traits);
       ~GDirichletBdy();
        GDirichletBdy(const GDirichletBdy &bu) = default;
        GDirichletBdy &operator=(const GDirichletBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);
        
private:


        Traits              traits_;        // Traits structure

};

#include "gdirichlet_bdy.ipp"

#endif
