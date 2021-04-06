//==================================================================================
// Module       : gnoslip_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                no-slip boundary conditions. Acts on kinetic
//                vector.
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_GNOSLIP_BDY_HPP)
#define _GNOSLIP_BDY_HPP

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
class GNoSlipBdy : public UpdateBdyBase<TypePack>
{
public:
        using Types      = TypePack;
        using Base       = UpdateBdyBase<Types>;
        using EqnBase    = EquationBase<TypePack>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // GNoSlipBdy solver traits:
        struct Traits {
          GINT             bdyid;    // bdy id
          GTVector<GINT>   istate;   // state indices to operate on
          GTVector<GSIZET> ibdyvol;  // indir. inidices into comput volume

        };

        GNoSlipBdy() = delete; 
        GNoSlipBdy(typename GNoSlipBdy<Types>::Traits &traits);
       ~GNoSlipBdy();
        GNoSlipBdy(const GNoSlipBdy &bu) = default;
        GNoSlipBdy &operator=(const GNoSlipBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);
        
private:


        Traits              traits_;        // Traits structure
        GINT                nstate_;        // size of traits_.istate

};

#include "gnoslip_bdy.ipp"

#endif
