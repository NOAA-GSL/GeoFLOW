//==================================================================================
// Module       : gsimple_outflow_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                simple outflow boundaries
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_GSIMPLE_OUTFLOW_BDY_HPP)
#define _GSIMPLE_OUTFLOW_BDY_HPP

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
class GSimpleOutflowBdy : public UpdateBdyBase<TypePack>
{
public:
        using Interface  = UpdateBdyBaseBase<TypePack>;
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

        // GSimpleOutflowBdy solver traits:
        struct Traits {
          GBOOL     compute_once=FALSE; // compute bdy cond once?
          GTVector<GINT>  istate;    // state indices to operate on
        };

        GSimpleOutflowBdy() = delete; 
        GSimpleOutflowBdy(GSimpleOutflowBdy<TypePack>::Traits &traits);
       ~GSimpleOutflowBdy();
        GSimpleOutflowBdy(const GSimpleOutflowBdy &bu) = default;
        GSimpleOutflowBdy &operator=(const GSimpleOutflowBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
        
private:

        GBOOL               bcomputed_;     // was computation done?

        Traits              traits_;        // Traits structure

};

#include "gsimple_outflow_bdy.ipp"

#endif
