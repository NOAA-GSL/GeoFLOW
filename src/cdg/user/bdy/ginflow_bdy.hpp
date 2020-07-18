//==================================================================================
// Module       : ginflow_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                inflow boundaries that are set by initialization method,
//                of by external function call.
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_GINFLOW_BDY_HPP)
#define _GINFLOW_BDY_HPP

#include "gtypes.h"
#include <functional>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ggfx.hpp"
#include "ginitstate_factory.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GInflowBdy : public UpdateBdyBase<TypePack>
{
public:
        using Types      = UpdateBdyBase<TypePack>;
        using Base       = Types;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Value;
        using Time       = typename Types::Time;
        using StateInfo  = typename Types::StateInfo;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // GInflowBdy solver traits:
        struct Traits {
          GBOOL     compute_once=FALSE; // compute bdy cond once?
          GBOOL        buse_init=FALSE; // set via state initialzation method?
          GINT                   bdyid; // bdy id
          GTVector<GINT>        istate; // state indices to operate on
          GTVector<GSIZET>     ibdyvol; // indir. inidices into comput volume

          std::function<void(Grid       &grid, 
                             StateInfo  &stinfo,
                             Time       &time,
                             State      &utmp,
                             State      &u,
                             State      &ub)> callback = NULLPTR;
        };

        GInflowBdy() = delete; 
        GInflowBdy(typename GInflowBdy<Types>::Traits &traits);
       ~GInflowBdy();
        GInflowBdy(const GInflowBdy &bu) = default;
        GInflowBdy &operator=(const GInflowBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
        
private:
        GBOOL               update_from_init (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
        GBOOL               update_from_user (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);


        Traits              traits_;        // Traits structure
        GBOOL               bcomputed_;     // tells us that operation was called
        State               unew_;          // helper vector
        State               tmpnew_;        // helper vector

};

#include "ginflow_bdy.ipp"

#endif
