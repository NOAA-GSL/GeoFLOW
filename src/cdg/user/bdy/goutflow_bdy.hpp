//==================================================================================
// Module       : g_outflow_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                outflow boundaries
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
#include "ggfx.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GOutflowBdy : public UpdateBdyBase<TypePack>
{
public:
        using Types      = TypePack;
        using Base       = UpdateBdyBase<Types>;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

        // GOutflowBdy solver traits:
        struct Traits {
          GINT           bdyid;    // bdy id
          vector<GINT>   istate;   // state indices to operate on
          vector<GSIZET> ibdyvol;  // indir. inidices into comput volume

        };

        GOutflowBdy() = delete; 
        GOutflowBdy(typename GOutflowBdy<Types>::Traits &traits);
       ~GOutflowBdy();
        GOutflowBdy(const GOutflowBdy &bu) = default;
        GOutflowBdy &operator=(const GOutflowBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);

        std::vector<int>&   get_istate_impl() { return traits_.istate; }


        
private:


        Traits              traits_;        // Traits structure

};


#include "goutflow_bdy.ipp"


#endif
