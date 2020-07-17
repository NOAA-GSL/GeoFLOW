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
#include "ggrid.hpp"
#include "ggfx.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GOutflowBdy : public UpdateBdyBase<TypePack>
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

        // GOutflowBdy solver traits:
        struct Traits {
          GBOOL     compute_once=FALSE; // compute bdy cond once?
          GINT             bdyid;    // bdy id
          GTVector<GINT>   istate;   // state indices to operate on
          GTVector<GSIZET> ibdyvol;  // indir. inidices into comput volume

        };

        GOutflowBdy() = delete; 
        GOutflowBdy(GOutflowBdy<TypePack>::Traits &traits);
       ~GOutflowBdy();
        GOutflowBdy(const GOutflowBdy &bu) = default;
        GOutflowBdy &operator=(const GOutflowBdy &bu) = default;


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

#include "goutflow_bdy.ipp"

#endif
