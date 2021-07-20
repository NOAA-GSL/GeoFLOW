//==================================================================================
// Module       : g0flux_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                0-flux boundary conditions. Acts on kinetic
//                vector.
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_G0FLUX_BDY_HPP)
#define _G0FLUX_BDY_HPP

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
class G0FluxBdy : public UpdateBdyBase<TypePack>
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
        using VVecFtype  = GTVector<GTVector<Ftype>>;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

        // G0FluxBdy solver traits:
        struct Traits {
          GINT           bdyid;    // bdy id
          vector<GINT>   istate;   // state indices to operate on
          vector<GSIZET> ibdyvol;  // indir. inidices into comput volume
          vector<GSIZET> ibdyloc;  // inidices of ubdyvol in global bdy array, grid->igbdy
          vector<GUINT>  ibdydsc;  // bdy descriptor
        };

        G0FluxBdy() = delete; 
        G0FluxBdy(typename G0FluxBdy<Types>::Traits &traits);
       ~G0FluxBdy();
        G0FluxBdy(const G0FluxBdy &bu) = default;
        G0FluxBdy &operator=(const G0FluxBdy &bu) = default;


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

#include "g0flux_bdy.ipp"

#endif
