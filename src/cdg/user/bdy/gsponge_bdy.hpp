//==================================================================================
// Module       : gsponge_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for 
//                'sponge' boundaries.
//                
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================
#if !defined(_GSPONGE_UPDATE_HPP)
#define _GSPONGE_UPDATE_HPP

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
class GSpongeBdy : public UpdateBdyBase<TypePack>
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

        // GSpongeBdy solver traits:
        struct Traits {
          GINT            idir  = GDIM;  
                                     // canonical coord direction definining surfaces
          GTVector<GINT>  istate;    // state indices to operate on
          GTVector<Ftype> farfield;  // far-field solution for each istate
          GTVector<Ftype> exponent;  // fall-off exponent for solution
          GTVector<Ftype> sigma;     // 'diffusion' factor in sponge layer
          GTVector<Ftype> rs;        // vector defining sponge surface
          GTVector<Ftype> ro;        // vector defining outer-most surface
        };

        GSpongeBdy() = delete; 
        GSpongeBdy(GSpongeBdy<TypePack>::Traits &traits);
       ~GSpongeBdy();
        GSpongeBdy(const GSpongeBdy &bu) = default;
        GSpongeBdy &operator=(const GSpongeBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
        
private:

//      void                init(GSpongeBdy::Traits &);                     // initialize 
        GBOOL               update_cart (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);

        GBOOL               update_sphere(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub);
       

        GBOOL               bcomputed_;     // already computed?
        GBOOL               bcomput_once_;  // compute once??
        Traits              traits_;        // Traits structure

};

#include "gsponge_bdy.ipp"

#endif
