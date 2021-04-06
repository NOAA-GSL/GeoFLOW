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
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "ggfx.hpp"
#include "pdeint/update_bdy_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GSpongeBdy : public UpdateBdyBase<TypePack>
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
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // GSpongeBdy solver traits:
        struct Traits {
           
          GINT             idir  = GDIM;  
                                     // canonical coord direction definining surfaces
          GINT             bdyid;    // bdy id
          GTVector<GSIZET> ibdyvol;  // indir. inidices into comput volume
          GTVector<GINT>   istate;   // state indices to operate on
          GTVector<Ftype>  farfield; // far-field solution for each istate
          GTVector<Ftype>  exponent; // fall-off exponent for solution
          GTVector<Ftype>  sigma;    // 'diffusion' factor in sponge layer
          GTVector<Ftype>  rs;       // vector defining sponge surface
          Ftype            ro;       // outer-most surface (may be negative)
        };

        GSpongeBdy() = delete; 
        GSpongeBdy(typename GSpongeBdy<Types>::Traits &traits);
       ~GSpongeBdy();
        GSpongeBdy(const GSpongeBdy &bu) = default;
        GSpongeBdy &operator=(const GSpongeBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);
        
private:

//      void                init(GSpongeBdy::Traits &);                     // initialize 
        GBOOL               update_cart (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);

        GBOOL               update_sphere(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);
       

        Traits              traits_;        // Traits structure

};

#include "gsponge_bdy.ipp"

#endif
