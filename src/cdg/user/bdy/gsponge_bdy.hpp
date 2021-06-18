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
        using GridBox    = typename Types::GridBox;
        using GridIcos   = typename Types::GridIcos;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

        // GSpongeBdy solver traits:
        struct Traits {
           
          GINT             idir  = GDIM;  
                                     // canonical coord direction definining surfaces
          GINT             bdyid;    // bdy id
          GTVector<GSIZET> ibdyvol;  // indir. inidices into comput volume
          GTVector<GINT>   istate;   // state indices to operate on
          GTVector<Ftype>  farfield; // far-field solution for each istate
          GTVector<Ftype>  falloff;  // fall-off rate for solution
          GTVector<Ftype>  exponent; // decay exponents 
          vector<GSIZET>   isponge;  // contains grid indices of sponge layer points; 

          Ftype            xstart;   // number defining sponge surface start
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

        void                init(Grid &grid);                     // initialize 
        GBOOL               update_box(
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
       

        GBOOL               binit_;         // object initialized?
        GFTYPE              xmax_;          // max coord in direction idir
        Traits              traits_;        // Traits structure
        GTVector<GSIZET>    isponge_;       // contains grid indices of sponge layer points; 

};

#include "gsponge_bdy.ipp"

#endif
