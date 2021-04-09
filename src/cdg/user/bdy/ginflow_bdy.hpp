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
#include "ggfx.hpp"
#include "ginitstate_factory.hpp"
#include "pdeint/update_bdy_base.hpp"
#include "tbox/property_tree.hpp"


using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GInflowBdy : public UpdateBdyBase<TypePack>
{
public:
        using Types      = TypePack;
        using Base       = UpdateBdyBase<Types>;
        using EqnBase    = EquationBase<TypePack>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using GridBox    = typename Types::GridBox;
        using GridIcos   = typename Types::GridIcos;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

        // GInflowBdy solver traits:
        struct Traits {
          GBOOL     compute_once=FALSE; // compute bdy cond once?
          GBOOL         use_init=FALSE; // set via state initialzation method?
          GINT                   bdyid; // bdy id
          GTVector<GINT>        istate; // state indices to operate on
          GTVector<GSIZET>     ibdyvol; // indir. inidices into comput volume
          GString              smethod; // init method if use_init==TRUE

          std::function<GBOOL(EqnBasePtr &eqn,
                              Grid       &grid, 
                              Time       &time,
                              const GINT  id,
                              State      &utmp,
                              State      &u,
                              State      &ub)> callback = NULLPTR;
          geoflow::tbox::PropertyTree ptree;
        };

        GInflowBdy() = delete; 
        GInflowBdy(typename GInflowBdy<Types>::Traits &traits);
       ~GInflowBdy();
        GInflowBdy(const GInflowBdy &bu) = default;
        GInflowBdy &operator=(const GInflowBdy &bu) = default;


protected:
        GBOOL               update_impl (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u);
        
private:
        GBOOL               compute_bdy_data (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u,
                              State      &ub);


        Traits              traits_;        // Traits structure
        GBOOL               bcomputed_;     // tells us that bdydata was computed
        State               bdydata_;       // bdy data arrays

};

#include "ginflow_bdy.ipp"

#endif
