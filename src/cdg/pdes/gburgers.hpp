//==================================================================================
// Module       : gburgers.hpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GBURGERS_HPP)
#define _GBURGERS_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "equation_base.hpp"

class GBurgers :: public EquationBase<TypePack>
{

public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Grid       = typename Interface::Grid;

        GBurgers() = delete; 
//      GBurgers(GGrid &grid, GTvector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers() = default;
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const Burgers &bu) = default;


//friend  std::ostream&       operator<<(std::ostream&, GBurgers &) {};    // Output stream operator

protected:
        GBOOL               has_dt_impl() const {return FALSE;}          // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);  // Get dt
        void                dudt_impl(const Time &t, State& u, 
                                      Derivative& dudt);                 // Compute RHS
        void                dfdu_impl(const Time& t, State &u, 
                                      Jacobian &dfdu);                   // Compute Jacobian dF/du
        void                set_bdy_callback(
                            std::function<void(GGrid &)> &callback);     // set bdy-set callback

private:
        void                init2d();                                    // initialize for 2d grid
        void                init3d();                                    // initialize for 3d grid
       

std::function<void(GGrid&)>
                       *bdycallback_ ; // callback object+method to set bdy conditions
State                   tmp_;

};

#endif
