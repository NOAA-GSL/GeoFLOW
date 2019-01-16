//==================================================================================
// Module       : gburgers_equation.hpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GBURGERS_EQUATION_HPP)
#define _GBURGERS_EQUATION_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "equation_base.hpp"
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
//#include "gflux.hpp"

class GBurgers_equation :: public EquationBase<TypePack>
{
        static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVectorGFTYPE>>>::value,
               "Derivative is of incorrect type");

public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        GBurgers_equation() = delete; 
        GBurgers_equation(GGrid &grid, State &u, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers_equation();
        GBurgers_equation(const GBurgers_equation &bu) = default;
        GBurgers_equation &operator=(const Burgers &bu) = default;


//friend  std::ostream&       operator<<(std::ostream&, GBurgers_equation &) {};    // Output stream operator

protected:
        GBOOL               has_dt_impl() const {return FALSE;}          // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);  // Get dt
        void                dudt_impl(const Time &t, State& u, 
                                      Derivative& dudt);                 // Compute RHS
        void                dfdu_impl(const Time& t, State &u, 
                                      Jacobian &dfdu);                   // Compute Jacobian dF/du
        void                apply_bc_impl();                             // Apply bdy conditions

private:

        GBOOL               bconserved_;

        void                init();                                    // initialize 
        void                init2d();                                    // initialize for 2d grid
        void                init3d();                                    // initialize for 3d grid
       

        GTVector<GTVector<GFLOAT>*>  
                            tmp_;
        GMassop            *gmass_;
        GAdvect            *gadvect_;
        GHelmholtz         *ghelm_;
        GpdV               *gpdv_;
//      GFlux              *gflux_;

};

#endif
