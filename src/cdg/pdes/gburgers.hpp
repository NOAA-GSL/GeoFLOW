//==================================================================================
// Module       : gburgers.hpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE:
//                     du/dt + u . Del u = nu Del^2 u
//                This solver can be built in 2D or 3D, and can be configured to
//                remove the nonlinear terms so as to solve only the heat equation.
//
//                The State variable must always be of specific type
//                   GTVector<GTVector<GFTYPE>*>, but the elements rep-
//                resent different things depending on whether
//                the equation is doing nonlinear advection, heat only, or 
//                pure linear advection. If solving with nonlinear advection or 
//                the heat equation, the State consists of elements [*u1, *u2, ....]
//                which is a vector solution. If solving the pure advection equation,
//                the State consists of [*u, *c1, *c2, *c3], where u is the solution
//                desired (there is only one, and it's a scalar), and ci 
//                are the constant-in-time Eulerian velocity components.
// 
// 
//                The dissipation coefficient may also be provided as a spatial field.
// 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
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
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
//#include "gflux.hpp"


template<typename TypePack>
class GBurgers :: public EquationBase<TypePack>
{
        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");

public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        // Burgers solver traits:
        struct Traits {
          GBOOL        doheat;
          GBOOL        bpureadv;
          GBOOL        bconserved;
          GStepperType steptype;
          GINT         itorder;
          GINT         inorder;
        };

        GBurgers() = delete; 
        GBurgers(GGrid &grid, State &u, Traits &traits, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers();
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const Burgers &bu) = default;

protected:
        GBOOL               has_dt_impl() const {return FALSE;}           // Has dynamic dt?
        GBOOL               step_impl(const Time &t, State &uin, 
                                      Time &dt, State &uout);             // Take a step
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                set_nu(GTVector<GFTYPE> &nu);                 // Set nu
        void                apply_bc_impl(const Time &t, State &u, 
                                          State &ub);                     // Apply bdy conditions
       void                 set_bdy_callback(
                            std::function<void(Time &t, State &u,
                                          State &ub)> &callback)          // set bdy-update callback
                              {bdy_update_callback_ = callback;}

private:

        void                init(State &u, Traits &);      // initialize 
        void                step_exrk  (const Time &t, State &uin,
                                        Time &dt, State &uout);
        void                dudt_impl  (const Time &t, State &u,
                                        Time &dt, State &dudt);
        void                step_multistep(const Time &t, State &uin,
                                           Time &dt, State &uout);
        void                cycle_keep(State &u);
       

        GBOOL               doheat_;        // flag to do heat equation alone
        GBOOL               bpureadv_;      // do pure (linear) advection?
        GBOOL               bconserved_;    // use conservation form?
        GStepperType        isteptype_;     // stepper type
        GINT                nsteps_         // num steps taken
        GINT                itorder_;       // time deriv order
        GINT                inorder_;       // nonlin term order
        GTVector<GFTYPE>    tcoeffs_;       // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;       // coeffs for NL adv term
        GButcherRK          butcher_;       // Butcher tableau for EXRK
        GTVector<GTVector<GFLOAT>*>  
                            tmp_;
        GTVector<GTVector<GFLOAT>*>  
                            c_;             // linear velocity if bpureadv = TRUE
        GTVector<State>     u_keep_;        // state at prev. time levels
        GTVector<GStepperType>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>   *nu_   ;         // dissipoation
        GExRKstepper        gexrk_;         // ExRK stepper, if needed
        GMassop            *gmass_;         // mass op
        GMassop            *gimass_;        // inverse mass op
        GAdvect            *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GpdV               *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op
        std::function<void(Time &t, State &u, State &ub)>
                           *bdy_update_callback_; // bdy update callback function


};

#endif
