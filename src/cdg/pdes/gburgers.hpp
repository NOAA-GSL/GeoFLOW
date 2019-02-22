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
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gtvector.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gab.hpp"
#include "gext.hpp"
#include "gbdf.hpp"
#include "gpdv.hpp"
#include "gmass.hpp"
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
#include "gbc.hpp"
//#include "gflux.hpp"
#include "gexrk_stepper.hpp"
#include "gbutcherrk.hpp"
#include "ggfx.hpp"
#include "pdeint/equation_base.hpp"

using namespace geoflow::pdeint;
using namespace std;




template<typename TypePack>
class GBurgers : public EquationBase<TypePack>
{
public:
        using Interface  = EquationBase<TypePack>;
        using Base       = Interface;
        using State      = typename Interface::State;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");

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
        GBurgers(GGFX &ggfx, GGrid &grid, State &u, GBurgers<TypePack>::Traits &traits, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers();
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const GBurgers &bu) = default;

        void                 set_bdy_update_callback(
                             std::function<void(const Time &t, State &u,
                                           State &ub)> *callback) 
                             { update_bdy_callback_ = callback;
                               gbc_->set_update_callback(callback); }     // set bdy-update callback


protected:
        void                step_impl(const Time &t, State &uin, State &ub, 
                                      const Time &dt){};                  // Take a time step
        void                step_impl(const Time &t, const State &uin, State &ub,
                                      const Time &dt, State &uout);       // Take a step
        GBOOL               has_dt_impl() const {return FALSE;}           // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                set_nu(GTVector<GFTYPE> &nu);                 // Set nu
        void                apply_bc_impl(const Time &t, State &u, 
                                          const State &ub);               // Apply bdy conditions
private:

        void                init(State &u, GBurgers::Traits &);           // initialize 
        GINT                req_tmp_size();                               // required tmp size
        void                step_exrk  (const Time &t, const State &uin, State &ub,
                                        const Time &dt, State &uout);
        void                dudt_impl  (const Time &t, const State &u,
                                        const Time &dt, Derivative &dudt);
        void                step_multistep(const Time &t, const State &uin, State &ub,
                                           const Time &dt, State &uout);
        void                cycle_keep(State &u);
       

        GBOOL               doheat_;        // flag to do heat equation alone
        GBOOL               bpureadv_;      // do pure (linear) advection?
        GBOOL               bconserved_;    // use conservation form?
        GStepperType        isteptype_;     // stepper type
        GINT                nsteps_ ;       // num steps taken
        GINT                itorder_;       // time deriv order
        GINT                inorder_;       // nonlin term order
        GTVector<GFTYPE>    tcoeffs_;       // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;       // coeffs for NL adv term
        GTVector<GFTYPE>    dthist_;        // coeffs for NL adv term
        GTVector<GTVector<GFTYPE>*>  
                            utmp_;
        GTVector<GTVector<GFTYPE>*>  
                            urhstmp_;       // helper arrays set from utmp
        GTVector<GTVector<GFTYPE>*>  
                            urktmp_;        // helper arrays set from utmp
        GTVector<GTVector<GFTYPE>*>  
                            c_;             // linear velocity if bpureadv = TRUE
        GTVector<State>     ukeep_;         // state at prev. time levels
        GTVector<GStepperType>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>   *nu_   ;         // dissipoation
        GGrid              *grid_;          // GGrid object
        GExRKStepper<GFTYPE>
                           *gexrk_;         // ExRK stepper, if needed
        GMass              *gmass_;         // mass op
        GMass              *gimass_;        // inverse mass op
        GAdvect            *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GpdV               *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op
        GBC                *gbc_;           // bdy conditions operator
        GC_COMM             comm_;          // communicator
        GGFX               *ggfx_;          // gather-scatter operator
        std::function<void(const Time &t, State &u,
                           State &ub)> *update_bdy_callback_;


};

#include "gburgers.ipp"

#endif
