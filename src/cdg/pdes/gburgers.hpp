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
// Derived From : EquationBase.
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
        using Grid       = typename Interface::Grid;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");
        GTVector<GStateCompType> 
                         icomptype_; // descriptions of each state component

        // Burgers solver traits:
        struct Traits {
          GBOOL          doheat      = FALSE;
          GBOOL          bpureadv    = FALSE;
          GBOOL          bconserved  = FALSE;
          GBOOL          bforced     = FALSE;
          GBOOL          variabledt  = FALSE;
          GINT           nstate      = GDIM; // no. vars in state vec
          GINT           nsolve      = GDIM; // no. vars to solve for
          GINT           ntmp        = 8;
          GINT           itorder     = 2;
          GINT           inorder     = 2;
          GFTYPE         courant     = 0.5;
          GTVector<GFTYPE> 
                         nu          = 0.0;
          GTVector<GINT> iforced;
          GString        ssteptype   = "GSTEPPER_EXRK";
        };

        GBurgers() = delete; 
        GBurgers(GGFX<GFTYPE> &ggfx, Grid &grid, GBurgers<TypePack>::Traits &traits, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers();
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const GBurgers &bu) = default;

        GTVector<GFTYPE>    &get_nu() { return nu_; };                       // Set nu/viscosity
        void                 set_bdy_update_callback(
                             std::function<void(const Time &t, State &u,
                                           State &ub)> callback) 
                             { update_bdy_callback_ = callback; bupdatebc_ = TRUE;
                               if ( gbc_ != NULLPTR ) gbc_->set_update_callback(callback); 
                               if ( gexrk_ != NULLPTR ) 
                                 gexrk_->set_update_bdy_callback(callback);} // set bdy-update callback

        void                set_steptop_callback(
                            std::function<void(const Time &t, State &u, 
                                               const Time &dt)> callback) 
                             { steptop_callback_ = callback; bsteptop_ = TRUE;}
                                            

protected:
        void                step_impl(const Time &t, State &uin, State &uf, State &ub, 
                                      const Time &dt);                    // Take a step
        void                step_impl(const Time &t, const State &uin, State &uf, State &ub,
                                      const Time &dt, State &uout);       // Take a step
        GBOOL               has_dt_impl() const {return bvariabledt_;}    // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                apply_bc_impl(const Time &t, State &u, 
                                          const State &ub);               // Apply bdy conditions
private:

        void                init(GBurgers::Traits &);                     // initialize 
        GINT                req_tmp_size();                               // required tmp size
        void                dudt_impl  (const Time &t, const State &u, const State &uf, const State &ub,
                                        const Time &dt, Derivative &dudt);
        void                step_exrk  (const Time &t, State &uin, State &uf, State &ub,
                                        const Time &dt, State &uout);
        void                step_multistep(const Time &t, State &uin, State &uf, State &ub,
                                           const Time &dt);
        void                cycle_keep(State &u);
       

        GBOOL               doheat_;        // flag to do heat equation alone
        GBOOL               bpureadv_;      // do pure (linear) advection?
        GBOOL               bconserved_;    // use conservation form?
        GBOOL               bforced_;       // use forcing vectors
        GBOOL               bupdatebc_;     // bdy update callback set?
        GBOOL               bsteptop_;      // is there a top-of-step callback?
        GBOOL               bvariabledt_;   // is dt allowed to vary?
        GStepperType        isteptype_;     // stepper type
        GINT                nsteps_ ;       // num steps taken
        GINT                itorder_;       // time deriv order
        GINT                inorder_;       // nonlin term order
        GFTYPE              courant_;       // Courant number if dt varies
        GTVector<GFTYPE>    tcoeffs_;       // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;       // coeffs for NL adv term
        GTVector<GFTYPE>    dthist_;        // coeffs for NL adv term
        GTVector<GTVector<GFTYPE>*>  
                            uevolve_;       // helper array to specify evolved sstate components
        State               utmp_;
        State               uold_;          // helper arrays set from utmp
        State               urhstmp_;       // helper arrays set from utmp
        State               uoptmp_;        // helper arrays set from utmp
        State               urktmp_;        // helper arrays set from utmp
        State               c_;             // linear velocity if bpureadv = TRUE
        GTVector<State>     ukeep_;         // state at prev. time levels
        GTVector<GString>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>    nu_   ;         // dissipoation
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
        GGFX<GFTYPE>       *ggfx_;          // gather-scatter operator
        std::function<void(const Time &t, State &u,
                           State &ub)> update_bdy_callback_;
        std::function<void(const Time &t, State &u, const Time &dt)>
                           steptop_callback_;


};

#include "gburgers.ipp"

#endif
