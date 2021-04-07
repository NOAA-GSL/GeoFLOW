//==================================================================================
// Module       : gburgers.hpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE:
//                     du/dt + c(u) . Del u = nu Del^2 u
//                This solver can be built in 2D or 3D, and can be configured to
//                remove the nonlinear terms so as to solve only the heat equation.
//
//                The State variable must always be of specific type
//                   GTVector<GTVector<GFTYPE>*>, but the elements rep-
//                resent different things depending on whether
//                the equation is doing nonlinear advection, heat only, or 
//                pure linear advection. If solving with nonlinear advection 
//                equation, the State consists of elements [*u1, *u2, ....]
//                which is a vector solution. If solving the pure advection equation,
//                the State consists of [*u, *c1, *c2, *c3], where u is the solution
//                desired (there is only one, and it's a scalar), and ci 
//                are the unevolved Eulerian velocity components. If solving the 
//                heat equation, the state consists of [*u] alone. PRESCRIBED quantities
//                in state vector, such as adv. velocities, must be placed at the end
//                of the state vector.
// 
// 
//                The dissipation coefficient may also be provided as a spatial field.
// 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : EquationBase.
//==================================================================================
#if !defined(_GBURGERS_HPP)
#define _GBURGERS_HPP

#include "gtypes.h"
#include <functional>
#include <cstdlib>
#include <memory>
#include <cmath>
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
        using Types      = TypePack;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Ftype;
        using Derivative = typename Types::Derivative;
        using Time       = typename Types::Time;
        using CompDesc   = typename Types::CompDesc;
        using Jacobian   = typename Types::Jacobian;
        using Size       = typename Types::Size;
        using FilterBasePtr = std::shared_ptr<FilterBase<Types>>;
        using FilterList    = std::vector<FilterBasePtr>;


        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // Burgers solver traits:
        struct Traits {
          GBOOL          doheat      = FALSE;
          GBOOL          bpureadv    = FALSE;
          GBOOL          bconserved  = FALSE;
          GBOOL          bforced     = FALSE;
          GBOOL          variabledt  = FALSE;
          GBOOL          bSSP        = FALSE;// use strong stab preserv. RK?
          GINT           nstate      = GDIM; // no. vars in state vec
          GINT           nsolve      = GDIM; // no. vars to solve for
          GINT           ntmp        = 8;
          GINT           itorder     = 2;
          GINT           nstage      = 2;
          GINT           inorder     = 2;
          GFTYPE         courant     = 0.5;
          GFTYPE         nu          = 0.0;
          GTVector<GINT> iforced;
          GString        ssteptype;
        };

        GBurgers() = delete; 
        GBurgers(Grid &grid, GBurgers<TypePack>::Traits &traits);
       ~GBurgers();
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const GBurgers &bu) = default;


        void                set_steptop_callback(
                            std::function<void(const Time &t, State &u, 
                                               const Time &dt)> callback) 
                             { steptop_callback_ = callback; bsteptop_ = TRUE;}
        GBurgers<TypePack>::Traits &get_traits() { return traits_; }      // Get traits
                                            

protected:
        void                step_impl(const Time &t, State &uin, State &uf, 
                                      const Time &dt);                    // Take a step
        void                step_impl(const Time &t, const State &uin, State &uf, 
                                      const Time &dt, State &uout);       // Take a step
        void                compute_derived_impl(const State &uin, GString sop, 
                                      State &utmp, State &uout, 
                                      std::vector<GINT> &iuout);          // Compute derived quantity
        GINT                tmp_size_impl();                              // required tmp size
        GINT                solve_size_impl()                             // required solve size
                            {return traits_.nsolve;}
        GINT                state_size_impl()                             // required state size
                            {return traits_.nstate;}
        std::vector<GINT>  &iforced_impl()                                // required force comps
                            {return stdiforced_;}

        void                init_impl(State &u, State &tmppool);          // initialize 
        GTVector<GFTYPE>    &get_nu() { return nu_; };                    // Set nu/viscosity
        GBOOL               has_dt_impl() const {return bvariabledt_;}    // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                apply_bc_impl(const Time &t, State &u);       // Apply bdy conditions

private:

        void                dudt_impl  (const Time &t, const State &u, const State &uf, 
                                        const Time &dt, Derivative &dudt);
        void                step_exrk  (const Time &t, State &uin, State &uf, 
                                        const Time &dt, State &uout);
        void                step_multistep(const Time &t, State &uin, State &uf, 
                                           const Time &dt);
        void                cycle_keep(State &u);
       

        GBOOL               bInit_;         // solver initialized?
        GBOOL               doheat_;        // flag to do heat equation alone
        GBOOL               bpureadv_;      // do pure (linear) advection?
        GBOOL               bconserved_;    // use conservation form?
        GBOOL               bforced_;       // use forcing vectors
        GBOOL               bsteptop_;      // is there a top-of-step callback?
        GBOOL               bvariabledt_;   // is dt allowed to vary?
        GBOOL               bSSP_;          // use strong stab. preserv. RK?
        GStepperType        isteptype_;     // stepper type
        GINT                nsteps_ ;       // num steps taken
        GINT                itorder_;       // time deriv order
        GINT                inorder_;       // nonlin term order
        GINT                nstage_;        // no. stages in time deriv RK
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
        std::vector<GINT>   stdiforced_;    // std::verctor for traits_.iforced
        GGrid              *grid_;          // GGrid object
        GExRKStepper<GFTYPE>
                           *gexrk_;         // ExRK stepper, if needed
        GMass              *gmass_;         // mass op
        GMass              *gimass_;        // inverse mass op
        GAdvect<TypePack>  *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GpdV<TypePack>     *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op
        GC_COMM             comm_;          // communicator
        GGFX<GFTYPE>       *ggfx_;          // gather-scatter operator
        GBurgers<TypePack>::Traits         
                            traits_;        // solver traits

        std::function<void(const Time &t, State &u, const Time &dt)>
                           steptop_callback_;


};

#include "gburgers.ipp"

#endif
