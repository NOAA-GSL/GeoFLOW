//==================================================================================
// Module       : gmconv.hpp
// Date         : 6/11/20 (DLR)
// Description  : Object defining a moist convection solver:
// 
//                PDEs:
//                  d_t rhoT + Div (rhoT v)   = -Ltot, total mass
//                  d_t q_v + v.Grad q_v      = q_v Ltot/rhoT + dot(s_v)/rhoT
//                  d_t q_h + v.Grad q_i      = q_h Ltot/rhoT - Div(rhoT q_h W_i)/rhoT
//                                            + dot(s_h)/rhoT
//
//                where 
//                  Ltot = Sum_h Div(rhoT q_h vec{W}_i), 
//                is the total mass loss due to hydrometeor fallout, and              
//                and q_h are the hydrometeor (liquid and ice) mass
//                fractions. The dry mass fraction is:
//                  q_d = 1 - Sum_h q_h.
//                Note:
//                  Sum_i dot(s_i) = 0, where sum is over all densities,
//                and dot(s_i) are the mass sources for vapor, and hydrometeors.
//                This solver can be built in 2D or 3D for box or Icos grids, 
//                but for the 2D embedded sphere, no 'preferred directions'
//                are allowed (gravity, fallout).
//
//                The State vector consists of the following:
//                  vx/sx   : x-velocity or momentum density 
//                  vy/sy   : y-velocity or momentum density 
//                  vz/sz   : z-velocity or momentum density 
//                  e       : internal or total energy density
//                  rho_tot : total density or dry density, if no moisture
//                  qvapor  : water vapor mass fraction
//                  q_liq_0 : liquid substance 0 mass fraction |
//                  q_liq_1 : liquid substance 1 mass fraction |  'liquid' mass sector
//                  q_liq_2 : liquid substance 2 mass fraction |
//                   ...
//                  q_ice_0 : 'ice' substance 0 mass fraction  |
//                  q_ice_1 : 'ice' substance 1 mass fraction  |  'ice' mass sector
//                  q_ice_2 : 'ice' substance 2 mass fraction  |
//                   ...
//                  rho_s   : base state density               |
//                  p_s     : base state pressure              | base state sector, prescibed
//
//                  w_liq_0 : liquid substance 2 term velocity |
//                  w_liq_1 : liquid substance 2 term velocity | 'liquid' term vel. sector,
//                  w_liq_2 : liquid substance 2 term velocity |  prescribed

//                  w_ice_0 : 'ice' substance 2 term velocity  |
//                  w_ice_1 : liquid substance 2 term velocity | 'ice' term vel. sector,
//                  w_ice_2 : 'ice' substance 2 term velocity  |  prescribed
//                   ...
//
//                The base state & terminal velocities in this state are prescribed,
//                if used; all others are evolved. If hydrometeor fallout
//                is specified, then terminal velocities for all hydrometeors
//                must be provided. Each vector may be either of length 1, 
//                for which it's assumed to be a constant, or it may be 
//                space-dependent, for which it provides the terminal 
//                velocity in the preferred diection at all spatial locaions.
//                PRESCRIBED quantities in state vector, such as base state 
//                or terminal velocities, must be placed at at the end of the 
//                state vector.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : EquationBase.
//==================================================================================
#if !defined(_GMCONV_HPP)
#define _GMCONV_HPP

#include "gtypes.h"
#include "gphysics.h"
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
#include "gstressen.hpp"
#include "gdiv.hpp"
//#include "gflux.hpp"
#include "gexrk_stepper.hpp"
#include "gbutcherrk.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "ggfx.hpp"
#include "gutils.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"


using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GMConv : public EquationBase<TypePack>
{
public:
        using Interface  = EquationBase<TypePack>;
        using Base       = Interface;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Ftype      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;
        using FilterList = typename Interface::FilterList;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<GFTYPE>>::value,
               "StatCompe is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // MConv solver traits:
        struct Traits {
          GBOOL           dodry       = TRUE;   // do dry dynamics?
          GBOOL           docoriolis  = FALSE;  // use Coriolis force?
          GBOOL           dograv      = TRUE;   // use gravitational force?
          GBOOL           dofallout   = FALSE;  // allow precip fallout?
          GBOOL           bconserved  = FALSE;  // use conserved form?
          GBOOL           bforced     = FALSE;  // use forcing?
          GBOOL           usemomden   = TRUE;   // use momentum density form?
          GBOOL           usebase     = TRUE;   // use hydrostatic base state?
          GBOOL           variabledt  = FALSE;  // use variable timestep?
          GBOOL           bvarvterm   = FALSE;  // time dep term velocities?
          GBOOL           divopcolloc = TRUE;   // use collocation in GDivOp?
          GBOOL           Stokeshyp   = FALSE;  // use Stokes hypothesis
          GBOOL           bindepdiss  = FALSE;  // indep. mom & energy diss?
          GBOOL           bSSP        = FALSE;  // use strong stab pres RK?
          GINT            nstate      = GDIM+2; // no. vars in state vec
          GINT            nsolve      = GDIM+2; // no. vars to solve for
          GINT            nlsector    = 0;      // no. vars in liq-sector
          GINT            nisector    = 0;      // no. vars in ice-sector
          GINT            nbase       = 2;      // no. vars in base state (p, d)
          GINT            itorder     = 2;      // formal  time iorder
          GINT            nstage      = 2;      // no. stages for time integ.
          GINT            inorder     = 2;      // formal nonlin. extrap order
          GStepperType    isteptype   = GSTEPPER_EXRK;
          GFTYPE          Ts_base     = 300.0;  // base state surf temp (K)
          GFTYPE          P0_base     = 1000.0; // base state ref pressure (mb)
          GFTYPE          courant     = 0.5;    // Courant factor
          GFTYPE          nu          = 0.0;    // shear viscosity constant
          GFTYPE          kappa       = 0.0;    // energy-shear visc constant
          GFTYPE          zeta        = 0.0;    // mom bulk viscosity constant
          GFTYPE          lambda      = 0.0;    // energy bulk shear visc const
          GTVector<GINT>  iforced;              // state comps to force
          GTVector<Ftype> omega;                // rotation rate vector
          GString         ssteptype;            // stepping method
        };

        GMConv() = delete; 
        GMConv(Grid &grid, GMConv<TypePack>::Traits &traits);
       ~GMConv();
        GMConv(const GMConv &bu) = default;
        GMConv &operator=(const GMConv &bu) = default;

        StateComp           &get_nu() { return nu_; };                       // Set nu/viscosity
        StateComp           &get_kappa() { return kappa_; };                 // Set nu/viscosity

        void                set_steptop_callback(
                            std::function<void(const Time &t, State &u, 
                                               const Time &dt)> callback) 
                             { steptop_callback_ = callback; bsteptop_ = TRUE;}

        GMConv<TypePack>::Traits &get_traits() { return traits_; }           // Get traits
                                            

protected:
        void                step_impl(const Time &t, State &uin, State &uf, State &ub, 
                                      const Time &dt);                    // Take a step
        void                step_impl(const Time &t, const State &uin, State &uf, State &ub,
                                      const Time &dt, State &uout);       // Take a step
        GBOOL               has_dt_impl() const {return traits_.variabledt;}  
        GINT                tmp_size_impl();                              // required tmp size
        GINT                solve_size_impl()                             // required solve size
                            {return traits_.nsolve;}
        GINT                state_size_impl()                             // required tmp size
                            {return traits_.nstate;}
        std::vector<GINT>  &iforced_impl()                                // required tmp size
                            {return stdiforced_;}
        void                init_impl(State &u, State &tmppool);                    // initialize 
                                                                          // Has dynamic dt?
        void                compute_derived_impl(const State &uin, GString sop,
                                      State &utmp, State &uout,
                                      std::vector<GINT> &iuout);          // Compu

        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                apply_bc_impl(const Time &t, State &u, 
                                          State &ub);                     // Apply bdy conditions
private:

        void                dudt_impl  (const Time &t, const State &u, const State &uf, const State &ub,
                                        const Time &dt, Derivative &dudt);
        void                step_exrk  (const Time &t, State &uin, State &uf, State &ub,
                                        const Time &dt, State &uout);
        void                step_multistep(const Time &t, State &uin, State &uf, State &ub,
                                           const Time &dt);
        void                cycle_keep   (const State &u);
inline  void                compute_cv   (const State &u, StateComp &utmp, StateComp &cv);
inline  void                compute_qd   (const State &u, StateComp &qd);
inline  void                compute_falloutsrc
                                         (StateComp &g, State &qi, State &v, GINT jexcl, State &utmp, StateComp &r );
inline  void                compute_v    (const State &u, StateComp &id, State &v);
inline  void                compute_vpref(StateComp &tv, State &W);
inline  void                compute_vpref(StateComp &tv, GINT idir, StateComp &W);
inline  void                assign_helpers(const State &u, const State &uf);
inline  void                compute_pe   (StateComp &rhoT, State &qi, State &tvi, State &utmp, StateComp      &r);
        void                compute_base  (State &u);
inline  GINT                szrhstmp();
 

        GBOOL               bInit_;         // object initialized?
        GBOOL               bforced_;       // use forcing vectors
        GBOOL               bsteptop_;      // is there a top-of-step callback?
        GBOOL               bvterm_;        // teminal vel. computed?
        GINT                istage_;        // RK stage number
        GINT                nevolve_;       // num StateComp's evolved
        GINT                nhydro_;        // num hydrometeors
        GINT                nmoist_;        // number of moist components
        GSIZET              icycle_;        // internal cycle number
        GStepperType        isteptype_;     // stepper type
        GTVector<GFTYPE>    tcoeffs_;       // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;       // coeffs for NL adv term
        GTVector<GFTYPE>    dthist_;        // coeffs for NL adv term
        State               uold_;          // helper arrays set from utmp
        State               uevolve_;       // helper array to specify evolved state components
        State               ubase_;         // helper array pointing to base state components
        State               utmp_;          // tmp pool
        State               urhstmp_;       // helper arrays set from utmp
        State               urktmp_;        // helper arrays set from utmp
        State               qi_;            // full mass fraction vector
        State               qice_;          // ice mass fraction vector
        State               qliq_;          // liquid mass fraction vector
        State               tvi_;           // term vel. vector for all qi
        State               tvice_;         // term vel. vector for all qice
        State               tvliq_;         // term vel. vector for all qliq
        State               fv_;            // state velocity forcing components
        State               s_;             // state momentum components
        State               v_;             // state velocity components
        State               W_;             // terminal velocity components
        StateComp           dtmp_;          // density base state temp array
        StateComp           ptmp_;          // pressure base state temp array
        GTVector<State>     ukeep_;         // state at prev. time levels
        GTVector<GString>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>    nu_   ;         // KE dissipoation
        GTVector<GFTYPE>    kappa_;         // internal energy dissipoation
        GTVector<GFTYPE>    dxmin_ ;        // element face mins
        GTVector<GFTYPE>    maxbyelem_ ;    // element-based maxima for dt
        std::vector<GINT>   stdiforced_;    // traits_.iforced as a std::vector
        GGrid              *grid_;          // GGrid object
        GExRKStepper<GFTYPE>
                           *gexrk_;         // ExRK stepper, if needed
        GMass              *gmass_;         // mass op
        GMass              *gimass_;        // inverse mass op
        GAdvect<TypePack>  *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GStressEnOp<TypePack>
                           *gstressen_;     // viscous stress-energy
        GpdV<TypePack>     *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op
        GDivOp<TypePack>   *gdiv_;          // volumetric divergence op
        GC_COMM             comm_;          // communicator
        GGFX<GFTYPE>       *ggfx_;          // gather-scatter operator
        GMConv<TypePack>::Traits 
                            traits_;        // traits structure

        std::function<void(const Time &t, State &u, const Time &dt)>
                           steptop_callback_;

        // Starting indices for each sector:
        GINT               MOMENTUM;
        GINT               ENERGY;
        GINT               DENSITY;
        GINT               VAPOR;
        GINT               LIQMASS;
        GINT               ICEMASS;
        GINT               PRESCRIBED;
        GINT               BASESTATE;
        GINT               LIQTERMV;
        GINT               ICETERMV;


};

#include "gmconv.ipp"

#endif
