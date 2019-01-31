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
#include "stepper_base.hpp"
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
//#include "gflux.hpp"

enum GICOSPTYPE {GICOS_BASE, GICOS_ELEMENTAL}; 


class GBurgers :: public StepperBase<TypePack>
{
        static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVectorGFTYPE>>>::value,
               "Derivative is of incorrect type");

public:
        using Interface  = StepperBase<TypePack>;
        using State      = typename Interface::State;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        GBurgers() = delete; 
        GBurgers(GGrid &grid, State &u, GStepperType isteptype, GTVector<GTVector<GFTYPE>*> &tmp);
       ~GBurgers();
        GBurgers(const GBurgers &bu) = default;
        GBurgers &operator=(const Burgers &bu) = default;

protected:
        GBOOL               has_dt_impl() const {return FALSE;}           // Has dynamic dt?
        GBOOL               step_impl(const Time &t, State &uin, 
                                      Time &dt, State &uout);             // Take a step
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                apply_bc_impl();                              // Apply bdy conditions

private:

        GBOOL               bconserved_;
        GStepperType        isteptype_;                                  // stepper type
        GINT                nsteps_                                      // num steps taken
        GINT                itorder_;                                    // time deriv order
        GINT                inorder_;                                    // nonlin term order
        GTVector<GFTYPE>    tcoeffs_;                                    // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;                                    // coeffs for NL adv term
        GButcherRK          butcher_;                                    // Butcher tableau for EXRK
        void                init(State &u);                              // initialize 
        void                step_exrk  (const Time &t, State &uin,
                                        Time &dt, State &uout);
        void                dudt_impl  (const Time &t, State &u,
                                        Time &dt, State &dudt);
        void                step_extbdf(const Time &t, State &uin,
                                        Time &dt, State &uout);
       

        GTVector<GTVector<GFLOAT>*>  
                            tmp_;
        GTVector<State *>   u_keep_;        // state at prev. time levels
        GTVector<GStepperType>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>    butcher_alpha_; // Butcher tableau alpha (time) coeffs
        GTMatrix<GFTYPE>    butcher_beta_;  // Butcher tableau beta coeffs
        GTVector<GFTYPE>    butcher_c_;     // Butcher tableau c coeffs
        GMassop            *gmass_;         // mass op
        GAdvect            *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GpdV               *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op

};

#endif
