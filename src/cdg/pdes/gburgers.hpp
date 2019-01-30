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
        GBurgers(GGrid &grid, State &u, GTVector<GTVector<GFTYPE>*> &tmp);
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

        void                init();                                    // initialize 
        void                init2d();                                    // initialize for 2d grid
        void                init3d();                                    // initialize for 3d grid
        void                step_rk23(const Time &t, State &uin,
                                      Time &dt, State &uout);
        void                step_abbdf(const Time &t, State &uin,
                                      Time &dt, State &uout);
       

        GTVector<GTVector<GFLOAT>*>  
                            tmp_;
        GTVector<GStepperType>
                            valid_types_;  
        GMassop            *gmass_;
        GAdvect            *gadvect_;
        GHelmholtz         *ghelm_;
        GpdV               *gpdv_;
//      GFlux              *gflux_;

};

#endif
