//==================================================================================
// Module       : gexrk_stepper.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GGExRKSTEPPER_HPP)
#define _GGExRKSTEPPER_HPP

#include <functional>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gbutcherrk.hpp"
#include "ggfx.hpp"
#include "gmtk.hpp"


template <typename Grid, typename T>
class GExRKStepper
{
         using State     = typename Grid::State;
         using StateComp = typename Grid::StateComp;
         using Ftype     = typename Grid::Ftype;
         using Time      = typename Grid::Time;

public:
        // MConv solver traits:
        struct Traits {
          GBOOL           bSSP        = FALSE;  // do strong stability-pres?
          GINT            norder      = 2;      // order
          GINT            nstage      = 2;      // no. stages
        };
                           GExRKStepper() = delete;
                           GExRKStepper(Traits &traits, Grid &grid);
                           GExRKStepper(Traits &traits);
                          ~GExRKStepper();
                           GExRKStepper(const GExRKStepper &a) = default;
                           GExRKStepper &operator=(const GExRKStepper &bu) = default;
        void               setOrder(GINT order, GINT nstage) {
                             norder_ = order; nstage_=nstage; 
                             if ( !bSSP_ ) butcher_ .setOrder(norder_); }

        void               step(const Time &t, const State &uin,
                                State &uf,
                                const Time &dt, State &tmp,
                                State &uout);

        void               step(const Time &t, State &uin, 
                                State &uf,
                                const Time &dt, State &tmp);

        void               setRHSfunction(std::function<void(
                                          const Time &t, 
                                          const State &uin,
                                          const State &uf,
                                          const Time &dt, 
                                          State &dudt)> callback)
                                          { rhs_callback_ = callback; 
                                            bRHS_ = TRUE; }           // RHS callback, required

        void               set_apply_bdy_callback(
                           std::function<void(const Time &t, State &u
                                             )> callback)
                                         { bdy_apply_callback_ = callback;
                                           bapplybc_ = TRUE; }        // set bdy-application callback
        void                set_ggfx(GGFX<Ftype> *ggfx)
                            {ggfx_ = ggfx;}                           // set geom-free exchange op



private:
// Private methods:
        void               resize(GINT nstate);              // resize member data 
        void               step_b(const Time &t, const State &uin, State &uf,
                                  const Time &dt, State &tmp,
                                  State &uout);

        void               step_b(const Time &t, State &uin, State &uf, 
                                  const Time &dt, State &tmp);  // Butcher-form

        void               step_ssp(const Time &t, const State &uin, State &uf, 
                                    const Time &dt, State &tmp,
                                    State &uout);                 // SSP-form
        void               step_ssp22(const Time &t, const State &uin,
                                    State &uf, 
                                    const Time &dt, State &tmp,
                                    State &uout);                 // SSP-form
        void               step_ssp33(const Time &t, const State &uin,
                                    State &uf, 
                                    const Time &dt, State &tmp,
                                    State &uout);                 // SSP-form
        void               step_ssp34(const Time &t, const State &uin,
                                    State &uf, 
                                    const Time &dt, State &tmp,
                                    State &uout);                 // SSP-form

        void               step_ssp(const Time &t, State &uin, 
                                    State &uf, 
                                    const Time &dt, State &tmp);

        void               step_euler(const Time &t, const State &uin, 
                                      State &uf, 
                                      const Time &dt, State &uout);

// Private data:
        GBOOL              bRHS_;
        GBOOL              bapplybc_;
        GBOOL              bSSP_;                            // is strong-stability-preserving?
        GINT               norder_;                          // order
        GINT               nstage_;                          // no stages (not nec. 'order'!)
        GButcherRK<T>      butcher_;                         // Butcher tableau
        GTVector<State>    K_;                               // RK stage update vectors
        Grid              *grid_;                            // grid object
        GGFX<Ftype>       *ggfx_;                            // geom-free exchange op
        std::function<void(const Time &t,                    
                           const State  &uin,
                           const State  &uf,
                           const Time &dt, 
                           State &dudt)>
                           rhs_callback_;                   // RHS callback function
        std::function<void(const Time &t, State &u)>
                           bdy_apply_callback_;             // bdy apply callback

};

#include "gexrk_stepper.ipp"

#endif

