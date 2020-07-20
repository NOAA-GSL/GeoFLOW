//==================================================================================
// Module       : gexrk_stepper.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGExRKSTEPPER_HPP)
#define _GGExRKSTEPPER_HPP

#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gbutcherrk.hpp"
#include "ggrid.hpp"
#include "ggfx.hpp"

#include <functional>

template <typename T>
class GExRKStepper
{
typedef GTVector<GTVector<T>*> State;
typedef GTVector<T> StateElem;
typedef T  Time;

public:
                           GExRKStepper() = delete;
                           GExRKStepper(GGrid &grid, GSIZET nstage);
                          ~GExRKStepper();
                           GExRKStepper(const GExRKStepper &a) = default;
                           GExRKStepper &operator=(const GExRKStepper &bu) = default;
        void               setOrder(GINT nstage);

        void               step(const Time &t, const State &uin,
                                State &uf,
                                State &ub,
                                const Time &dt, State &tmp,
                                State &uout);

        void               step(const Time &t, State &uin, 
                                State &uf,
                                State &ub,
                                const Time &dt, State &tmp);

        void               setRHSfunction(std::function<void(
                                          const Time &t, 
                                          const State &uin,
                                          const State &uf,
                                          const State &ub,
                                          const Time &dt, 
                                          State &dudt)> callback)
                                          { rhs_callback_ = callback; 
                                            bRHS_ = TRUE; }           // RHS callback, required

        void               set_apply_bdy_callback(
                           std::function<void(const Time &t, State &u,
                                         State &ub)> callback)
                                         { bdy_apply_callback_ = callback;
                                           bapplybc_ = TRUE; }        // set bdy-application callback
        void                set_ggfx(GGFX<GFTYPE> *ggfx)
                            {ggfx_ = ggfx;}                           // set geom-free exchange op



private:
// Private methods:
        void               resize(GINT nstate);              // resize member data 

// Private data:
        GBOOL              bRHS_;
        GBOOL              bapplybc_;
        GINT               nstage_;                          // no stages (not nec. 'order'!)
        GButcherRK<T>      butcher_;                         // Butcher tableau
        GTVector<State>    K_;                               // RK stage update vectors
        GGrid             *grid_;                            // grid object
        GGFX<GFTYPE>      *ggfx_;                            // geom-free exchange op
        std::function<void(const Time &t,                    
                           const State  &uin,
                           const State  &uf,
                           const State  &ub,
                           const Time &dt, 
                           State &dudt)>
                           rhs_callback_;                   // RHS callback function
        std::function<void(const Time &t, State &u, State &ub)>
                           bdy_apply_callback_;             // bdy apply callback

};

#include "gexrk_stepper.ipp"

#endif

