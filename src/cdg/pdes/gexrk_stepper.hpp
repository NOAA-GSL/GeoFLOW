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

#define GTVector<GTVector<T>*> State


template <typename T>
class GExRKStepper
{
public:
                           GExRKStepper() = delete;
                           GExRKStepper(GSIZET nstage);
                          ~GExRKStepper();
                           GExRKStepper(const GExRKStepper &a) = default;
                           GExRKStepper &operator=(const GExRKStepper &bu) = default;
        void               setOrder(GINT nstage);

        void               set_apply_bdy_callback(
                           std::function<void(const T &t, State &u,
                                         State &ub)> *callback)
                                         { bdy_apply_callback_ = callback; }   // set bdy-application callback

        void               set_update_bdy_callback(
                           std::function<void(const T &t, State &u,
                                         GTVector<GTVector<T>*> &ub)> *callback)
                                         { bdy_update_callback_ = callback; }   // set bdy-update callback

        void               setRHSfunction(std::function<void(
                                          const T &t, 
                                          const State &uin,
                                          const T &dt, 
                                          State &dudt)> *callback)
                                          { rhs_callback_ = callback; }          // RHS callback

        void               step(const T &t, const State &uin,
                                State &ub,
                                const T &dt, State &tmp,
                                State &uout);


private:
// Private data:
        GINT               nstage_;                          // no stages (not nec. 'order'!)
        GButcherRK<T>      butcher_;                         // Butcher tableau
        std::function<void(const T &t,                  // RHS callback function
                           const State  &uin,
                           const T &dt, 
                           State &dudt)>
                           *rhs_callback_;
        std::function<void(const T &t, State &u, State &ub)>
                           *bdy_update_callback_;            // bdy update callback
        std::function<void(const T &t, State &u, State &ub)>
                           *bdy_apply_callback_;             // bdy apply callback

};

#include "gexrk_stepper.ipp"

#endif

