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
                           std::function<void(const GFTYPE &t, GTVector<GTVector<T>*> &u,
                                         GTVector<GTVector<T>*> &ub)> *callback)
                                          { bdy_apply_callback_ = callback; }   // set bdy-application callback

        void               set_update_bdy_callback(
                           std::function<void(const GFTYPE &t, GTVector<GTVector<T>*> &u,
                                         GTVector<GTVector<T>*> &ub)> *callback)
                                          { bdy_update_callback_ = callback; }   // set bdy-update callback

        void               setRHSfunction(std::function<void(
                                          const T &t, 
                                          const GTVector<GTVector<T>*> &uin,
                                          const T &dt, 
                                          GTVector<GTVector<T>*> &dudt)> *callback)
                                          { rhs_callback_ = callback; }          // RHS callback

        void               step(const T &t, const GTVector<GTVector<T>*> &uin,
                                GTVector<GTVector<T>*> &ub,
                                const T &dt, GTVector<GTVector<T>*> &tmp,
                                GTVector<GTVector<T>*> &uout);


private:
// Private data:
        GINT               nstage_;                          // no stages (not nec. 'order'!)
        GButcherRK<T>      butcher_;                         // Butcher tableau
        std::function<void(const GFTYPE &t,                  // RHS callback function
                           const GTVector<GTVector<GFTYPE>*> &uin,
                           const GFTYPE &dt, 
                           GTVector<GTVector<GFTYPE>*> &dudt)>
                           *rhs_callback_;
        std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u,
        GTVector<GTVector<GFTYPE>*> &ub)>
                           *bdy_update_callback_;            // bdy update callback
        std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u,
        GTVector<GTVector<GFTYPE>*> &ub)>
                           *bdy_apply_callback_;             // bdy apply callback

};

#include "gexrk_stepper.ipp"

#endif

