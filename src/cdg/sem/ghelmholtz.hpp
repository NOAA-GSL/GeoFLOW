//==================================================================================
// Module       : ghelmholtz.hpp
// Date         : 11/1/18 (DLR)
// Description  : Represents the SEM generalized Helmholtz operator:
//                  H = qM + pL,
//                where M = mass operator, L is Laplacian operator, and
//                q, p are scalars that may or may not be constant. 
//                The mass term is added only if calls are made to 
//                set_mass_scalar and p is applied only if set_Lap_scalar, 
//                respectively. In fact, only if the mass operator is set
//                is this operator a real Helmholtz operator; otherwise, it's
//                really just a weak Laplacian operator. 
//                Note: this operator will fail if the grid contains more 
//                      than a single element type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#if !defined(_GHELMHOLTZOP_HPP)
#define _GHELMHOLTZOP_HPP
#include "gtvector.hpp"
#include "gmass.hpp"
#include "gelem_base.hpp"
#include "pdeint/equation_base.hpp"


using namespace geoflow::pdeint;
using namespace std;


template<typename TypePack>
class GHelmholtz 
{

public:
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using Grid       = typename Types::Grid;
        using Mass       = typename Types::Mass;
        using Ftype      = typename Types::Ftype;
        using Derivative = typename Types::Derivative;
        using Size       = typename Types::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
               "StateComp is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Ftype>*>>::value,
               "Derivative is of incorrect type");



                          GHelmholtz(Grid &grid);
                          GHelmholtz(const GHelmholtz &);
                         ~GHelmholtz();

        void              opVec_prod(StateComp &in, 
                                     State     &utmp,
                                     StateComp &out);                         // Operator-vector product 
        void              set_Lap_scalar(GTVector<Ftype> &p);                // Scalar multipliying Laplacian
        void              set_mass_scalar(GTVector<Ftype> &q);               // Scalar multiplying Mass
        void              init();                                             // must call after all 'sets'
        void              use_metric(GBOOL flag) {buse_metric_ = flag;}       // set flag to use metric & Jacobian
//      void              set_tmp(GTVector<GTVector<Ftype>*> &utmp) 
//                        { utmp_.resize(utmp.size()); utmp_ = utmp; }       // Set temp space 

private:
        void              def_init();
        void              reg_init();
        void              def_prod(StateComp &in, 
                                   State     &utmp,
                                   StateComp &out);
        void              reg_prod(StateComp   &in, 
                                   State       &utmp,
                                   StateComp   &out);
        void              embed_prod(StateComp  &in, 
                                     State      &utmp,
                                     StateComp  &out);
        void              compute_refderivs(GTVector<Ftype> &, 
                                            GTVector<GTVector<Ftype>*> &, GBOOL btrans=FALSE);
        void              compute_refderivsW(GTVector<Ftype> &, 
                                             GTVector<GTVector<Ftype>*> &, GBOOL btrans=FALSE);

        void              compute_div(State       &, 
                                      StateComp   &, GBOOL btrans=TRUE); 

        GBOOL                         bInitialized_;  
        GBOOL                         buse_metric_;   // use metric terms?
        GBOOL                         bown_q_;        // does object own q?
        GBOOL                         bown_p_;        // does object own p?
        GBOOL                         bown_mass_;     // does object own massop?
        GBOOL                         bcompute_helm_; // compute full Helm, not just Lap?
        GTVector<Ftype>              *p_;    // scalar multiplying L
        GTVector<Ftype>              *q_;    // scalar multiplying M
        GTVector<Ftype>               etmp1_;// elem-based (not global) tmp vector
        State                         utmp_; // global array of temp vectors
        GTMatrix<GTVector<Ftype>*>   G_;    // metric components
        Grid                         *grid_; // grid set on construction


};


#include "ghelmholtz.ipp"


#endif
