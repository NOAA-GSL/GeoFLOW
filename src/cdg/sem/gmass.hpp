//==================================================================================
// Module       : gmass.hpp
// Date         : 10/19/18 (DLR)
// Description  : Represents the SEM mass operator.  Mass-lumping is assumed.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(_GMASSOP_HPP)
#define _GMASSOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "pdeint/equation_base.hpp"


template<typename TypePack>
class GMass: 
{

public:

        using Types      = EquationBase<TypePack>;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using Grid       = typename Types::Grid;
        using Mass       = typename Types::Mass;
        using Ftype      = typename Types::Ftype;
        using Derivative = typename Types::Derivative;
        using Time       = typename Types::Time;
        using CompDesc   = typename Types::CompDesc;
        using Jacobian   = typename Types::Jacobian;
        using Size       = typename Types::Size;


                          GMass(Grid &grid, GBOOL doinverse=FALSE);
                          GMass(const GMass &);
                         ~GMass();

        void              opVec_prod(GTVector<Ftype> &in, 
                                     GTVector<GTVector<Ftype>*> &utmp,
                                     GTVector<Ftype> &out);                       // Operator-vector product
        GTVector<Ftype>  *data() { return &mass_; }
//      void              do_mass_lumping(GBOOL bml);                              // Set mass lumping flag

private:
        void              init();
        void              init1d();
        void              init2d();
        void              init3d();


        GBOOL             bdoinverse_;
        GBOOL             bmasslumped_;
        GTVector<Ftype>   mass_;
        Grid             *grid_;


};


#include "gmass.ipp"


#endif
