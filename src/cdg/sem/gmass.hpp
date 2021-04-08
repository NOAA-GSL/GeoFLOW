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
#include "ggfx.hpp"
#include "pdeint/equation_base.hpp"


using namespace geoflow::pdeint;
using namespace std;


template<typename TypePack>
class GMass
{

public:

        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using Grid       = typename Types::Grid;
        using Ftype      = typename Types::Ftype;


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


        GBOOL             bInitialized_;
        GBOOL             bdoinverse_;
        GBOOL             bmasslumped_;
        GTVector<Ftype>   mass_;
        Grid             *grid_;


};


#include "gmass.ipp"


#endif
