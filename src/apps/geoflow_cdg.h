//==================================================================================
// Module       : geoflow_cdg.h
// Date         : 7/7/19 (DLR)
// Description  : GeoFLOW main driver for CG and DG initial value,
//                boundary value problems
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#if !defined(_GEOFLOW_CDG_MAIN_H)
#define _GEOFLOW_CDG_MAIN_H

#if defined(_OPENMP)
   #include <omp.h>
#endif

#include "gtypes.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include <typeinfo>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "gburgers.hpp"
#include "gmconv.hpp"
#include "gio.hpp"
#include "gstateinfo.hpp"
//#include "ggrid.hpp"
//#include "ggrid_box.hpp"
//#include "ggrid_icos.hpp"
#include "ggrid_factory.hpp"
#include "ginitstate_factory.hpp"
#include "ginitforce_factory.hpp"
#include "gupdatebdy_factory.hpp"
#include "gspecterrain_factory.hpp"
#include "pdeint/equation_factory.hpp"
#include "pdeint/integrator_factory.hpp"
#include "pdeint/mixer_factory.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/filter_factory.hpp"
#include "tbox/command_line.hpp"
#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/tracer.hpp"


template< // Complete typepack
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateCompType = GTVector<GFTYPE>,
typename StateInfoType = GStateInfo,
typename FloatType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = FloatType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET
>
struct MyTypePack {
        using EqnBase    = EquationBase<MyTypePack>;   // Equation Base type
        using EqnBasePtr = std::shared_ptr<EqnBase>;   // Equation Base ptr
        using State      = StateType;
        using StateComp  = StateCompType;
        using Grid       = GGrid<MyTypePack>;
        using GridBox    = GGridBox<MyTypePack>;
        using GridIcos   = GGridIcos<MyTypePack>;
        using StateInfo  = StateInfoType;
        using Mass       = GMass<MyTypePack>;
        using Ftype      = FloatType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using IBdyVol    = GTVector<GSIZET>;
        using TBdyVol    = GTVector<GBdyType>;
        using Size       = SizeType;
        using FilterBasePtr = std::shared_ptr<FilterBase<MyTypePack>>;
        using FilterList    = std::vector<FilterBasePtr>;
};
using StateInfo     = GStateInfo;
using MyTypes       = MyTypePack<>;               // Define grid types used
using EqnBase       = typename MyTypes::EqnBase;           
using EqnBasePtr    = typename MyTypes::EqnBasePtr;
using Grid          = typename MyTypes::Grid;           
using IOBaseType    = IOBase<MyTypes>;            // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using MixBase       = MixerBase<MyTypes>;         // Stirring/mixing Base type
using MixBasePtr    = std::shared_ptr<MixBase>;   // Stirring/mixing Base ptr
using IntegratorBase= Integrator<MyTypes>;        // Integrator
using IntegratorPtr = std::shared_ptr<IntegratorBase>; // Integrator ptr
                                                  // Integrator ptr
using ObsBase       = ObserverBase<EqnBase>;      // Observer Base type
using ObsTraitsType = ObserverBase<MyTypes>::Traits;
using BasisBase     = GTVector<GNBasis<GCTYPE,GFTYPE>*>; 
                                                  // Basis pool type

// 'Member' data:
GBOOL            bench_=FALSE;
GINT             nsolve_;      // number *solved* 'state' arrays
GINT             nstate_;      // number 'state' arrays
GINT             ntmp_  ;      // number tmp arrays
GINT             irestobs_;    // index in pObservers_ that writes restarts
Grid            *grid_=NULLPTR;// grid object
State            u_;           // state array
State            c_;           // advection velocity, if used
State            uf_;          // forcing tendency
State            utmp_;        // temp array
GTVector<GFTYPE> nu_(3);       // viscosity
BasisBase        gbasis_(GDIM);// basis vector
EqnBasePtr       pEqn_;        // equation pointer
IntegratorPtr    pIntegrator_; // integrator pointer
MixBasePtr       pMixer_;      // mixer object
PropertyTree     ptree_;       // main prop tree

GGFX<GFTYPE>    *ggfx_=NULLPTR;// DSS operator
IOBasePtr        pIO_=NULLPTR; // ptr to IOBase operator
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>>
                 pObservers_(new std::vector<std::shared_ptr<ObserverBase<MyTypes>>>()); // observer array
GC_COMM          comm_ ;       // communicator


// Callback functions:
void update_forcing   (const Time &t, State &u, State &uf);     // forcing vec update
void steptop_callback (const Time &t, State &u, const Time &dt);// backdoor function

// Public methods:
void init_state       (const PropertyTree &ptree, Grid &, EqnBasePtr &pEqn, Time &t, State &utmp, State &u);
void init_force       (const PropertyTree &ptree, Grid &, EqnBasePtr &pEqn, Time &t, State &utmp, State &u, State &uf);
void allocate         (const PropertyTree &ptree);
void deallocate       ();
void create_observers (EqnBasePtr &eqn_ptr, PropertyTree &ptree, GSIZET icycle, Time time, std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers);
void create_equation  (const PropertyTree &ptree, EqnBasePtr &pEqn);
void create_mixer     (PropertyTree &ptree, MixBasePtr &pMixer);
void create_basis_pool(PropertyTree &ptree, BasisBase &gbasis);
void do_terrain       (const PropertyTree &ptree, Grid &grid);
void init_ggfx        (PropertyTree &ptree, Grid &grid, GGFX<GFTYPE> *&ggfx);
void gresetart        (PropertyTree &ptree);
void compare          (const PropertyTree &ptree, Grid &, EqnBasePtr &pEqn, Time &t, State &utmp, State &u);
void do_restart       (const PropertyTree &ptree, Grid &, State &u, GTMatrix<GINT>&p,  GSIZET &cycle, Time &t);
void do_bench         (GString sbench, GSIZET ncyc);

//#include "init_pde.h"

#endif
