//==================================================================================
// Module       : geoglow.h
// Date         : 7/7/19 (DLR)
// Description  : GeoFLOW main driver header
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#if !defined(GEOFLOW_MAIN_H)
#define GEOFLOW_MAIN_H

#include "gexec.h"
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
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "gposixio_observer.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/equation_factory.hpp"
#include "pdeint/integrator_factory.hpp"
#include "pdeint/mixer_factory.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/null_observer.hpp"
#include "ginitstate_direct_user.hpp"
#include "ginitforce_direct_user.hpp"
#include "ginitbdy_user.hpp"
#include "gupdatebdy_user.hpp"
#include "gspecbdy_user.hpp"
#include "ginitstate_factory.hpp"
#include "ginitforce_factory.hpp"
#include "ginitbdy_factory.hpp"
#include "gupdatebdy_factory.hpp"
#include "gspecbdy_factory.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"

using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;

template< // default template arg types
typename StateType  = GTVector<GTVector<GFTYPE>*>,
typename GridType   = GGrid,
typename ValueType  = GFTYPE,
typename DerivType  = StateType,
typename TimeType   = ValueType,
typename CompType   = GTVector<GStateCompType>,
typename JacoType   = StateType,
typename SizeType   = GSIZET
>
struct EquationTypes {
        using State      = StateType;
        using Grid       = GridType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};

using MyTypes     = EquationTypes<>;           // Define types used
using EqnBase     = EquationBase<MyTypes>;     // Equation Base Type
using EqnBasePtr  = std::shared_ptr<EqnBase>;  // Equation Base Ptr
using MixBase     = MixerBase<MyTypes>;        // Stirring/mixing Base Type
using MixBasePtr  = std::shared_ptr<MixBase>;  // Stirring/mixing Base Ptr
using ObsBase     = ObserverBase<EqnBase>;     // Observer Base Type
using BasisBase   = GTVector<GNBasis<GCTYPE,GFTYPE>*>; // Basis pool type


// 'Member' data:
GBOOL  bench_=FALSE;
GINT   nsolve_;  // number *solved* 'state' arrays
GINT   nstate_;  // number 'state' arrays
GINT   ntmp_  ;  // number tmp arrays
GGrid *grid_;
State  u_;
State  c_;
State  ub_;
State  uf_;
State  utmp_;
GTVector<GFTYPE> nu_(3);
BasisBase        gbasis_(GDIM);
PropertyTree     ptree_;      // main prop tree
GGFX<GFTYPE>     ggfx_;       // DSS operator
GC_COMM          comm_ ;      // communicator


// Callback functions:
void update_boundary  (const Time &t, State &u, State &ub);     // bdy vector update
void update_forcing   (const Time &t, State &u, State &uf);     // forcing vec update
void steptop_callback (const Time &t, State &u, const Time &dt);// backdoor function

// Public methods:
void init_state       (const PropertyTree &ptree, GGrid &, EqnBasePtr &pEqn, Time &t, State &utmp, State &u, State &ub);
void init_force       (const PropertyTree &ptree, GGrid &, EqnBasePtr &pEqn, Time &t, State &utmp, State &u, State &uf);
void init_bdy         (const PropertyTree &ptree, GGrid &, Time &t, State &utmp, State &u, State &ub);
void allocate         (const PropertyTree &ptree);
void deallocate       ();
void create_observers (EqnBasePtr &eqn_ptr, PropertyTree &ptree, GSIZET icycle, Time time,     std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers);
void create_equation  (const PropertyTree &ptree, EqnBasePtr &pEqn);
void create_mixer     (PropertyTree &ptree, MixBasePtr &pMixer);
void create_basis_pool(PropertyTree &ptree, BasisBase &gbasis);
void init_ggfx        (PropertyTree &ptree, GGrid &grid, GGFX<GFTYPE> &ggfx);
void gresetart        (PropertyTree &ptree);
void do_bench         (GString sbench, GSIZET ncyc);

//#include "init_pde.h"

#endif
