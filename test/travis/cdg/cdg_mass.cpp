//==================================================================================
// Module       : gtest_mass.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of GMass classes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "pdeint/observer_base.hpp"
#include "pdeint/io_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/null_observer.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/error_handler.hpp"
#include "gcomm.hpp"
#include "gmass.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "gio_observer.hpp"
#include "gio.hpp"


using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;

struct stStateInfo {
  GINT        sttype  = 1;       // state type index (grid=0 or state=1)
  GINT        gtype   = 0;       // check src/cdg/include/gtypes.h
  GSIZET      index   = 0;       // time index
  GSIZET      nelems  = 0;       // num elems
  GSIZET      cycle   = 0;       // continuous time cycle
  GFTYPE      time    = 0.0;     // state time
  std::vector<GString>
              svars;             // names of state members
  GTVector<GStateCompType>
              icomptype;         // encoding of state component types    
  GTMatrix<GINT>
              porder;            // if ivers=0, is 1 X GDIM; else nelems X GDIM;
  GString     idir;              // input directory
  GString     odir;              // output directory
};

template< // default template arg types
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateInfoType = stStateInfo,
typename GridType      = GGrid,
typename ValueType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = ValueType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET 
>
struct EquationTypes {
        using State      = StateType;
        using StateInfo  = StateInfoType;
        using Grid       = GridType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};
using MyTypes       = EquationTypes<>;           // Define types used
using EqnBase       = EquationBase<MyTypes>;     // Equation Base type
using EqnBasePtr    = std::shared_ptr<EqnBase>;  // Equation Base ptr
using IOBaseType    = IOBase<MyTypes>;           // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using Grid          = GGrid;


GGrid       *grid_;
GC_COMM      comm_=GC_COMM_WORLD ;      // communicator

void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<GFTYPE> &ggfx);


int main(int argc, char **argv)
{
    GString   serr ="main: ";
    GINT      errcode, gerrcode;
    GFTYPE    err, radiusi;
    IOBasePtr pIO;         // ptr to IOBase operator

    typename ObserverBase<MyTypes>::Traits
                   binobstraits;

    EH_MESSAGE("main: Starting...");

    if ( argc > 1 ) {
      cout << "No arguments accepted" << endl;
      exit(1);
    }
    errcode = 0;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    mpixx::environment  env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv);
    GlobalManager::startup();



    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    EH_MESSAGE("main: Read prop tree...");

    // Get minimal property tree:
    PropertyTree gtree, ptree; 
    GString sgrid;
    std::vector<GINT> pstd(GDIM);

    ptree.load_file("mass_input.jsn");
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gtree       = ptree.getPropertyTree(sgrid);

    assert(sgrid == "grid_icos" && "Must use ICOS grid for now");

    radiusi = gtree.getValue<GFTYPE>("radius",1.0);

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }

    EH_MESSAGE("main: Generate grid...");

    ObserverFactory<MyTypes>::get_traits(ptree, "gio_observer", binobstraits);
    grid_ = GGridFactory<MyTypes>::build(ptree, gbasis, pIO, binobstraits, comm_);

    EH_MESSAGE("main: Initialize gather-scatter...");

    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(ptree, *grid_, ggfx);

    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();

    GTVector<GFTYPE> f(grid_->ndof());
    GTVector<GFTYPE> g(grid_->ndof());

    GTVector<GFTYPE> imult_ref(grid_->ndof());
    GTVector<GFTYPE> *imult = &imult_ref;
    ggfx.get_imult(imult_ref);

//  GTVector<GFTYPE> *jac = &(grid_->Jac());
    GTVector<GTVector<GFTYPE>*> utmp(1);
    GMass massop(*grid_,FALSE);

    f = 1.0;

    EH_MESSAGE("main: Compute integral...");

    GFTYPE integral;
    GFTYPE gintegral;

#if 0
    massop.opVec_prod(f,utmp,g);
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
    integral = g.sum();
    GComm::Allreduce(&integral, &gintegral, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
#else
    gintegral = grid_->integrate(f, g);
#endif

    EH_MESSAGE("main: Write to file...");

    std::ifstream itst;
    std::ofstream ios;
    itst.open("mass_err.txt");
    ios.open("mass_err.txt",std::ios_base::app);

    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
      ios << "#elems" << "  ";
      for ( GSIZET j=0; j<GDIM; j++ ) ios << "p" << j+1 << "  ";
      ios <<  "ndof    area_computed    area_analytic    diff " << std::endl;
    }
    itst.close();

    // Get global no elements and dof:
    GTVector<GSIZET> lsz(2), gsz(2);
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);


    GFTYPE aintegral;
    #if defined(_G_IS2D)
    aintegral = 4.0*PI*pow(radiusi,2.0);
    #elif defined(_G_IS3D)
    aintegral = 4.0*PI*pow(radiusi,3.0)/3.0;
    #endif
    std::cout << "main: integral=" << gintegral << "; area=" << aintegral << std::endl;

    err = fabs(gintegral-aintegral)/aintegral;

    ios <<   gsz[0] << "  ";
    for ( GSIZET j=0; j<GDIM; j++ ) ios << pstd[j] << "  ";
    ios << gsz[1]     << "  "
        << gintegral  << "  " 
        << aintegral  << "  "
        << fabs(gintegral-aintegral) << std::endl;
    ios.close();
    
    errcode = err < 1e-10 ? 0 : 1;

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);

 
    if ( gerrcode != 0 ) {
      cout << serr << " Error: code=" << errcode << endl;
      cout << serr << " Error: " << err << endl;
    }
    else {
      cout << serr << " Success!" << endl;
    }

    GComm::TermComm();

    return(gerrcode);

} // end, main


//==================================================================================
// Module       : gtools.cpp
// Date         : 5/5/19 (DLR)
// Description  : GeoFLOW initialization tools
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
//#include "gtools.h"

//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : Grid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(PropertyTree& ptree, GGrid& grid, GGFX<GFTYPE>& ggfx)
{
  const auto ndof = grid_->ndof();
  std::vector<std::array<GFTYPE,GDIM>> xyz(ndof);
  for(std::size_t i = 0; i < ndof; i++){
	  for(std::size_t d = 0; d < GDIM; d++){
		  xyz[i][d] = grid.xNodes()[d][i];
	  }
  }

  // Create GGFX
  ggfx.init(0.25*grid.minnodedist(), xyz);

} // end method init_ggfx

