//==================================================================================
// Module       : gtest_mass.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of GMass classes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

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

struct TypePack {
        using State      = GTVector<GTVector<GFTYPE>*>;
        using StateComp  = GTVector<GFTYPE>;
        using StateInfo  = GStateInfo;
        using Grid       = GGrid<TypePack>;
        using GridBox    = GGridBox<TypePack>;
        using GridIcos   = GGridIcos<TypePack>;
        using Mass       = GMass<TypePack>;
        using Ftype      = GFTYPE;
        using Derivative = State;
        using Time       = Ftype;
        using CompDesc   = GTVector<GStateCompType>;
        using Jacobian   = State;
        using Size       = GSIZET;
        using EqnBase    = EquationBase<TypePack>;     // Equation Base type
        using EqnBasePtr = std::shared_ptr<EqnBase>;   // Equation Base ptr
        using IBdyVol    = GTVector<GSIZET>;
        using TBdyVol    = GTVector<GBdyType>;
        using Operator   = GHelmholtz<TypePack>;
        using GElemList  = GTVector<GElem_base*>;
        using Preconditioner = GHelmholtz<TypePack>;
        using ConnectivityOp = GGFX<Ftype>;
        using FilterBasePtr  = std::shared_ptr<FilterBase<TypePack>>;
        using FilterList     = std::vector<FilterBasePtr>;

};
using Types         = TypePack;                   // Define types used
using EqnBase       = TypePack::EqnBase;          // Equation Base type
using EqnBasePtr    = TypePack::EqnBasePtr;       // Equation Base ptr
using IOBaseType    = IOBase<Types>;              // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using Grid          = TypePack::Grid;
using Mass          = TypePack::Mass;
using Ftype         = TypePack::Ftype;
using Time          = TypePack::Time;


Grid        *grid_;
GC_COMM      comm_=GC_COMM_WORLD ;      // communicator

GINT szMatCache_ = _G_MAT_CACHE_SIZE;
GINT szVecCache_ = _G_VEC_CACHE_SIZE;


void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<Ftype> &ggfx);


int main(int argc, char **argv)
{
    GString   serr ="main: ";
    GINT      errcode, gerrcode;
    Ftype     err, radiusi;
    IOBasePtr pIO;         // ptr to IOBase operator

    typename ObserverBase<Types>::Traits
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

    radiusi = gtree.getValue<Ftype>("radius",1.0);

    // Create basis:
    GTVector<GNBasis<GCTYPE,Ftype>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,Ftype>(pstd[k]);
    }

    EH_MESSAGE("main: Generate grid...");

    ObserverFactory<Types>::get_traits(ptree, "gio_observer", binobstraits);
    grid_ = GGridFactory<Types>::build(ptree, gbasis, pIO, binobstraits, comm_);

    EH_MESSAGE("main: Initialize gather-scatter...");

    // Initialize gather/scatter operator:
    GGFX<Ftype> ggfx;
    init_ggfx(ptree, *grid_, ggfx);

    // Test some of the coord transformation methods:
    Ftype xlat, xlatc, xlong, xlongc;
    Ftype eps = std::numeric_limits<Ftype>::epsilon();

    GTVector<Ftype> f(grid_->ndof());
    GTVector<Ftype> g(grid_->ndof());

    GTVector<Ftype> imult_ref(grid_->ndof());
    GTVector<Ftype> *imult = &imult_ref;
    ggfx.get_imult(imult_ref);

//  GTVector<Ftype> *jac = &(grid_->Jac());
    GTVector<GTVector<Ftype>*> utmp(1);
    Mass massop(*grid_,FALSE);

    f = 1.0;

    EH_MESSAGE("main: Compute integral...");

    Ftype integral;
    Ftype gintegral;

#if 0
    massop.opVec_prod(f,utmp,g);
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
    integral = g.sum();
    GComm::Allreduce(&integral, &gintegral, 1, T2GCDatatype<Ftype>() , GC_OP_SUM, comm_);
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


    Ftype aintegral;
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
void init_ggfx(PropertyTree& ptree, Grid& grid, GGFX<Ftype>& ggfx)
{
  const auto ndof = grid_->ndof();
  std::vector<std::array<Ftype,GDIM>> xyz(ndof);
  for(std::size_t i = 0; i < ndof; i++){
	  for(std::size_t d = 0; d < GDIM; d++){
		  xyz[i][d] = grid.xNodes()[d][i];
	  }
  }

  // Create GGFX
  ggfx.init(0.25*grid.minnodedist(), xyz);

} // end method init_ggfx

