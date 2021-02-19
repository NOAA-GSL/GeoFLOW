//==================================================================================
// Module       : gtest_div.cpp
// Date         : 2/19/21 (DLR)
// Description  : GeoFLOW test for GDiv operator
//                improvement
// Copyright    : Copyright 2021. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================


#include "gtypes.h"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gdiv.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/tracer.hpp"

#include <cassert>
//#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <unistd.h>


using namespace geoflow::tbox;
using namespace std;



template< // Complete typepack
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateCompType = GTVector<GFTYPE>,
typename StateInfoType = GStateInfo,
typename GridType      = GGrid,
typename MassOpType    = GMass,
typename ValueType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = ValueType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET
>
struct TypePack {
        using State      = StateType;
        using StateComp  = StateCompType;
        using StateInfo  = StateInfoType;
        using Grid       = GridType;
        using Mass       = MassOpType;
        using Ftype      = ValueType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};
using MyTypes       = TypePack<>;           // Define grid types used
using Grid          = GGrid;
using IOBaseType    = IOBase<MyTypes>;          // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using ObsTraitsType = ObserverBase<MyTypes>::Traits;

Grid      *grid_ = NULLPTR;
IOBasePtr  pIO_  = NULLPTR; // ptr to IOBase operator




int main(int argc, char **argv)
{
	GEOFLOW_TRACE_INITIALIZE();

    GINT    errcode=0 ;
    GINT    nc=GDIM; // no. coords
    GFTYPE  eps=100.0*std::numeric_limits<GFTYPE>::epsilon();
    GFTYPE  dnorm, time=0.0;
    GString sgrid;// name of JSON grid object to use
    GTVector<GINT>
            pvec;
    std::vector<GINT> 
            pstd(GDIM);  
    GTVector<GFTYPE> *mass;
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();
    pio::initialize(myrank);
    GEOFLOW_TRACE();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;     // main prop tree
    PropertyTree gridptree; // grid prop tree
    PropertyTree divptree;  // GDiv prop tree

    GDivOp<MyTypes> *gdiv;
    typename GDivOp<MyTypes>::Traits trdiv;

    /////////////////////////////////////////////////////////////////
    /////////////////////// Initialize system ///////////////////////
    /////////////////////////////////////////////////////////////////
    ptree.load_file("gdiv_test.jsn");    // main param file structure

    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    divptree    = ptree.getPropertyTree("div_test");
    nc          = GDIM;
    assert(sgrid == "grid_box"); // use only Cartesian grids
    assert(pstd.size() >= GDIM); 

    pvec.resize(pstd.size()); pvec = pstd; pvec.range(0,GDIM-1);

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( auto k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }
    
    // Create grid:
    ObsTraitsType binobstraits;
    grid_ = GGridFactory<MyTypes>::build(ptree, gbasis, pIO_, binobstraits, comm);
    mass = grid_->massop().data();

    /////////////////////////////////////////////////////////////////
    ////////////////////// Allocate arrays //////////////////////////
    /////////////////////////////////////////////////////////////////

    // Create state and tmp space:
    State     utmp (4);
    State     v    (nc);
    StateComp da, diff, div, u;
    
    for ( auto j=0; j<utmp .size(); j++ ) utmp [j] = new StateComp(grid_->size());
    for ( auto j=0; j<v    .size(); j++ ) v    [j] = new StateComp(grid_->size());
    div  .resize(grid_->size());
    diff .resize(grid_->size());
    da   .resize(grid_->size());

    /////////////////////////////////////////////////////////////////
    ////////////////////// Compute solutions/////////////////////////
    /////////////////////////////////////////////////////////////////

    // Initialize u: set p, q, r exponents;
    //   u = x^p y^q z^r:
    std::vector<std::vector<GFTYPE>> 
            vpoly = divptree.getArray2D<GFTYPE>("poly");
    std::vector<std::vector<GFTYPE>> 
            pi(GDIM);
    trdiv.docollocation = divptree.getValue<GBOOL>("docolloc");
    gdiv = new GDivOp<MyTypes>(trdiv, *grid_);
    
    GINT    ncyc  = divptree.getValue<GINT>("ncycles",100);
    GINT    idir  = divptree.getValue<GINT>("idir",2);
    GFTYPE  x, y, z;
    GTVector<GTVector<GFTYPE>> 
           *xnodes = &grid_->xNodes();   

    assert(vpoly.size() >= GDIM); 
    assert(idir >= 1 && idir <= GDIM);
    for ( auto j=0; j<GDIM; j++) {
      pi[j].resize(3); 
      for ( auto i=0; i<3   ; i++) pi[j][i] = 0.0;
      for ( auto i=0; i<GDIM; i++) pi[j][i] = vpoly[j][i];
    }

    // Set scalar field, and analytic derivative:
    for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      z = GDIM == 3 ? (*xnodes)[2][j] : 1.0;
      (*v[0])[j]   = pow(x,pi[0][0])*pow(y,pi[0][1])*pow(z,pi[0][2]); 
      (*v[1])[j]   = pow(x,pi[1][0])*pow(y,pi[1][1])*pow(z,pi[1][2]); 
      if ( GDIM == 3 )
      (*v[2])[j]   = pow(x,pi[2][0])*pow(y,pi[2][1])*pow(z,pi[2][2]); 
      da     [j]   = pi[0][0]*pow(x,pi[0][0]-1)*pow(y,pi[0][1])*pow(z,pi[0][2])
                   + pi[1][1]*pow(x,pi[1][0])*pow(y,pi[1][1]-1)*pow(z,pi[1][2]);
      if ( GDIM == 3 ) da[j] +=
                   + pi[2][2]*pow(x,pi[1][0])*pow(y,pi[1][1])*pow(z,pi[1][2]-1);

    } // end, loop over grid points

    for ( auto j=0; j<nc; j++ ) dnorm = da.amax();

    // Compute numerical divergence:
    GEOFLOW_TRACE_START("gdiv");
    for ( auto n=0; n<ncyc; n++ ) {
       gdiv->apply(v, utmp, div);
    }
    GEOFLOW_TRACE_STOP();

    // divide by MJ
    for ( auto j=0; j<div.size(); j++ ) div[j] /= (*mass)[j];

cout << "da  =" <<  da << endl;
cout << "div =" <<  div<< endl;

    // Find inf-norm and L2-norm errors for each method::
    GTVector<GFTYPE> errs(2); // for each method, Linf and L2 errs
    StateComp        lnorm(2), gnorm(2);

    /////////////////////////////////////////////////////////////////
    //////////////////////// Compute Errors /////////////////////////
    /////////////////////////////////////////////////////////////////
    diff     = da - div;
   *utmp [0] = diff;                   // for inf-norm
   *utmp [1] = diff; utmp[1]->rpow(2); // for L2 norm

    lnorm[0] = utmp[0]->infnorm (); 
    gnorm[1] = sqrt(grid_->integrate(*utmp[1],*utmp[2]))/dnorm;

    // Accumulate to find global inf-norm:
    GComm::Allreduce(lnorm.data(), gnorm.data(), 1, T2GCDatatype<GFTYPE>(), GC_OP_MAX, comm);
    errs[0] = gnorm[0]; errs[1] = gnorm[1];
    if ( myrank == 0 ) {
      if ( errs[1] > eps ) {
        std::cout << "main: ---------------------------GDivOp FAILED: " << errs[1]<< std::endl;
      errcode += 1;
      } else {
        std::cout << "main: ---------------------------GDivOp OK: " << errs[1] << std::endl;
        errcode += 0;
      }
    }



    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("div_err.txt");
    ios.open("div_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "#  idir   p      num_elems  ncyc  inf_err  L2_err  tw" << std::endl;
    }
    itst.close();

    ios << idir << "  " << pvec[0] ;
        for ( auto j=1; j<GDIM; j++ ) ios << " " <<  pvec[j]; 
    ios << "  " << grid_->ngelems() 
        << "  " << ncyc
        << "  " << errs[0] << "  " << errs[1] << "  " << time
        << std::endl;
    ios.close();
 
    pio::finalize();
    GEOFLOW_TRACE_STOP();
    GComm::TermComm();
    GEOFLOW_TRACE_FINALIZE();

    return( errcode );

} // end, main

