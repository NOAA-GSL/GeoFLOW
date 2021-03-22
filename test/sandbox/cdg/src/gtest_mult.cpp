//==================================================================================
// Module       : gtest_mult.cpp
// Date         : 3/19/21 (DLR)
// Description  : GeoFLOW independent test for GGFX 
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

    GBOOL   usebdy;
    GINT    errcode=0 ;
    GFTYPE  eps=100.0*std::numeric_limits<GFTYPE>::epsilon();
    GString sgrid;// name of JSON grid object to use
    GTVector<GINT>
            pvec;
    std::vector<GINT> 
            pstd(GDIM);  
    GTVector<GINT>    dim(2);
    GTVector<GFTYPE> *mass;
    GTPoint<GFTYPE> dP, P0, P1;
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

    /////////////////////////////////////////////////////////////////
    /////////////////////// Initialize system ///////////////////////
    /////////////////////////////////////////////////////////////////
    ptree.load_file("mult_test.jsn");    // main param file structure



    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    assert(sgrid == "grid_box"); // use only Cartesian grids
    assert(pstd.size() >= GDIM); 


    pvec.resize(pstd.size()); pvec = pstd; pvec.range(0,GDIM-1);

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( auto k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }
    

    // Create grid:
    pstd = gridptree.getArray<GINT>("num_elems");
    assert(GDIM==2 && pstd[0] >=3 && pstd[1] >= 3);

    ObsTraitsType binobstraits;
    grid_ = GGridFactory<MyTypes>::build(ptree, gbasis, pIO_, binobstraits, comm);
    mass = grid_->massop().data();


    // Create GGFX operator:
    GGFX<GFTYPE> ggfx;
    std::vector<std::array<GFTYPE,GDIM>> xyz(grid_->ndof());
    const auto ndof = grid_->ndof(); 
    for(std::size_t i = 0; i < ndof; i++){
      for(std::size_t d = 0; d < GDIM; d++){
        xyz[i][d] = grid_->xNodes()[d][i];
      }
    }
  
    // Create GGFX 
    ggfx.init(0.25*grid_->minnodedist(), xyz);

    // Compute multiplicity:
    GTVector<GFTYPE> amult(grid_->ndof());
    GTVector<GFTYPE> dmult(grid_->ndof());
    GTVector<GFTYPE> mult(grid_->ndof());

    mult = 1.0;
    ggfx.doOp(mult,GGFX<GFTYPE>::Sum());

    ggfx.doOp(mult,GGFX<GFTYPE>::Smooth());

    // Assume grid of at least 3x3 elements:


    // Compute 'analytic' multiplicity:
    GSIZET          n=0, ne=0;
    GTVector<GINT> *testty = &grid_->testty();
    GTVector<GINT> *testid = &grid_->testid();
    GElemList      *gelems     = &grid_->elems(); 
    assert(testid->size() == gelems->size());

    amult = 1.0;
    for ( auto e=0; e<gelems->size(); e++ ) {
      dim = (*gelems)[e]->dim();
      for ( auto j=0; j<dim[1]; j++ ) {
        for ( auto i=0; i<dim[0]; i++ ) {
          if      ( (*testty)[ne] == 0  ) { // vert elem
            if      ( (*testid)[ne] == 0 ) {
              if      ( i == 0        && j == 0        ) amult[n] = 1;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 2;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == dim[0]-1                  ) amult[n] = 2;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 1 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 2;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 1;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0                         ) amult[n] = 2;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 2 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 4;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 2;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 1;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0                         ) amult[n] = 2;  
              else if ( j == 0                         ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 3 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 2;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 4;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 1;  
              else if ( i == dim[0]-1                  ) amult[n] = 2;  
              else if ( j == 0                         ) amult[n] = 2;  
            }
          }
          else if ( (*testty)[ne] == 1  ) { // non-vert bdy
            if      ( (*testid)[ne] == 0 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 2;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 2;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 4;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
              else if ( i == 0                         ) amult[n] = 2;  
              else if ( i == dim[0]-1                  ) amult[n] = 2;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 1 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 4;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 2;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0                         ) amult[n] = 2;  
              else if ( j == 0                         ) amult[n] = 2;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 2 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 4;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 4;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == 0                         ) amult[n] = 2;  
              else if ( i == dim[0]-1                  ) amult[n] = 2;  
              else if ( j == 0                         ) amult[n] = 2;  
            }
            else if ( (*testid)[ne] == 3 ) { 
              if      ( i == 0        && j == 0        ) amult[n] = 2;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 4;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 2;  
              else if ( i == dim[0]-1                  ) amult[n] = 2;  
              else if ( j == 0                         ) amult[n] = 2;  
              else if ( j == dim[1]-1                  ) amult[n] = 2;  
            }
          }
          else if ( (*testty)[ne] == 2  ) { // interior elem
              if      ( i == 0        && j == 0        ) amult[n] = 4;
              else if ( i == dim[0]-1 && j == 0        ) amult[n] = 4;  
              else if ( i == dim[0]-1 && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0        && j == dim[1]-1 ) amult[n] = 4;  
              else if ( i == 0                         ) amult[n] = 2; 
              else if ( i == dim[0]-1                  ) amult[n] = 2; 
              else if ( j == 0                         ) amult[n] = 2; 
              else if ( j == dim[1]-1                  ) amult[n] = 2; 
          }
          n++;
        } // end, i-loop
      } // end, j-loop
      ne++;
   
    } // end, e-loop
    


    // Find inf-norm and L2-norm errors for each method::
    GTVector<GFTYPE> errs(1); // for each method, Linf and L2 errs
    StateComp        lnorm(1), gnorm(1);

    /////////////////////////////////////////////////////////////////
    //////////////////////// Compute Errors /////////////////////////
    /////////////////////////////////////////////////////////////////
    dmult = amult - mult;

cout << "main:  mult=" << mult << endl;
cout << "main: amult=" << amult << endl;
cout << "main: dmult=" << dmult << endl;

    lnorm[0] = dmult.amax();

    // Accumulate to find global inf-norm:
    GComm::Allreduce(lnorm.data(), gnorm.data(), 1, T2GCDatatype<GFTYPE>(), GC_OP_MAX, comm);
    errs[0] = dmult.amax();


    if ( myrank == 0 ) {
      if ( errs[0] > eps ) {
        std::cout << "main: ---------------------------FAILED: " << errs[0]<< std::endl;
      errcode += 1;
      } else {
        std::cout << "main: ---------------------------OK: " << errs[0] << std::endl;
        errcode += 0;
      }
    }



    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("mult_err.txt");
    ios.open("mult_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "#  p      num_elems  inf_err " << std::endl;
    }
    itst.close();

    ios << "  " << pvec[0] ;
        for ( auto j=1; j<GDIM; j++ ) ios << " " <<  pvec[j]; 
    ios << "  " << grid_->ngelems() 
        << "  " << errs[0] 
        << std::endl;
    ios.close();
 
    pio::finalize();
    GEOFLOW_TRACE_STOP();
    GComm::TermComm();
    GEOFLOW_TRACE_FINALIZE();

    return( errcode );

} // end, main

