//==================================================================================
// Module       : gtest_derivs.cpp
// Date         : 2/24/19 (DLR)
// Description  : GeoFLOW test of derivatives
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"

using namespace geoflow::tbox;
using namespace std;


GGrid *grid_ = NULLPTR;
void init_ggfx(GGrid &grid, GGFX &ggfx);

int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   np=1;    // elem 'order'
    GINT   nstate=GDIM;  // number 'state' arrays
    GSIZET maxSteps;
    GFTYPE radiusi=1, radiuso=2;
    std::vector<GINT> ne(3); // # elements in each direction in 3d
    GString sgrid;// name of JSON grid object to use
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;     // main prop tree
    PropertyTree gridptree; // grid prop tree

    ptree.load_file("input.jsn");       // main param file structure
    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    np          = ptree.getValue<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);

    ne          = gridptree.getArray<GINT>("num_elems");  // may be modified by command line

#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "i:j:k:l:m:p:h")) != -1) {
      switch (iopt) {
      case 'i': // handled by InputManager
          break;
      case 'j': // get # elements in r/x
          ne[0] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'k': // get # elements in lat/y
          ne[1] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'm': // get # elements in long/z
          ne[2] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          gridptree.setValue<GINT>("ilevel",ilevel);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          ptree.setValue<GINT>("exp_order",np);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-j #Elems in x/r] [-k #Elems in y/lat]  [-m #Elems in z/long] [-l refine level] -p expansion order] " << std::endl;
          std::cout << "Note: Icos grid only requires refine level to specify number of elements. " << std::endl;
          exit(1); 
          break;
      case ':': // missing option argument
          std::cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
          exit(1); 
          break;
      case '?':
      default: // invalid option
          std::cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
          exit(1);
          break;
      }
    }

#endif

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(np);
    }
    
    GPTLstart("gen_grid");
    // Create grid:
    grid_ = GGridFactory::build(gridptree, gbasis, comm);
    GPTLstop("gen_grid");


    GPTLstart("do_gather_op");
    // Initialize gather/scatter operator:
    GGFX ggfx;
    init_ggfx(*grid_, ggfx);
    GPTLstop("do_gather_op");


    // Create state and tmp space:
    GTVector<GTVector<GFTYPE>*> utmp(4);
    GTVector<GTVector<GFTYPE>*> u   (1);
    GTVector<GTVector<GFTYPE>*> du  (GDIM);
    GTVector<GTVector<GFTYPE>*> da  (GDIM);
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<du  .size(); j++ ) du  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<da  .size(); j++ ) da  [j] = new GTVector<GFTYPE>(grid_->size());

    // Initialize u:
    GFTYPE x, y, z;
    GTVector<GFTYPE> etmp1;
    GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();   
    GSIZET nxy = nxy = (*xnodes)[0].size();
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      if ( xnodes->size() > 2 ) z = (*xnodes)[2][j];
      (*u [0])[j] = x*y*y;
      (*da[0])[j] = y*y;
      (*da[1])[j] = 2*x*y;
      if ( GDIM > 2 ) (*da[2])[j] = 0.0;
    }
    GMTK::compute_grefderivs(*grid_, *u[0], etmp1, FALSE, du);
    for ( GSIZET j=0; j<GDIM; j++ ) {
      grid_->deriv(*u[0], j+1, *utmp[0], *du[j]);
    }


#if 1
    // Compute analytic solution, do comparisons:
    GTVector<GFTYPE> lerrnorm(3), gerrnorm(3), maxerror(3);
    maxerror = 0.0;
    for ( GSIZET j=0; j<du.size(); j++ ) { //local errors
cout << "main: da[" << j << "]=" << *da[j] << endl;
cout << "main: du[" << j << "]=" << *du[j] << endl;
      *utmp[0] = *du[j] - *da[j];
       lerrnorm[0]  = utmp[0]->L1norm (); // inf-norm
       lerrnorm[1]  = utmp[0]->Eucnorm(); // Euclidean-norm
       lerrnorm[1] *= lerrnorm[1]; // square it, so it can be added 
       gerrnorm[2]  = grid_->integrate(*utmp[0],*utmp[1]) /
                      grid_->integrate(*da  [j],*utmp[1]) ; // Global L2 error norm 
      // Accumulate to find global errors for this field:
      GComm::Allreduce(lerrnorm.data()  , gerrnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
      GComm::Allreduce(lerrnorm.data()+1, gerrnorm.data()+1, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm);
      gerrnorm[1] = sqrt(gerrnorm[1]);
      // now find max errors of each type for each field:
      for ( GSIZET i=0; i<3; i++ ) maxerror[i] = MAX(maxerror[i],gerrnorm[j]);
    }

    cout << "main: maxerror = " << maxerror << endl;
   
    GSIZET lnelems=grid_->nelems();
    GSIZET gnelems;
    GComm::Allreduce(&lnelems, &gnelems, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm);

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("deriv_err.txt");
    ios.open("deriv_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p  num_elems   L1_err   Eucl_err   L2_err" << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  "
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();

#endif
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();

    GComm::TermComm();
    if ( grid_ != NULLPTR ) delete grid_;

    return( maxerror.max()  > 10*std::numeric_limits<GFTYPE>::epsilon());

} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather scatter operator
// ARGS  : grid    : GGrid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(GGrid &grid, GGFX &ggfx)
{
  GFTYPE                         delta;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;

  delta  = grid.minnodedist().min();
  dX     = 0.5*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);
  gmorton.key(glob_indices, *xnodes);

  // Initialize gather/scatter operator:
  GBOOL bret = ggfx.init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

} // end method init_ggfx

