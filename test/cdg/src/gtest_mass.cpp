//==================================================================================
// Module       : gtest_mass.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of GMass classes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"

using namespace geoflow::tbox;
using namespace std;

void init_ggfx(GGrid &grid, GGFX<GFTYPE> &ggfx);

GGrid       *grid_;
GC_COMM      comm_=GC_COMM_WORLD ;      // communicator


int main(int argc, char **argv)
{

    GString  serr ="main: ";
    GINT     errcode, gerrcode;
    GFTYPE   radiusi;

    errcode = 0;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);


    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();

    // Get minimal property tree:
    PropertyTree gtree, ptree; 
    GString sgrid;
    std::vector<GINT> pstd(GDIM);

    ptree.load_file("input.jsn");
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

    // Generate grid:
    GPTLstart("gen_grid");
    grid_ = GGridFactory::build(ptree, gbasis, comm_);
    GPTLstop("gen_grid");

    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(*grid_, ggfx);

    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();

    GTVector<GFTYPE> f(grid_->ndof());
    GTVector<GFTYPE> g(grid_->ndof());
    GTVector<GTVector<GFTYPE>*> utmp(1);
    GMass massop(*grid_);

    f = 1.0;
    ggfx.doOp(f, GGFX_OP_SMOOTH);

    GPTLstart("massop_prod");
    massop.opVec_prod(f,utmp,g);
    GPTLstop("massop_prod");
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
#if 0
    // Multiply f by inverse multiplicity:
    g.pointProd(imult);
#endif

    GFTYPE integral=g.sum();
    GFTYPE gintegral;
    GComm::Allreduce(&integral, &gintegral, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);

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

    ios <<   gsz[0] << "  ";
    for ( GSIZET j=0; j<GDIM; j++ ) ios << pstd[j] << "  ";
    ios << gsz[1]     << "  "
        << gintegral  << "  " 
        << aintegral  << "  "
        << fabs(gintegral-aintegral) << std::endl;
    ios.close();
    

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);

 
    GPTLpr_file("timing.txt");
    GPTLfinalize();


term:
    if ( gerrcode != 0 ) {
      GPP(comm_,serr << " Error: code=" << errcode);
    }
    else {
      GPP(comm_,serr << "     Success!");
    }

    GComm::TermComm();
    return(gerrcode);
} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : grid    : GGrid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(GGrid &grid, GGFX<GFTYPE> &ggfx)
{
  GFTYPE                         delta;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;

  // First, periodize coords if required to, 
  // before labeling nodes:
  if ( typeid(grid) == typeid(GGridBox) ) {
    static_cast<GGridBox*>(&grid)->periodize();
  }

  delta  = grid.minnodedist();
  dX     = 0.1*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);

  gmorton.key(glob_indices, *xnodes);
  GComm::Synch(comm_);

#if 0
  GPP(comm_,"init_ggfx: glob_indices=" << glob_indices);
#endif

  // Initialize gather/scatter operator:
  GBOOL bret;
  bret = ggfx.init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

  // Unperiodize nodes now that connectivity map is
  // generated, so that coordinates mean what they should:
  if ( typeid(grid) == typeid(GGridBox) ) { // periodize coords
    static_cast<GGridBox*>(&grid)->unperiodize();
  }

} // end method init_ggfx

