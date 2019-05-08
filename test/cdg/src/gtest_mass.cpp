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
#include "gtools.h"

using namespace geoflow::tbox;
using namespace std;


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
    init_ggfx(ptree, *grid_, ggfx);

    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();

    GTVector<GFTYPE> f(grid_->ndof());
    GTVector<GFTYPE> g(grid_->ndof());
    GTVector<GFTYPE> *imult = &ggfx.get_imult();
//  GTVector<GFTYPE> *jac = &(grid_->Jac());
    GTVector<GTVector<GFTYPE>*> utmp(1);
    GMass massop(*grid_);

    f = 1.0;
#if 1
//  ggfx.doOp(f, GGFX_OP_SMOOTH);
    // Multiply f by inverse multiplicity:
    f.pointProd(*imult);
cout << serr << "imult=" << *imult << endl;
#endif

    GFTYPE integral;
    GFTYPE gintegral;

#if 0
    GPTLstart("massop_prod");
    massop.opVec_prod(f,utmp,g);
    GPTLstop("massop_prod");
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
    integral = g.sum();
    GComm::Allreduce(&integral, &gintegral, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
#else
    gintegral = grid_->integrate(f, g);
#endif

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


