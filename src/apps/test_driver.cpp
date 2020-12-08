//==================================================================================
// Module       : test_driver
// Date         : 12/07/20 (DLR)
// Description  : Simplified driver
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================


#include "gexec.h"
#include "gtypes.h"
#include "gcomm.hpp"
#include "gmass.hpp"
#include "gmtk.hpp"
#include "gmtk.hpp"
#include "gcblas.hpp"
#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/tracer.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <unistd.h>


using namespace geoflow::tbox;
using namespace std;


int main(int argc, char **argv)
{
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();
    pio::initialize(myrank);
    GEOFLOW_TRACE();


    /////////////////////////////////////////////////////////////////
    ////////////////////// Allocate arrays //////////////////////////
    /////////////////////////////////////////////////////////////////

    // Create state and tmp space:
    GINT Ne=1;
    GINT M=10, N=M*Ne, K=M, lda=M, ldb=K, ldc=M;
    GTVector<double> u(Ne*M*M);
    GTVector<double> y(Ne*M*M);
    GTMatrix<double> D(M,M);
    GCBLAS::GBlasHandle h;
    
    GCBLAS::handle_create(h);
    u = 1.0;
    for ( auto i=0; i<D.dim(1); i++) {
       for ( auto j=0; j<D.dim(2); j++ ) {
	  D(i,j) = 2*i+ j;
       }
    }
    cout << "main: D=" << D << endl;
    cout << "main: u=" << u << endl;
    GCBLAS::gemm<double>(h, GCBLAS::CblasRowMajor, GCBLAS::CblasNoTrans, GCBLAS::CblasNoTrans,
                    M, N, K, 1.0, (double*)D.data().data(), lda, (double*)u.data(), ldb, 0.0, (double*)y.data(), ldc);

    cout << "main: y = D u=" << y << endl;

    GCBLAS::handle_destroy(h);

    pio::finalize();
    GComm::TermComm();

    return(0);
} // end, main

