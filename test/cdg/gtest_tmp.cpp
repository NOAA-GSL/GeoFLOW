//==================================================================================
// Module       : gtest_blas.cpp
// Date         : 6/30/18 (DLR)
// Description  : GeoFLOW BLAS test
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include "gtvector.hpp"
#include "gllbasis.hpp"
#include "gelem_base.hpp"

int main(int argc, char **argv)
{
    GString serr ="main: ";
    GBOOL   bret;
    GINT    errcode=0;
    GINT    np     = 2;
    GINT    nelems = 10;

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) gbasis[k] = new GLLBasis<GCTYPE,GFTYPE>(np);

    // Create element list:
    GTVector<GElem_base*> gelems;
    for ( GSIZET j=0; j<nelems; j++ ) {
        gelems.push_back(new GElem_base(GE_2DEMBEDDED,gbasis));
    }

term:
    if ( errcode != 0 ) {
      std::cout << serr << " Error: code=" << errcode <<  std::endl;
    }
    else {
      std::cout << serr << "     Success!" << std::endl;
    }

    return(errcode);
} // end, main
