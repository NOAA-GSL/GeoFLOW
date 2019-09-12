/*
 * gtest_gvector.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: bflynt
 */


#include "gtmatrix.hpp"


int main(int argc, char **argv) {

    GTMatrix<GFLOAT> A(4,4);
    GTMatrix<GFLOAT> B(4,4);
    GTMatrix<GFLOAT> C(4,4);
    GTMatrix<GFLOAT> D(4,4);

    A(0,0) = 1 ; A(0,1) = 2 ; A(0,2) = 3 ; A(0,3) = 4 ;
    A(1,0) = 5 ; A(1,1) = 6 ; A(1,2) = 7 ; A(1,3) = 8 ;
    A(2,0) = 9 ; A(2,1) = 10; A(2,2) = 11; A(2,3) = 12;
    A(3,0) = 13; A(3,1) = 14; A(3,2) = 15; A(3,3) = 16;

    B(0,0) = 1 ; B(0,1) = 2 ; B(0,2) = 3 ; B(0,3) = 4 ;
    B(1,0) = 5 ; B(1,1) = 6 ; B(1,2) = 7 ; B(1,3) = 8 ;
    B(2,0) = 9 ; B(2,1) = 10; B(2,2) = 11; B(2,3) = 12;
    B(3,0) = 13; B(3,1) = 14; B(3,2) = 15; B(3,3) = 16;

    C = A;

    return 0;
}


