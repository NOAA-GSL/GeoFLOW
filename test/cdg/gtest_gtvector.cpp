/*
 * gtest_gvector.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: bflynt
 */


#include "gtvector.hpp"



int main(int argc, char **argv) {


    GTVector<GFLOAT> d(4);
    GTVector<GFLOAT> e(4);
    GTVector<GFLOAT> u(4);
    GTVector<GFLOAT> v(4);


    d[0] = 3;
    d[1] = 5;
    d[2] = 7;
    d[3] = 9;

    e[0] = 3;
    e[1] = 5;
    e[2] = 7;
    e[3] = 9;

    u[0] = 3;
    u[1] = 5;
    u[2] = 7;
    u[3] = 9;

    v[0] = 3;
    v[1] = 5;
    v[2] = 7;
    v[3] = 9;

    v = d + e;

    return 0;
}

