//==================================================================================
// Module       : gtest_rk.cpp
// Date         : 2/24/19 (DLR)
// Description  : GeoFLOW test of GExRK solver
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
#include "gexrl_stepper.hpp"

using namespace geoflow::tbox;
using namespace std;

typedef State GTVector<GTVector<GFTYPE>*>;
typedef Time  GFTYPE;

GFTYPE omega=1.0;

void init_ggfx(GGrid &grid, GGFX &ggfx);
void update_dirichlet(const Time &t, State &u, State &ub);
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt);


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GINT    errcode, iopt;
    GINT    nstate=GDIM;  // number 'state' arrays
    GSIZET  nstage=2, maxSteps=100;
    GFTYPE  dt=1.0e-2, t, maxerror;

    /r option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "d:m:n:w:h")) != -1) {
      switch (iopt) {
      case 'd': // set dt
          dt = atof(optarg);
          break;
      case 'm': // set max timesteps
          maxSteps = atoi(optarg);
          break;
      case 'n': // set no. RK stages
          nstage = atoi(optarg);
          break;
      case 'w': // set max timesteps
          omega = atof(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-d dt] [-m maxSteps] [-n #stages] [-w omega]" << std::endl;
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

    // Create GExRK object:
    GExRK gexrk(nstage);

   std::function<void(const Time &t,                    // RHS callback function
                       const State  &uin,
                       const Time &dt,
                       State &dudt)> rhs
                    = [this](const Time &t,
                       const State  &uin,
                       const Time &dt,
                       State &dudt){dudt_impl(t, uin, dt, dudt);};

    gexrk.setRHSfunction(rhs);
    

    // Create state and tmp space:
    GFTYPE tmax = maxSteps*dt;
    GTVector<GTVector<GFTYPE>*> utmp(4);
    GTVector<GTVector<GFTYPE>*> u   (1);
    GTVector<GTVector<GFTYPE>*> uout(1);
    GTVector<GTVector<GFTYPE>*> ua  (1);
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<du  .size(); j++ ) du  [j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<da  .size(); j++ ) da  [j] = new GTVector<GFTYPE>(1);

    // Find solution to 
    // du/dt = sin(omega t); 
    ua[0] = -(cos(omega*tmax) - cos(0));
    u [0] = 0.0;

    for ( GSIZET k=0; k<maxSteps; k++ ) {
      gexrk.step(u, ub, dt, utmp, uout);
      u = uout;
      t += dt;
    }

    maxerror = fabs(ua[0]-u[0])/ua[0];
    cout << "main: error ua=" << ua 
         <<            " u =" << u << endl;
   
    if ( maxerror > 10*std::numeric_limits<GFTYPE>::epsilon() ) {
      std::cout << "main: -------------------------------------derivative FAILED" << std::endl;
      errcode = 1;
    } else {
      std::cout << "main: -------------------------------------derivative OK" << std::endl;
      errcode = 0;
    }

    return( errcode );

} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD: dudt
// DESC  : RHS function
// ARGS  : 
//**********************************************************************************
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt);
{
  dudt[0] = sin(omega * t); 

} // end method dudt
