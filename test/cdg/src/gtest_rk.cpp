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
#include "gexrk_stepper.hpp"

typedef GTVector<GTVector<GFTYPE>*> State;
typedef GFTYPE Time;

GFTYPE omega=1.0;

void update_dirichlet(const Time &t, State &u, State &ub);
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt);


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GINT    errcode, iopt;
    GINT    nstate=GDIM;  // number 'state' arrays
    GSIZET  nstage=2, maxSteps=1;
    GFTYPE  dt=1.0e-2, t, t0=0.0, tmax=3.14159, maxerror;

    // : option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "d:t:n:w:h")) != -1) {
      switch (iopt) {
      case 'd': // set dt
          dt = atof(optarg);
          break;
      case 't': // set max time
          tmax = atof(optarg);
          break;
      case 'n': // set no. RK stages
          nstage = atoi(optarg);
          break;
      case 'w': // set max timesteps
          omega = atof(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-d dt] [-t max time] [-n #stages] [-w omega]" << std::endl;
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
    GExRKStepper<GFTYPE> gexrk(nstage);

   std::function<void(const Time &t,                    // RHS callback function
                       const State  &uin,
                       const Time &dt,
                       State &dudt)> rhs = dudt;

    gexrk.setRHSfunction(rhs);
    

    // Create state and tmp space:
    maxSteps = static_cast<GSIZET>((tmax - t0)/dt);
    GTVector<GTVector<GFTYPE>*> u   (1);
    GTVector<GTVector<GFTYPE>*> uout(1);
    GTVector<GTVector<GFTYPE>*> ua  (1);
    GTVector<GTVector<GFTYPE>*> ub  (1);
    GTVector<GTVector<GFTYPE>*> utmp(u.size()*(nstage+2)+1);
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<uout.size(); j++ ) uout[j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<ua  .size(); j++ ) ua  [j] = new GTVector<GFTYPE>(1);
    for ( GSIZET j=0; j<ub  .size(); j++ ) ub  [j] = new GTVector<GFTYPE>(1);

    // Find solution to 
    // du/dt = sin(omega t) at tmax: 
    *ua[0] = -(cos(omega*tmax) - cos(omega*t0))/omega;
    *u [0] = 0.0; // initialize numerical soln

    t = t0;
    for ( GSIZET k=0; k<maxSteps; k++ ) {
      gexrk.step(t, u, ub, dt, utmp);
//    *u[0] = *uout[0];
      t += dt;
    }

    maxerror = fabs((*ua[0])[0]-(*u[0])[0])/(*ua[0])[0];
    cout << "main: ua=" << *ua[0] << ";  u =" << *u[0] << endl;
   
    cout << "main: relative error=" << maxerror << endl;
    cout << "main: maxSteps      =" << maxSteps << endl;

    if ( maxerror > 10*std::numeric_limits<GFTYPE>::epsilon() ) {
      std::cout << "main: -------------------------------------RK stepper FAILED" << std::endl;
      errcode = 1;
    } else {
      std::cout << "main: -------------------------------------RK stepper OK" << std::endl;
      errcode = 0;
    }

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("rk_err.txt");
    ios.open("rk_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# nstages     dt     tmax     err  " << std::endl;
    }
    itst.close();

    ios << nstage << "  "  <<  dt  << "  " << tmax << "  "  << maxerror  << endl;
    ios.close();

    return(errcode);

} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD: dudt
// DESC  : RHS function
// ARGS  : 
//**********************************************************************************
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt)
{
  (*dudt[0])[0] = sin(omega * t); 

} // end method dudt
