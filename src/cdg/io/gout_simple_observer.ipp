//==================================================================================
// Module       : ggio_simple.cpp
// Date         : 3/18/19 (DLR)
// Descrption   : Observer object fr carrying out simple non-MPI-based 
//                I/O
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : nstage: number stages not necessarily == truncation order
//**********************************************************************************
template<typename T>
GGIOSimple<T>::GGIOSimple(GSIZET nstage)
:
bgrid_printed_        (FALSE)
{
} // end of constructor (1) method



//**********************************************************************************
//**********************************************************************************
// METHOD     : observer
// DESCRIPTION: Prints state to files specified by traits. Format is:
//                  var1.CCCCCC.TTTT.out,
//              where CCCCCC represents a cycle number, and TTTT represents
//              the mpi task doing the writing.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGIOSimple<T>::observe(const Traits &traits, const Time t, const State &u)
{
  
  gio(*grid_, u, traits.istate, traits.icycle, traits.state_names, bprgrid_);
  
} // end of method observe


