//==================================================================================
// Module       : ginitfv.cpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity force inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitfv.hpp"

namespace ginitfv {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : ptree  : initial condition property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_rand(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{

  return FALSE;

} // end, method impl_rand


} // end, ginitfv  namespace
