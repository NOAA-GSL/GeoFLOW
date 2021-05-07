//==================================================================================
// Module       : ginitfs.hpp
// Date         : 7/16/19 (DLR)
// Description  : Scalar forcing inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : ptree  : initial condition property tree
//          sconfig: ptree block name containing variable config
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          uf     : state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
template<typename Types>
GBOOL ginitfs<Types>::impl_rand(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{

  return FALSE;

} // end, method impl_rand


