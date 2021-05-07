//==================================================================================
// Module       : ginitforce_direct_user.cpp
// Date         : 7/10/19 (DLR)
// Description  : Direct user force initialization function implementations. These
//                methods are called directly during configuration, and can set 
//                forcing for entire state (v+b+s) etc. The 'component' types
//                in which component groups (v, b, s, etc) are set individually
//                are contained in separate namespaces.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================




//**********************************************************************************
//**********************************************************************************
// METHOD : initf_impl_null
// DESC   : Initialize for zero force function
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          eqn    : eqnation implementation
//          grid   : grid object
//          t      : time
//          u      : current state
//          uf    : force vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
template<typename Types>
GBOOL ginitforce<Types>::impl_null(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &t, State &utmp, State &u, State &uf)
{

  for ( auto j=0; j<uf.size(); j++ ) {
    if ( uf[j] != NULLPTR ) *uf[j] = 0.0;
  }

  return TRUE;

} // end of method impl_null



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Initialize for random force function
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          eqn    : eqnation implementation
//          grid   : grid object
//          t      : time
//          u      : current state
//          uf    : force vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
template<typename Types>
GBOOL ginitforce<Types>::impl_rand(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &t, State &utmp, State &u, State &uf)
{

  assert(FALSE);

  for ( auto j=0; j<uf.size(); j++ ) {
    if ( uf[j] != NULLPTR ) *uf[j] = 0.0;
  }

  return TRUE;

} // end of method impl_rand



