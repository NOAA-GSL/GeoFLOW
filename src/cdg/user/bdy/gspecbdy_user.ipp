//==================================================================================
// Module       : gspecbdy_user.cpp
// Date         : 7/11/19 (DLR)
// Description  : User-specified boundary specification function implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================





//**********************************************************************************
//**********************************************************************************
// METHOD : impl_my_mixed_bdy
// DESC   : Template for a mixed bdy-specification method to return
//          array of volume indices that define boundary
// ARGS   : sptree : specification prop tree
//          grid   : grid
//          id     : canonical bdy id
//          ibdy   : indirection array into state indicating global bdy, returned
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
template<typename Types>
GBOOL gepecbdy<Types>::impl_my_mixed_bdy(const PropertyTree &sptree, Grid &grid, const GINT id, GTVector<GSIZET> &ibdy)
{
  GBdyType btype = geoflow::str2bdytype(sptree.getValue<GString>("base_type", "GBDY_NONE"));


  return FALSE;

} // end, method impl_my_mixed_bdy


