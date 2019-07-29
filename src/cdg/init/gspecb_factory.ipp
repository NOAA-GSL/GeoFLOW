//==================================================================================
// Module       : gspecb_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do bdy initialization
// ARGS   : ptree  : main property tree
//          grid   : GGrid object
//          ibdy   : indirection array into state indicating global bdy
//          tbdy   : array of size ibdy.size giving bdy condition type, returned
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GSpecBFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, const GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("specb_block","");

  if ( "specb_none" == sinit
    || "none"       == sinit
    || ""           == init 
    || ibdy.size()  == 0   ) {
    return;
  }
  else if ( "specb_box_uniform"    == sinit ) {
    bret = gspecb::impl_box_uniform      (ptree, grid, ibdy, tbdy);
  }
  else if ( "specb_icos_uniform"   == sinit ) {
    bret = gspecb::impl_icos_uniform     (ptree, grid, ibdy, tbdy);
  }
  else                                        {
    bret = assert(FALSE && "Boundary specification method unknown");
  }

  return bret;

} // end, init method



} // namespace pdeint
} // namespace geoflow

