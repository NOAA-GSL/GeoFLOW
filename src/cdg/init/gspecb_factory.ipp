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
// METHOD : dospec
// DESC   : Do bdy specification (ids and types)
// ARGS   : ptree : specification property tree; not main prop tree
//          grid  : GGrid object
//          ibdy  : indirection array into state indicating global bdy
//          tbdy  : array of size ibdy.size giving bdy condition type, returned
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GSpecBFactory<EquationType>::dospec(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("specb_block","");

  // ibdy and tbdy should not come in empty. But they
  // generally refer to only individual canonical boundaries 
  // (individual faces for boxes, individual surfaces for 
  // icos spheres, etc.), and not the complete list of global bdys:

  if ( "specb_none" == sinit
    || "none"       == sinit
    || ""           == sinit 
    || ibdy.size()  == 0     ) {
    bret = TRUE;
  }
  else if ( "specb_uniform"    == sinit ) {
    bret = gspecb::impl_uniform      (ptree, grid, ibdy, tbdy);
  }
  else if ( "mybdyspec    "    == sinit ) {
    bret = gspecb::impl_mybdyspec    (ptree, grid, ibdy, tbdy);
  }
  else                                        {
    bret = assert(FALSE && "Boundary specification method unknown");
  }

  return bret;

} // end, init method



} // namespace pdeint
} // namespace geoflow

