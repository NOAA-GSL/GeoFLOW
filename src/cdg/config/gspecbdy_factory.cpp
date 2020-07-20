//==================================================================================
// Module       : gspecbdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include "gspecbdy_factory.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : dospec
// DESC   : Do bdy specification (ids and types)
// ARGS   : sptree: specification property tree; not main prop tree
//          grid  : GGrid object
//          id    : can serve as boundary id (which canonical bdy)
//          ibdy  : indirection array into state indicating global bdy
// RETURNS: none.
//**********************************************************************************
GBOOL GSpecBdyFactory::dospec(const geoflow::tbox::PropertyTree& sptree, GGrid &grid, const GINT id, GTVector<GSIZET> &ibdy)
{
  GBOOL         bret    = FALSE;
  GSIZET        nbdy;
  GString       sclass  = sptree.getValue<GString>("bdy_class");
  GString       sinit   = sptree.getValue<GString>("bdy_config_method","");

  // ibdy should not come in empty. But they
  // generally refer to only individual canonical boundaries 
  // (individual faces for boxes, individual radial surfaces for 
  // icos spheres, etc.), and not the complete list of global bdys:

/*
  // If bdy_class is uniform, don't need config method:
  if ( "uniform"    == sclass) { 
    bret = gspecbdy::impl_uniform    (sptree, grid, id, ibdy);
    return bret;
  }
*/


  // Else, is mixed; call specified config method.
  // These methods allow user to _both_ configure
  // bdy and to update it, if desired:
  if ( "specb_none" == sinit
    || "none"       == sinit
    || ""           == sinit 
    || ibdy.size()  == 0     ) {
    bret = TRUE;
  }
  else if ( "mixed"        == sinit ) {
    
//  bret = gspecbdy::mybdyspec(sptree, grid, id, ibdy);
    bret = FALSE;
  }
  else                                        {
    assert(FALSE && "Boundary specification method unknown");
  }

  return bret;

} // end, dospec method


