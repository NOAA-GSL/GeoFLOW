//==================================================================================
// Module       : ginitstate_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : spec
// DESC   : Do specification of terrain
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          utmp   : tmp arrays
//          xb     : terrain coordinates (of size of a State vector)
//          bterr  : flag telling if terrain was loaded
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GSpecTerrainFactory<Types>::spec(const PropertyTree& ptree, Grid &grid, State &utmp, State &xb, GBOOL &bterr)
{
  GBOOL         bret = FALSE;
  GString       stype;  
  GGridBox     *box = dynamic_cast<GGridBox*>(&grid);
  GGridIcos    *icos= dynamic_cast<GGridIcos*>(&grid);


  // Get type of initialization: by-name or by-block:
  stype = ptree.getValue<GString>("terrain_type","");
  if ( "none"   == stype 
    || ""       == stype ) {
    bterr = FALSE;         // terrain not loaded into xb
    return TRUE;
  }

  // Terrain makes no sense if we don't have deformed elements,
  // so check that here:
  assert(grid.gtype() != GE_REGULAR && "Invalid element type");


  if      ( box ) {
    bret = spec_box   (ptree, grid, utmp, stype, xb);
  }
  else if ( icos ) {
    bret = spec_sphere(ptree, grid, utmp, stype, xb);
  }
  else {
    assert(FALSE && "Invalid specification class or grid type");
  }

  if ( bret ) bterr = TRUE; // terrain loaded

  return bret;

} // end, init method


//**********************************************************************************
//**********************************************************************************
// METHOD : spec_box
// DESC   : Do terrain specification for box grid types
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          utmp   : tmp arrays
//          stype  : terrain type block name
//          xb     : terrain coordinates
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GSpecTerrainFactory<Types>::spec_box(const PropertyTree& ptree, Grid &grid, State &utmp, GString stype, State &xb)
{
  GBOOL         bret    = FALSE;
  GString       sname;
  PropertyTree  blktree = ptree.getPropertyTree(stype);

  sname = blktree.getValue<GString>("name");
  if      ( "boxgauss_range"   == sname ) {
    bret = gterrainSpecBox<Types>::impl_gauss_range   (ptree, stype, grid, utmp, xb);
  }
  else if ( "boxpoly_range"    == sname ) {
    bret = gterrainSpecBox<Types>::impl_poly_range    (ptree, stype, grid, utmp, xb);
  }
  else if ( "boxschar_range"   == sname ) {
    bret = gterrainSpecBox<Types>::impl_schar_range   (ptree, stype, grid, utmp, xb);
  }
  else                                        {
    assert(FALSE && "Terrain specification method unknown");

  }

  return bret;

} // end, spec_box method


//**********************************************************************************
//**********************************************************************************
// METHOD : spec_sphere
// DESC   : Do terrain specification for sphere grid types
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          utmp   : tmp arrays
//          stype  : terrain type block name
//          xb     : terrain coordinates
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GSpecTerrainFactory<Types>::spec_sphere(const PropertyTree& ptree, Grid &grid, State &utmp, GString stype, State &xb)
{
  GBOOL         bret    = FALSE;
  GString       sname;
  PropertyTree  blktree = ptree.getPropertyTree(stype);

  sname = blktree.getValue<GString>("name");

  if      ( "sphgauss_range"   == sname ) {
    bret = gterrainSpecSph<Types>::impl_gauss_range   (ptree, stype, grid, utmp, xb);
  }
  else if ( "sphpoly_range"    == sname ) {
    bret = gterrainSpecSph<Types>::impl_poly_range    (ptree, stype, grid, utmp, xb);
  }
  else                                        {
    assert(FALSE && "Terrain specification method unknown");

  }

  return bret;

} // end, spec_sphere method


