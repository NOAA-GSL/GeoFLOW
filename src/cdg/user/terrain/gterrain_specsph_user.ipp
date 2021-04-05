// Module       : gterrain_specsph.hpp
// Date         : 7/11/19 (DLR)
// Description  : Terrain specification function implementations provided
//                by user for spherical grids (2d and 3d).
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gterrain_specsph_user.hpp"


namespace gterrain_specsph {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_gauss_sphrange
// DESC   : Specify terrange as a 'Gaussian mountain range' with
//          info by ttree. Intended for sph grids
// ARGS   : ptree: main prop tree 
//          sblk : data block name
//          grid : grid
//          utmp : tmp arrays
//          xb   : terrain vectors
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_gauss_range(const PropertyTree &ptree, GString sblk, GGrid &grid, State &utmp, State &xb)
{

  assert(FALSE && "Under construction");

  PropertyTree ttree     = ptree.getPropertyTree(sblk);
  PropertyTree sphptree  = ptree.getPropertyTree("grid_sph");

  GSIZET nxy;
  GFTYPE dx, dy, eps;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  GTPoint<GFTYPE>  P0(3);
  
  nxy = (*xnodes)[0].size();

  std::vector<GFTYPE> x0   = ttree.getArray<GFTYPE>("x0");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> y0   = ttree.getArray<GFTYPE>("y0");
  std::vector<GFTYPE> xsig = ttree.getArray<GFTYPE>("xsigma");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> ysig = ttree.getArray<GFTYPE>("ysigma");
  std::vector<GFTYPE> h0   = ttree.getArray<GFTYPE>("h0");

  std::vector<GFTYPE> xyz0 = sphptree.getArray<GFTYPE>("xyz0");
//std::vector<GFTYPE> dxyz = sphptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; 

  eps =  100*std::numeric_limits<GFTYPE>::epsilon();

  // Build terrain height vector:
#if defined(_G_IS2D)
#elif defined(_G_IS3D)
#endif

  return TRUE;

} // end, impl_gauss_range


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_poly_range
// DESC   : Specify terrange as a inverse polynomials of the form:
//               h = h0/[(x-x0)^2/a^2 + (y-y0)^2/b^2]^(p/2)
//           Intended for sph grids.
//          
// ARGS   : ptree: main prop tree 
//          sblk : data block name
//          grid : grid
//          utmp : tmp arrays
//          xb   : terrain vectors
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_poly_range(const PropertyTree &ptree, GString sblk, GGrid &grid, State &utmp, State &xb)
{

  assert(FALSE && "Under construction");

  PropertyTree ttree     = ptree.getPropertyTree(sblk);
  PropertyTree sphptree  = ptree.getPropertyTree("grid_sph");

  GSIZET nxy;
  GFTYPE dx, dy, eps;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  GTPoint<GFTYPE>  P0(3);
  
  nxy = (*xnodes)[0].size();

  std::vector<GFTYPE> pexp = ttree.getArray<GFTYPE>("exponent");
  std::vector<GFTYPE> x0   = ttree.getArray<GFTYPE>("x0");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> y0   = ttree.getArray<GFTYPE>("y0");
  std::vector<GFTYPE> xsig = ttree.getArray<GFTYPE>("xsigma");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> ysig = ttree.getArray<GFTYPE>("ysigma");
  std::vector<GFTYPE> h0   = ttree.getArray<GFTYPE>("h0");

  std::vector<GFTYPE> xyz0 = sphptree.getArray<GFTYPE>("xyz0");
//std::vector<GFTYPE> dxyz = sphptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; 

  eps =  100*std::numeric_limits<GFTYPE>::epsilon();


  // Build terrain height vector:
#if defined(_G_IS2D)
#elif defined(_G_IS3D)
#endif

  return TRUE;

} // end, impl_poly_range




} // end, ginitbdy namespace
