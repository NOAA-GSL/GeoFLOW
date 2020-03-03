//==================================================================================
// Module       : ginitv.cpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitv.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"

#include <random>

namespace ginitv {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_abc_box
// DESC   : Inititialize velocity with Arnold-Beltrami-Childress (ABC)
//          initial conditions for box grids, 2d and 3d. May also be
//          used with 3D ICOS-based grids.
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_abc_box(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GGridBox *box = dynamic_cast<GGridBox*>(&grid);
  assert(box != NULLPTR && "Box grid required");

  GINT         kdn, kup, p, pdef;
  GSIZET       nn ;
  GFTYPE       A, B, C, E0, pi2, x, y, z;
  PropertyTree vtree ;
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();

#if defined(_G_IS3D)
  pdef = 2;
#else
  pdef = 3;
#endif

  vtree = ptree.getPropertyTree(sconfig);
  kdn   = vtree.getValue<GINT>("kdn");
  kup   = vtree.getValue<GINT>("kup");
  p     = vtree.getValue<GINT>("kpower",pdef);
  A     = vtree.getValue<GFTYPE>("A", 0.9);
  B     = vtree.getValue<GFTYPE>("B", 1.0);
#if defined(_G_IS3D)
  C     = vtree.getValue<GFTYPE>("C", 1.1);
#endif
  E0    = vtree.getValue<GFTYPE>("E0", 1.0);
  nn    = (*xnodes)[0].size();

#if defined(_G_IS2D)
  // Stream fcn is 
  //   psi = Sum_i { -A cos(2pi*ki*x) + B sin(2pi*ki*y) / ki^p }
  // Compute vel components s.t. ux = d psi / dy, uy = -d psi / dx

  *u[0] = 0.0;
  *u[1] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; 
    for ( GINT k=kdn; k<=kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  B*pi2*cos(pi2*y) / pow(k,p);
      (*u[1])[j] += -A*pi2*sin(pi2*x) / pow(k,p);
    }
  }
  
#elif defined(_G_IS3D)

  *u[0] = 0.0;
  *u[1] = 0.0;
  *u[2] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    for ( GINT k=kdn; k<kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  ( B*cos(pi2*y) + C*sin(pi2*z) ) / pow(k,p);
      (*u[1])[j] +=  ( A*sin(pi2*x) + C*cos(pi2*z) ) / pow(k,p);
      (*u[2])[j] +=  ( A*cos(pi2*x) + B*sin(pi2*y) ) / pow(k,p);
    }
  }

#endif

  GMTK::normalizeL2(grid, u, utmp, E0);

  return TRUE;

} // end, method impl_abc_box


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_abc_icos
// DESC   : Inititialize velocity with Arnold-Beltrami-Childress (ABC)
//          initial conditions for icos grids 2d, 3d. 
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_abc_icos(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GGridIcos *icos = dynamic_cast<GGridIcos*>(&grid);
  assert(icos != NULLPTR && "Icos grid required");

  GINT         kdn, kup, p, pdef;
  GSIZET       nn ;
  GFTYPE       A, B, C, E0, pi2, x, y, z;
  GFTYPE       alat, along, r;
  PropertyTree vtree ;
  GTVector<GTVector<GFTYPE>*> 
               usph(GDIM);
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();

#if defined(_G_IS3D)
  pdef = 2;
#else
  pdef = 3;
#endif

  vtree = ptree.getPropertyTree(sconfig);
  kdn   = vtree.getValue<GINT>("kdn");
  kup   = vtree.getValue<GINT>("kup");
  p     = vtree.getValue<GINT>("kpower",pdef);
  A     = vtree.getValue<GFTYPE>("A", 0.9);
  B     = vtree.getValue<GFTYPE>("B", 1.0);
#if defined(_G_IS3D)
  C     = vtree.getValue<GFTYPE>("C", 1.1);
#endif
  E0    = vtree.getValue<GFTYPE>("E0", 1.0);
  nn    = (*xnodes)[0].size();

  usph[0] = utmp[0]; // ulat storage
  usph[1] = utmp[1]; // ulong storage

#if defined(_G_IS2D)
  // Stream fcn is 
  //   psi = Sum_i { -A cos(ki*long) + B sin(ki*lat) / ki^p }
  // Compute vel components s.t. ulat = d psi / dlong, ulong = -d psi / dlat

  *usph[0] = 0.0;
  *usph[1] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    r = sqrt(x*x + y*y+z*z);
    alat = asin(z/r); along = atan2(y,x); 
    for ( GINT k=kdn; k<=kup; k++ ) {
#if 1
      (*usph[1])[j] +=  B*k*cos(k*along) / pow(k,p); // lat
      (*usph[0])[j] += -A*k*sin(k*alat) / pow(k,p);  // long
#else
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  ( B*cos(pi2*y) + C*sin(pi2*z) ) / pow(k,p);
      (*u[1])[j] +=  ( A*sin(pi2*x) + C*cos(pi2*z) ) / pow(k,p);
      (*u[2])[j] +=  ( A*cos(pi2*x) + B*sin(pi2*y) ) / pow(k,p);
#endif
    }
  }
  GMTK::vsphere2cart(grid, usph, GVECTYPE_PHYS, u);
  
#elif defined(_G_IS3D)

  *u[0] = 0.0;
  *u[1] = 0.0;
  *u[2] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    for ( GINT k=kdn; k<kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  ( B*cos(pi2*y) + C*sin(pi2*z) ) / pow(k,p);
      (*u[1])[j] +=  ( A*sin(pi2*x) + C*cos(pi2*z) ) / pow(k,p);
      (*u[2])[j] +=  ( A*cos(pi2*x) + B*sin(pi2*y) ) / pow(k,p);
    }
  }

#endif
 
  GMTK::constrain2sphere(grid, u);
  GMTK::normalizeL2(grid, u, utmp, E0);


  return TRUE;

} // end, method impl_abc_icos



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_simpsum1d_box
// DESC   : Inititialize velocity with simple sum of waves in x-dir only,
//          scaled by k^p. For box grids, 2d mimicking 1d
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_simpsum1d_box(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GGridBox *tgrid = dynamic_cast<GGridBox*>(&grid);
  assert(tgrid != NULLPTR && "Box grid required");

  GINT         kdn, kup, pdef;
  GSIZET       nn ;
  GFTYPE       E0, kn, knh, L, p, x, y, z;
  GFTYPE       knx, kny, knz, knxh, knyh, knzh;
  GFTYPE       phase1, phase2, phase3;

  GTPoint<GFTYPE>
               G0(2), G1(2);
  PropertyTree vtree ;
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();
  std::default_random_engine generator;
  std::normal_distribution<GFTYPE> *distribution;

#if defined(_G_IS3D)
  pdef = 2;
#else
  pdef = 3;
#endif

  G0 = tgrid->getP0();
  G1 = tgrid->getP1();

  vtree  = ptree.getPropertyTree(sconfig);
  kdn    = vtree.getValue<GINT>("kdn");
  kup    = vtree.getValue<GINT>("kup");
  p      = vtree.getValue<GFTYPE>("kpower",pdef);
  E0     = vtree.getValue<GFTYPE>("E0", 1.0);
  nn     = (*xnodes)[0].size();

//distribution = new normal_distribution<GFTYPE>(0,sqrt(2.0*E0));
  distribution = new normal_distribution<GFTYPE>(0,2.0*PI);

#if defined(_G_IS2D)
  // Stream fcn is 
  //   psi = Sum_i { -cos(2pi*ki*x) ) / ki^p }
  // Compute vel components s.t. ux = d psi / dy, uy = -d psi / dx

  L = G1.x1 - G0.x1;
  *u[0] = 0.0;
  *u[1] = 0.0;
  for ( GINT k=kdn; k<=kup; k++ ) {
    kn  = 2.0*PI*static_cast<GFTYPE>(k)/L;
    knh = 2.0*PI*static_cast<GFTYPE>(k+0.5)/L;
    phase1 = (*distribution)(generator);
    phase2 = (*distribution)(generator);
    for ( GSIZET j=0; j<nn; j++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; 
//    (*u[0])[j] +=  (cos(kn*x+phase1) + sin(kn*x+phase2)) / pow(kn,p);
      (*u[0])[j] +=  ( sin(kn*x+phase1) + 4.0*sin(knh*x+phase2) ) / pow(kn,p);
    }
  }
#elif defined(_G_IS3D)
  assert(FALSE && "method intended for 2d mimicking 1d only");
#endif

  GMTK::normalizeL2(grid, u, utmp, E0);

  delete distribution;

  return TRUE;
} // end, method impl_simpsum1d_box

  
//**********************************************************************************
//**********************************************************************************
// METHOD : impl_simpsum_box
// DESC   : Inititialize velocity with simple sum of waves in x-dir only,
//          scaled by k^p. For box grids, 2d mimicking 1d
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_simpsum_box(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GGridBox *tgrid = dynamic_cast<GGridBox*>(&grid);
  assert(tgrid != NULLPTR && "Box grid required");

  GINT         kdn, kup, pdef;
  GSIZET       nn ;
  GFTYPE       E0, L[GDIM], p, x, y, z;
  GFTYPE       kn, knx, kny, knz, knxh, knyh, knzh;
  GFTYPE       phase1, phase2, phase3;
  GTPoint<GFTYPE>
               G0(GDIM), G1(GDIM);
  PropertyTree vtree ;
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();
  std::default_random_engine
               generator;
  std::normal_distribution<GFTYPE>

              *distribution;

#if defined(_G_IS3D)
  pdef = 2;
#else
  pdef = 3;
#endif

  G0 = tgrid->getP0();
  G1 = tgrid->getP1();

  vtree  = ptree.getPropertyTree(sconfig);
  kdn    = vtree.getValue<GINT>("kdn");
  kup    = vtree.getValue<GINT>("kup");
  p      = vtree.getValue<GFTYPE>("kpower",pdef);
  E0     = vtree.getValue<GFTYPE>("E0", 1.0);
  nn     = (*xnodes)[0].size();

//distribution = new normal_distribution<GFTYPE>(0,sqrt(2.0*E0));
  distribution = new normal_distribution<GFTYPE>(0,2.0*PI);

#if defined(_G_IS2D)
  // Stream fcn is 
  //   psi = Sum_i { -cos(2pi*ki*x) ) / ki^p }
  // Compute vel components s.t. ux = d psi / dy, uy = -d psi / dx

  L[0] = G1.x1 - G0.x1;
  L[1] = G1.x2 - G0.x2;
  *u[0] = 0.0;
  *u[1] = 0.0;
  for ( GINT ky=kdn; ky<=kup; ky++ ) {
    kny  = 2.0*PI*static_cast<GFTYPE>(ky)/L[1];
    knyh = 2.0*PI*static_cast<GFTYPE>(ky+0.5)/L[1];
    phase2 = (*distribution)(generator);
    for ( GINT kx=kdn; kx<=kup; kx++ ) {
      knx    = 2.0*PI*static_cast<GFTYPE>(kx)/L[0];
      knxh   = 2.0*PI*static_cast<GFTYPE>(kx+0.5)/L[0];
      kn     = sqrt(knx*knx + kny*kny);
      phase1 = (*distribution)(generator);
      for ( GSIZET j=0; j<nn; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j];
        (*u[0])[j] +=  ( sin(knx*x+phase1) + 4.0*sin(knxh*x+phase1) ) / pow(kn,p);
        (*u[1])[j] +=  ( sin(kny*y+phase2) + 4.0*sin(knyh*y+phase2) ) / pow(kn,p);
      }
    }
  }

#elif defined(_G_IS3D)

  L[0] = G1.x1 - G0.x1;
  L[1] = G1.x2 - G0.x2;
  L[2] = G1.x3 - G0.x3;
  *u[0] = 0.0;
  *u[1] = 0.0;
  *u[2] = 0.0;
  for ( GINT m=0; m<3; m++ ) {
    for ( GINT kz=kdn; kz<=kup; kz++ ) {
      knz  = 2.0*PI*static_cast<GFTYPE>(kz)/L[2];
      knzh = 2.0*PI*static_cast<GFTYPE>(kz+0.5)/L[2];
      phase3 = (*distribution)(generator);
      for ( GINT ky=kdn; ky<=kup; ky++ ) {
        kny  = 2.0*PI*static_cast<GFTYPE>(ky)/L[1];
        knyh = 2.0*PI*static_cast<GFTYPE>(ky+0.5)/L[1];
        phase2 = (*distribution)(generator);
        for ( GINT kx=kdn; kx<=kup; kx++ ) {
          knx  = 2.0*PI*static_cast<GFTYPE>(kx)/L[0];
          knxh = 2.0*PI*static_cast<GFTYPE>(kx+0.5)/L[0];
          kn  = sqrt(knx*knx + kny*kny + knz*knz);
          phase1 = (*distribution)(generator);
          for ( GSIZET j=0; j<nn; j++ ) {
            x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
            (*u[m])[j] +=  ( sin(knx*x+phase1) + 4.0*sin(knyh*y+phase2) + 2.0* sin(knzh*z+phase3) ) / pow(kn,p);
          }
        }
      }
    }
  }

#endif


  GMTK::normalizeL2(grid, u, utmp, E0);

  delete distribution;

  return TRUE;

} // end, method impl_simpsum_box


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_simpsum_icos
// DESC   : Inititialize velocity with simple sum of waves in lat-long
//          scaled by k^p. For icos grids, 2d and 3d.
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_simpsum_icos(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GGridIcos *tgrid = dynamic_cast<GGridIcos*>(&grid);
  assert(tgrid != NULLPTR && "Icos grid required");

  GINT         kdn, kup, pdef;
  GSIZET       nn ;
  GFTYPE       E0, kn, p, r, x, y, z;
  GFTYPE       lat, lon;
  GFTYPE       phase1, phase2, phase3;
  PropertyTree vtree ;
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();
  GTVector<GTVector<GFTYPE>*> 
               usph(GDIM);
  std::default_random_engine generator;
  std::normal_distribution<GFTYPE> *distribution;

#if defined(_G_IS3D)
  pdef = 2;
  assert(FALSE && "Method not yet ready for 3d");
#else
  pdef = 3;
#endif


  vtree  = ptree.getPropertyTree(sconfig);
  kdn    = vtree.getValue<GINT>("kdn");
  kup    = vtree.getValue<GINT>("kup");
  p      = vtree.getValue<GFTYPE>("kpower",pdef);
  E0     = vtree.getValue<GFTYPE>("E0", 1.0);
  nn     = (*xnodes)[0].size();

  distribution = new normal_distribution<GFTYPE>(0,2.0*PI);

  usph[0] = utmp[0]; // ulat storage
  usph[1] = utmp[1]; // ulong storage

  *usph[0] = 0.0;
  *usph[1] = 0.0;
//*u[2] = 0.0;
  for ( GINT k=kdn; k<=kup; k++ ) {
    kn = static_cast<GFTYPE>(k);
    phase1 = (*distribution)(generator);
    phase2 = (*distribution)(generator);
    phase3 = (*distribution)(generator);
    for ( GSIZET j=0; j<nn; j++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j]; 
      r = sqrt(x*x + y*y + z*z);
      lat = asin(z/r); lon = atan2(y,x);
#if 0
      (*usph[0])[j] +=  (sin(kn*lat+phase1) + 4.0*sin((kn+0.5)*lat+phase1)) / pow(kn,p);
      (*usph[1])[j] +=  (sin(kn*lon+phase2) + 4.0*sin((kn+0.5)*lon+phase2)) / pow(kn,p);
#else
      (*usph[0])[j] +=  cos(kn*lon+phase3)*(sin(kn*lat+phase1) + 4.0*sin((kn+0.5)*lat+phase1)) / pow(kn,p);
      (*usph[1])[j] +=  sin(kn*lon+phase3)*(kn*cos(kn*lat+phase1) + 4.0*(kn+0.5)*cos((kn+0.5)*lat+phase1)) / pow(kn,p+1);
#endif
    } // end, j-loop
  } // end, k loop
  
  GMTK::vsphere2cart(grid, usph, GVECTYPE_PHYS, u);

  GMTK::constrain2sphere(grid, u);
  GMTK::normalizeL2(grid, u, utmp, E0);

  delete distribution;

  return TRUE;

} // end, method impl_simpsum_icos



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_rand(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  return FALSE;

} // end, method impl_rand




} // end, ginitv  namespace
