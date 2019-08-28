//==================================================================================
// Module       : ginitstate_user.cpp
// Date         : 7/10/19 (DLR)
// Description  : Direct user state initialization function implementations. These
//                methods are called directly during configuration, and can set 
//                forcing for entire state (v+b+s) etc. The 'component' types
//                in which component groups (v, b, s, etc) are set individually
//                are contained in separate namespaces.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitstate_direct_user.hpp"


namespace ginitstate {



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on box grids with
//          Dirichlet or periodic boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxnwaveburgers(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxnwaveburgers: ";
  GBOOL            bplanar=TRUE; // planar or circularized
  GBOOL            brot   =FALSE;
  GSIZET           i, j, idir, nxy;
  GFTYPE           A, K2, nu, Re, r2, t, t0, tdenom;
  GFTYPE           efact, sum, tfact, xfact;
  GTVector<GFTYPE> K(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);
  std::vector<GFTYPE> kprop;
  GString          snut;

  PropertyTree nwaveptree = ptree   .getPropertyTree("initstate_nwave");
  PropertyTree boxptree   = ptree   .getPropertyTree("grid_box");
  PropertyTree nuptree    = ptree.getPropertyTree("dissipation_traits");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
//    && "INFLOWT boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling. But we won't parameterize in terms of Re, 
  // but rather, nu.
  // Set some parameters:
  r0.x1  = nwaveptree.getValue<GFTYPE>("x0"); 
  r0.x2  = nwaveptree.getValue<GFTYPE>("y0"); 
  r0.x3  = nwaveptree.getValue<GFTYPE>("z0"); 
  A      = nwaveptree.getValue<GFTYPE>("ULparm",1.0);
//Re     = nwaveptree.getValue<GFTYPE>("Re",6.0);
  t0     = nwaveptree.getValue<GFTYPE>("t0",0.04);
  bplanar= nwaveptree.getValue<GBOOL>("planar",TRUE);
  kprop  = nwaveptree.getArray<GFTYPE>("prop_dir");
  nu     = nuptree   .getValue<GFTYPE>("nu",0.0833);
  snut   = nuptree   .getValue<GString>("nu_type","constant");
  K      = kprop;
  K     *= 1.0/K.Eucnorm();

  K2     = 0.0 ; for ( GSIZET i=0; i<GDIM; i++ ) K2 += K[i]*K[i];

  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // If prop direction has more than one component != 0. Then
  // front is rotated (but still planar):
  for ( i=0, brot=TRUE; i<GDIM; i++ ) brot = brot && K[i] != 0.0 ;
  for ( i=0, idir=0; i<GDIM; i++ ) if ( K[i] > 0 ) {idir=i; break;}

  if ( t == 0.0 ) t = K2 * t0;
  Re = A/(2.0*nu); // set Re from nu


  for ( j=0; j<nxy; j++ ) {
    for ( i=0; i<GDIM; i++ ) {
      xx[i] = (*xnodes)[i][j] - r0[i];
      (*u[i])[j] = 0.0;
    }
    if ( bplanar ) { // compute k.r for planar wave
      for ( i=0, sum=0.0; i<GDIM; i++ ) { 
        sum += K[i]*xx[i];
        xx[i] = 0.0;
      }
      xx[0] = sum;
    }
    for ( i=0, r2=0.0; i<GDIM; i++ ) r2 += xx[i]*xx[i];  

    // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t)) ]^-1
    tdenom = 1.0/(4.0*nu*time);
    tfact  = bplanar ? sqrt(t/t0): time/t0;
    efact  = tfact * exp(r2*tdenom) / ( exp(Re) - 1.0 );
    xfact  = 1.0 /( t * (  1.0 + efact ) );
    for ( i=0; i<GDIM; i++ ) (*u[i])[j] = xx[i]*xfact;
    // dU1max = 1.0 / ( time * (sqrt(time/A) + 1.0) );
    // aArea  = 4.0*nu*log( 1.0 + sqrt(A/time) );
  }

  return TRUE;

} // end, impl_boxnwaveburgers



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxdirgauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          Dirichlet boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          time   : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxdirgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxdirgauss: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp;
  GFTYPE           nxy, nu, sig0, u0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  State            c(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);
  GString          snut;

  PropertyTree heatptree = ptree.getPropertyTree("initstate_boxdirgauss");
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[GDIM+j];

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_INFLOWT") >= 2*GDIM
      && "GBDY_INFLOW boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  u0    = heatptree.getValue<GFTYPE>("u0");

  nu     = nuptree   .getValue<GFTYPE>("nu");
  snut   = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<GDIM; j++ ) *c[j] = 0.0;

  if ( bpureadv ) {
     for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];
  }


  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*time*nu); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0/sig[k],GDIM);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) changes, change this term accordingly:
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c[i])[j]*time;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*u[0])[j] = ufact[0]*exp(argxp);
  }

  return TRUE;
} // end, impl_boxdirgauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxpergauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          periodic boundaries
// ARGS   : stree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxpergauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxpergauss: ";
  GBOOL            bContin;
  GSIZET           i, j, k, n;
  GFTYPE           iargp, iargm ;
  GFTYPE           isum , irat , prod;
  GFTYPE           sumn , eps;
  GFTYPE           nxy, nu, pint, sig0, u0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);
  State            c(GDIM);
  GString          snut;

  PropertyTree heatptree = ptree.getPropertyTree("initstate_boxpergauss");
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = 1.0e-4*std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[GDIM+j];

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_PERIODIC") >= 2*GDIM
      && "Periodic boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  u0    = heatptree.getValue<GFTYPE>("u0");

  nu    = nuptree   .getValue<GFTYPE>("nu");
  snut  = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");


  // Set adv velocity components. Note:
  // First state elem is the scalar solution, and
  // the remainder are the velocity components:

  for ( j=0; j<GDIM; j++ ) {
    sig[j] = sqrt(sig0*sig0 + 4.0*time*nu);
    si [j] = 1.0/(sig[j]*sig[j]);
   *c  [j] = 0.0;
  }

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  if ( bpureadv ) for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];

  for ( n=0; n<nxy; n++ ) {

    prod = 1.0;
    for ( k=0; k<GDIM; k++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c(t) changes, change this term accordingly:
      f [k]  = modf((*c[k])[j]*time/gL[k],&pint);
//    f [k]  = (*c[k])[n]*t/gL[k];
      xx[k]  = (*xnodes)[k][n] - r0[k] - f[k]*gL[k];

      isum    = 0.0;
      i       = 0;
      irat    = 1.0;
      while ( irat > eps ) { // inner sum
        iargp   = pow((xx[k]+i*gL[k]),2)*si[k];
        iargm   = pow((xx[k]-i*gL[k]),2)*si[k];
        sumn    = i==0 ? exp(-iargp) : exp(-iargp) + exp(-iargm);
        isum   += sumn;
        irat    = sumn / isum ;
        i++;
      }
      prod *= isum;
    }
    (*u[0])[n] = u0*pow(sig0,GDIM)/pow(sig[0],GDIM)*prod;

  } // end, loop over grid points

  return TRUE;

} // end, impl_boxpergauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosgauss
// DESC   : Initialize state for Burgers with Gauss lump on ICOS grid
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GString          serr = "impl_icosgauss: ";
  GBOOL            bContin;
  GINT             j, k, n, nlumps;
  GSIZET           nxy;
  GFTYPE           alpha, argxp;
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r, s;
  GFTYPE           nu, rad, u0 ;
  GFTYPE           vtheta, vphi;
  GFTYPE           tiny = std::numeric_limits<GFTYPE>::epsilon();
  GTPoint<GFTYPE>           rt(3);
  GTVector<GFTYPE>          xx(3);
  GTVector<GFTYPE>          si(4), sig(4), ufact(4);
  GTVector<GFTYPE>          latp(4), lonp(4);
  std::vector<GFTYPE>       c0(4), sig0(4);
  std::vector<GFTYPE>       lat0(4), lon0(4); // up to 4 lumps
  std::vector<GFTYPE>       Omega;
  State            c(GDIM+1);
  GString          snut;

  PropertyTree lumpptree = ptree.getPropertyTree("initstate_icosgaussburgers");
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");


  for ( GSIZET j=0; j<GDIM+1; j++ ) c[j] = u[GDIM+j];

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  lat0  = lumpptree.getArray<GFTYPE>("latitude0"); // lat for each lump
  lon0  = lumpptree.getArray<GFTYPE>("longitude0"); // lon for each lump
  sig0  = lumpptree.getArray<GFTYPE>("sigma"); // sig for each lump
  c0    = lumpptree.getArray<GFTYPE>("c0");  // initial concentrations for each lump
  rad   = icosptree.getValue<GFTYPE>("radius");
  u0    = lumpptree.getValue<GFTYPE>("u0");
  alpha = lumpptree.getValue<GFTYPE>("alpha",0.0);
  nlumps= lumpptree.getValue<GINT>("nlumps",1);

  alpha *= PI/180.0;

  nu    = nuptree   .getValue<GFTYPE>("nu");
  snut  = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // Convert initial locations from degrees to radians,
  // & compute initial positions of lumps in Cart coords:
  for ( GSIZET k=0; k<nlumps; k++ ) {
    lat0[k] *= PI/180.0;
    lon0[k] *= PI/180.0;
  }

  // Set velocity here. Taken to be solid body rotation,
  // u = Omega X r, where rotation rate vector, Omega
  // is computed by rotating u0 x (0, 0, 1) about x axis by
  // an amount alpha. These point to components of state u_:
  if ( bpureadv ) {
    for ( k=0; k<nxy; k++ ) {
      x   = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r   = sqrt(x*x + y*y + z*z);
      // Compute lat & long:
      lat = asin(z/r);
      lon = atan2(y,x);
      // u_lat = u_theta = -u0 sin(lon) sin(alpha)
      // u_lon = u_phi   =  u0 (cos(theta) cos(alpha) + sin(theta)cos(lon)sin(alpha) )
      (*utmp[0])[k]  = -u0*sin(lon)*sin(alpha);
      (*utmp[1])[k]  =  u0*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha) );
    }
    GMTK::vsphere2cart(grid, utmp, GVECTYPE_PHYS, c);
//  GMTK::constrain2sphere(grid, c);
  }

  *u[0] = 0.0;
  for ( GSIZET k=0; k<nlumps; k++ ) {

    // Allow different sigma/concentration for each lump:
    sig  [k] = sqrt(sig0[k]*sig0[k] + 4.0*time*nu); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = c0[k]*pow(sig0[k]/sig[k],GDIM);

    // Find where lat/lon endpoint would be at t, if alpha=0:
    latp[k]  = lat0[k];
    lonp[k]  = lon0[k] + (u0/rad) * time;

    // Find where Cart endpoint would be at t, if alpha=0:
    rt[0] = rad*cos(latp[k])*cos(lonp[k]);
    rt[1] = rad*cos(latp[k])*sin(lonp[k]);
    rt[2] = rad*sin(latp[k]);

    // Now, rotate rt about x-axis by alpha to
    // find lat/lon of final position of lump:
    xx[0] = rt[0]; xx[1] = rt[1]; xx[2] = rt[2];
    if ( time > 0 ) {
      xx[1] =  cos(alpha)*rt[1] + sin(alpha)*rt[2];
      xx[2] = -sin(alpha)*rt[1] + cos(alpha)*rt[2];
    }
    latp[k]  = asin(xx[2]/rad);
    lonp[k]  = atan2(xx[1],xx[0]);

    rt[0] = rad*cos(latp[k])*cos(lonp[k]);
    rt[1] = rad*cos(latp[k])*sin(lonp[k]);
    rt[2] = rad*sin(latp[k]);

    for ( GSIZET j=0; j<nxy; j++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c becomes a function of t, this muct change:

      x   = (*xnodes)[0][j];
      y   = (*xnodes)[1][j];
      z   = (*xnodes)[2][j];
#if 1
      r   = sqrt(x*x + y*y + z*z);
      lat = asin(z/r);
      lon = atan2(y,x);
      // Compute arclength from point to where center _should_ be:
      s     = r*acos( sin(latp[k])*sin(lat) + cos(latp[k])*cos(lat)*cos(lon-lonp[k]) );
     (*u[0])[j] += ufact[k]*exp(-s*s*si[k]);
#else
     argxp = pow(x-rt[0],2) + pow(y-rt[1],2) + pow(z-rt[2],2);
     (*u[0])[j] += ufact[k]*exp(-argxp*si[k]);
#endif
    } // end, grid-point loop

  } // end, lump loop

  return TRUE;

} // end of method impl_icosgauss




} // end, ginitstate namespace
