#include "ginits.hpp"


namespace ginits {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on box grids with
//          Dirichlet or periodic boundaries
// ARGS   : stree: state prop tree
//          grid   : grid
//          t    : time
//          utmp : tmp arrays
//          ub   : bdy vectors (one for each state element)
//          u    : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxplaneburgers(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxnwaveburgers: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp;
  GFTYPE           nxy, sig0, u0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = grid.xNodes();
  GTVector<GTVector<GFTYPE>> c(GDIM) 

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_DIRICHLET") >= 2*GDIM
      && "Dirichlet boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  u0    = heatptree.getValue<GFTYPE>("u0");

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<GDIM; j++ ) *c [j] = 0.0;

  if ( bpureadv ) {
     for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];
  }

  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0/sig[k],GDIM);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) changes, change this term accordingly:
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c[i])[j]*t;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*ua[0])[j] = ufact[0]*exp(argxp);
  }

  return TRUE;
} // end, impl_boxnwaveburgers



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxdirgauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          Dirichlet boundaries
// ARGS   : stree: state prop tree
//          grid   : grid
//          t    : time
//          utmp : tmp arrays
//          ub   : bdy vectors (one for each state element)
//          u    : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxdirgauss(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxdirgauss: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp;
  GFTYPE           nxy, sig0, u0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = grid.xNodes();
  GTVector<GTVector<GFTYPE>> c(GDIM) 

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_DIRICHLET") >= 2*GDIM
      && "Dirichlet boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  u0    = heatptree.getValue<GFTYPE>("u0");

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<GDIM; j++ ) *c [j] = 0.0;

  if ( bpureadv ) {
     for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];
  }

  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0/sig[k],GDIM);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) changes, change this term accordingly:
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c[i])[j]*t;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*ua[0])[j] = ufact[0]*exp(argxp);
  }

  return TRUE;
} // end, impl_boxdirgauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxpergauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          periodic boundaries
// ARGS   : stree: state prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays
//          ub   : bdy vectors (one for each state element)
//          u    : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxpergauss((const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxpergauss: ";
  GBOOL            bContin;
  GSIZET           i, j, k, n;
  GFTYPE           iargp, iargm ;
  GFTYPE           isum , irat , prod;
  GFTYPE           sumn , eps;
  GFTYPE           nxy, pint, sig0, u0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = grid.xNodes();
  GTVector<GTVector<GFTYPE>> c(GDIM) 

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = 1.0e-4*std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

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

  // Set adv velocity components. Note:
  // First state elem is the scalar solution, and
  // the remainder are the velocity components:

  for ( j=0; j<GDIM; j++ ) {
    sig[j] = sqrt(sig0*sig0 + 4.0*t*nu_[0]);
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
      f [k]  = modf((*c[k])[j]*t/gL[k],&pint);
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
    (*ua[0])[n] = u0*pow(sig0,GDIM)/pow(sig[0],GDIM)*prod;

  } // end, loop over grid points

  return TRUE;

} // end, impl_boxpergauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosgauss
// DESC   : Initialize state for Burgers with Gauss lump on ICOS grid
// ARGS   : stree: state prop tree
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosgauss(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GString          serr = "impl_icosgauss: ";
  GBOOL            bContin;
  GINT             j, k, n, nlumps;
  GSIZET           nxy;
  GFTYPE           alpha, argxp;
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r, s;
  GFTYPE           rad, u0 ;
  GFTYPE           vtheta, vphi;
  GFTYPE           tiny = std::numeric_limits<GFTYPE>::epsilon();
  GTPoint<GFTYPE>           rt(3);
  GTVector<GFTYPE>          xx(3);
  GTVector<GFTYPE>          si(4), sig(4), ufact(4);
  GTVector<GFTYPE>          latp(4), lonp(4);
  std::vector<GFTYPE>       c0(4), sig0(4);
  std::vector<GFTYPE>       lat0(4), lon0(4); // up to 4 lumps

  PropertyTree lumpptree = ptree.getPropertyTree("init_icosgauss");
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = grid.xNodes();
  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");
  GTVector<GTVector<GFTYPE>> c(GDIM) 

  std::vector<GFTYPE> Omega;

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

  *ua[0] = 0.0;
  for ( GSIZET k=0; k<nlumps; k++ ) {

    // Allow different sigma/concentration for each lump:
    sig  [k] = sqrt(sig0[k]*sig0[k] + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = c0[k]*pow(sig0[k]/sig[k],GDIM);

    // Find where lat/lon endpoint would be at t, if alpha=0:
    latp[k]  = lat0[k];
    lonp[k]  = lon0[k] + (u0/rad) * t;

    // Find where Cart endpoint would be at t, if alpha=0:
    rt[0] = rad*cos(latp[k])*cos(lonp[k]);
    rt[1] = rad*cos(latp[k])*sin(lonp[k]);
    rt[2] = rad*sin(latp[k]);

    // Now, rotate rt about x-axis by alpha to
    // find lat/lon of final position of lump:
    xx[0] = rt[0]; xx[1] = rt[1]; xx[2] = rt[2];
    if ( t > 0 ) {
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
     (*ua[0])[j] += ufact[k]*exp(-s*s*si[k]);
#else
     argxp = pow(x-rt[0],2) + pow(y-rt[1],2) + pow(z-rt[2],2);
     (*ua[0])[j] += ufact[k]*exp(-argxp*si[k]);
#endif
    } // end, grid-point loop

  } // end, lump loop

  return TRUE;

} // end of method impl_icosgauss




} // end, ginits namespace