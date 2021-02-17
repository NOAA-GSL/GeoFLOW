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
#include "gmtk.hpp"
#include "ginitstate_direct_user.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include <random>


namespace ginitstate {



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on box grids with
//          Dirichlet or periodic boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxnwaveburgers(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxnwaveburgers: ";
  GBOOL            bret = TRUE, brot = FALSE;
  GINT             nlump=0;
  GSIZET           i, j, nxy;
  GFTYPE           K2, nu, Re, r2, tdenom;
  GFTYPE           efact, sum, tfact, tt, xfact;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), t0;
  GTPoint<GFTYPE>  kprop(3), r0(3), P0(3), gL(3);
  std::vector<GFTYPE>  kxprop, kyprop, kzprop;
  std::vector<GFTYPE>  xinit , yinit , zinit ;
  std::vector<GBOOL>   bplanar;
  std::vector<GFTYPE>  tinit;
  std::vector<GFTYPE>  ULparam;
  GString              snut;

  PropertyTree nwaveptree = ptree   .getPropertyTree(sconfig);
  PropertyTree boxptree   = ptree   .getPropertyTree("grid_box");
  PropertyTree nuptree    = ptree.getPropertyTree("dissipation_traits");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR 
      || grid.gtype() == GE_DEFORMED  && "Invalid element types");

  tt = time;

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
//    && "INFLOW boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling. But we won't parameterize in terms of Re, 
  // but rather, nu.
  // Get some parameters; xinit, tinit, ULparam, bplanar, kprop,
  // should have the same number, nlump, elements, one foreach 'wave':
  xinit      = nwaveptree.getArray<GFTYPE>("x0"); 
  yinit      = nwaveptree.getArray<GFTYPE>("y0"); 
  zinit      = nwaveptree.getArray<GFTYPE>("z0"); 
  tinit      = nwaveptree.getArray<GFTYPE>("t0");
  ULparam    = nwaveptree.getArray<GFTYPE>("ULparam");
  bplanar    = nwaveptree.getArray<GBOOL> ("planar");
  kxprop = nwaveptree.getArray<GFTYPE>("prop_dir_x");
  kyprop = nwaveptree.getArray<GFTYPE>("prop_dir_y");
  kzprop = nwaveptree.getArray<GFTYPE>("prop_dir_z");
//Re         = nwaveptree.getValue<GFTYPE>("Re",6.0);

  nu       = nuptree   .getValue<GFTYPE>("nu",0.0833);
  snut     = nuptree   .getValue<GString>("nu_type","constant");

  t0 = tinit;

  assert(yinit  .size() == xinit.size()
      && tinit  .size() == xinit.size()
      && ULparam.size() == xinit.size()
      && bplanar.size() == xinit.size()
      && kxprop .size() == xinit.size()
      && kyprop .size() == xinit.size()
      && "(1)Lump count must be consistent");

  if ( GDIM > 2 ) {
    assert(zinit  .size() == xinit.size()
        && kzprop .size() == xinit.size()
        && "(2)Lump count must be consistent");
  }
  nlump = xinit.size();

  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  if ( (time-t0.max()) <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = t0.max();


  for ( i=0; i<GDIM; i++ ) {
    *u[i] = 0.0;
  }

  tdenom  = 1.0/(4.0*nu*time);
  for ( auto ilump=0; ilump<nlump; ilump++ ) {
    r0[0]  = xinit[ilump]; r0[1]  = yinit[ilump]; 
    if ( GDIM > 2 ) r0[2]  = zinit[ilump]; 
    kprop[0] = kxprop[ilump]; kprop[1] = kyprop[ilump];
    if ( GDIM > 2 ) kprop[2]  = kzprop[ilump]; 
    Re  = ULparam[ilump]/nu; // set Re from nu and ULparam
//   cout << "impl_boxnwaveburgers: ilump=" << ilump << " nu=" << nu << " Re=" << Re << " tt=" << tt << " time=" << time << endl;
    for ( i=0, K2=0.0; i<GDIM; i++ ) K2 += kprop[i]*kprop[i];
    assert(K2 > 10.0*std::numeric_limits<GFTYPE>::epsilon() && "Prop direction, kprop, not set");
    if ( bplanar[ilump] ) kprop  *= 1.0/kprop.norm();
    tfact   = bplanar[ilump] ? sqrt(time/t0[ilump]): time/t0[ilump];
    // If prop direction has more than one component != 0. Then
    // front is rotated (but still planar):
//  for ( i=0, brot=TRUE; i<GDIM; i++ ) brot = brot && K[i] != 0.0 ;
//  for ( i=0, idir=0; i<GDIM; i++ ) if ( K[i] > 0 ) {idir=i; break;}
//  K2 = brot && K2 == 0.0 ? 1.0 : K2;
//  if ( time <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = K2 * t0;
    for ( j=0; j<nxy; j++ ) {
      for ( i=0; i<GDIM; i++ ) {
        xx[i] = (*xnodes)[i][j] - r0[i];
      }
      if ( bplanar[ilump] ) { // compute k.r for planar wave
        for ( i=0, sum=0.0; i<GDIM; i++ ) { 
          sum += kprop[i]*xx[i];
          xx[i] = 0.0;
        }
        xx[0] = sum;
      }
      for ( i=0, r2=0.0; i<GDIM; i++ ) r2 += xx[i]*xx[i];  
      efact   = tfact * exp(r2*tdenom) / ( exp(Re) - 1.0 );
      xfact   = 1.0 /( time * (  1.0 + efact ) );
      // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t)) ]^-1
      for ( i=0; i<GDIM; i++ ) (*u[i])[j] += xx[i]*xfact;
  //cout << "impl_boxnwaveburgers: ux[" << j << "]=" << (*u[0])[j] << endl;
      // dU1max = 1.0 / ( time * (sqrt(time/A) + 1.0) );
      // aArea  = 4.0*nu*log( 1.0 + sqrt(A/time) );
    } // end, coord loop
  } // end, lump loop

  return bret;

} // end, impl_boxnwaveburgers


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on icos grids
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosnwaveburgers(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_icosnwaveburgers: ";
  GBOOL            bret;
  GSIZET           i, j, nxy;
  GFTYPE           nu, Re, r, s, tdenom;
  GFTYPE           lat, lon;
  GFTYPE           efact, sum, tfact, xfact;
  GFTYPE           x, y, z;
  GTVector<GFTYPE>            t0, xx(GDIM+1);
  GTVector<GTPoint<GFTYPE>>   r0(GDIM+1);
  std::vector<GFTYPE>         lat0, lon0, st0, Uparam;

  PropertyTree nwaveptree = ptree   .getPropertyTree(sconfig);
  PropertyTree gridptree  = ptree   .getPropertyTree("grid_icos");
  PropertyTree nuptree    = ptree.getPropertyTree("dissipation_traits");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");
  assert(u.size() >= GDIM+1 && "Insufficient number of state members");


  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling. But we won't parameterize in terms of Re, 
  // but rather, nu.
  // Set some parameters:
  r      = gridptree.getValue <GFTYPE>("radius");
  lat0   = nwaveptree.getArray<GFTYPE>("latitude0"); 
  lon0   = nwaveptree.getArray<GFTYPE>("longitude0"); 
  Uparam = nwaveptree.getArray<GFTYPE>("Uparam");
//Re     = nwaveptree.getValue<GFTYPE>("Re",6.0);
  st0    = nwaveptree.getArray<GFTYPE>("t0");
  nu     = nuptree   .getValue<GFTYPE>("nu",0.0833);

  t0.resize(st0.size());
  t0     = st0;
  assert(lat0.size() == lon0.size() 
      && lat0.size() == t0.size()
      && "lat0, lon0, and t0 must be the same size");

  if ( (time-t0.max()) <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = t0.max();

  // Convert initial positions to radians:
  for ( GINT ilump=0; ilump<lat0.size(); ilump++) {
    lat0[ilump]    *= (PI/180.0);
    lon0[ilump]    *= (PI/180.0);
    r0  [ilump].x1  = r*cos(lat0[ilump])*cos(lon0[ilump]);
    r0  [ilump].x2  = r*cos(lat0[ilump])*sin(lon0[ilump]);
    r0  [ilump].x3  = r*sin(lat0[ilump]);
  }

  for ( i=0; i<GDIM+1; i++ ) *u[i] = 0.0;
   
  tdenom = 1.0/(4.0*nu*time);
  // Initialize each lump:
  for ( GINT ilump=0; ilump<lat0.size(); ilump++) {
    Re = Uparam[ilump]*r/nu; // set Re from nu, U, radius
    for ( j=0; j<nxy; j++ ) {
       x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
       lat = asin(z/r);
       lon = atan2(y,x);
      for ( i=0; i<GDIM+1; i++ ) { 
        xx[i] = (*xnodes)[i][j] - r0[ilump][i];
      }
  
      xx[0] = r*lat;
      xx[1] = r*lon;
  
      // find arclength from lump center
      s     = r*acos( sin(lat)*sin(lat0[ilump])
            +       cos(lat)*cos(lat0[ilump])*cos(lon - lon0[ilump]) );
  
      tfact  = time/t0[ilump];
      efact  = tfact * exp(s*s*tdenom) / ( exp(Re) - 1.0 );
      xfact  = 1.0 /( time * (  1.0 + efact ) );
      for ( i=0; i<GDIM+1; i++ ) {
        (*u[i])[j] += xx[i]*xfact;
        assert( std::isfinite((*u[i])[j]) );
      }
      // dU1max = 1.0 / ( time * (sqrt(time/A) + 1.0) );
      // aArea  = 4.0*nu*log( 1.0 + sqrt(A/time) );
    } // end, coord j-loop 
  } // end, ilump-loop

//GMTK::vsphere2cart(grid, usph, GVECTYPE_PHYS, u);
  GMTK::constrain2sphere(grid, u);

  bret = TRUE;
  for ( j=0; j<GDIM+1; j++ ) {
     bret = bret && u[j]->isfinite(i);
  }

  assert(bret && "Initial conditions not finite!");


  return bret;

} // end, impl_icosnwaveburgers


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxdirgauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          Dirichlet boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          time   : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxdirgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxdirgauss: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp;
  GFTYPE           nxy, nu, sig0, E0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  State            c(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);
  GString          snut;

cout << "impl_boxdirgauss: sconfig=" << sconfig << endl;

  PropertyTree initptree = ptree.getPropertyTree(sconfig);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = initptree.getArray<GFTYPE>("adv_vel");
  }

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[j+1];

  // Check bdy conditions:
  GTVector<GTVector<GBdyType>>
                           *igbdyt_face= &grid.igbdyt_bdyface();
  for ( auto j=0; j<igbdyt_face->size(); j++ ) { // for each face
cout << "boxpergauss: num=" << (*igbdyt_face)[j].size() << " igbdyt_face[" << j << "]=" << (*igbdyt_face)[j] << endl;
    assert( (*igbdyt_face)[j].onlycontains(GBDY_INFLOW) 
        &&  "Inflow conditions must be set on all boundaries");
  }

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = initptree.getValue<GFTYPE>("x0");
  r0.x2 = initptree.getValue<GFTYPE>("y0");
  r0.x3 = initptree.getValue<GFTYPE>("z0");
  sig0  = initptree.getValue<GFTYPE>("sigma");
  E0    = initptree.getValue<GFTYPE>("E0");

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
    ufact[k] = sqrt(2*E0)*pow(sig0/sig[k],GDIM);
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
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxpergauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxpergauss: ";
  GBOOL            bContin;
  GSIZET           i, j, k, n;
  GFTYPE           iargp, iargm ;
  GFTYPE           isum , irat , prod;
  GFTYPE           sumn , eps;
  GFTYPE           nxy, nu, pint, sig0, E0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);
  State            c(GDIM);
  GString          snut;

  PropertyTree initptree = ptree.getPropertyTree(sconfig);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");
  
  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = 1.0e4*std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[j+1];

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = initptree.getArray<GFTYPE>("adv_vel");
  }

  GTVector<GTVector<GBdyType>>
                           *igbdyt_face= &grid.igbdyt_bdyface();
  for ( auto j=0; j<igbdyt_face->size(); j++ ) { // for each face
cout << "boxpergauss: num=" << (*igbdyt_face)[j].size() << " igbdyt_face[" << j << "]=" << (*igbdyt_face)[j] << endl;
    assert( (*igbdyt_face)[j].onlycontains(GBDY_PERIODIC) 
        &&  "Periodic boundaries must be set on all boundaries");
  }

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = initptree.getValue<GFTYPE>("x0");
  r0.x2 = initptree.getValue<GFTYPE>("y0");
  r0.x3 = initptree.getValue<GFTYPE>("z0");
  sig0  = initptree.getValue<GFTYPE>("sigma");
  E0    = initptree.getValue<GFTYPE>("E0");

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
//    f [k]  = modf((*c[k])[j]*time/gL[k],&pint);
      f [k]  = (*c[k])[n]*time/gL[k];
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
    (*u[0])[n] = sqrt(2.0*E0)*pow(sig0,GDIM)/pow(sig[0],GDIM)*prod;

  } // end, loop over grid points

  return TRUE;

} // end, impl_boxpergauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosgauss
// DESC   : Initialize state for Burgers with Gauss lump on ICOS grid, based
//          on a the Cartesian method
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{

  GString             serr = "impl_icosgauss: ";
  GINT                nlumps;
  GSIZET              nxy;
  GFTYPE              alpha, lat, latc, lon, lonc;
  GFTYPE              x, y, z, r;
  GFTYPE              irad, nu, rad, rexcl;
  GFTYPE              cexcl, num, den;
  GFTYPE              c0, cosk, spc;
  GTVector<GFTYPE>    isig;
  GTVector<GTVector<GFTYPE>*>
                      c(3);
  std::vector<GFTYPE> u0, sig0;
  std::vector<GFTYPE> lat0, lon0; 
  GString             snut;

  PropertyTree lumpptree = ptree.getPropertyTree(sconfig);
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");

  assert(utmp.size() >= 7 );
  for (auto j=0; j<c.size(); j++ ) c[j] = u[j+1]; // adv vel. comp.


  nxy = (*xnodes)[0].size(); // same size for x, y, z

  lat0  = lumpptree.getArray<GFTYPE>("latitude0");       // lat for each lump
  lon0  = lumpptree.getArray<GFTYPE>("longitude0");      // lon for each lump
  sig0  = lumpptree.getArray<GFTYPE>("sigma");           // sig for each lump
  u0    = lumpptree.getArray<GFTYPE>("u0");              // initial concentrations for each lump
  c0    = lumpptree.getValue<GFTYPE>("c0");              // adv. vel. magnitude
  alpha = lumpptree.getValue<GFTYPE>("alpha");           // initial concentrations for each lump
  cexcl = lumpptree.getValue<GFTYPE>("excl_angle");      // exclusion angle (degrees)
  rad   = icosptree.getValue<GFTYPE>("radius");
  irad  = 1.0/rad;

  nlumps = lat0.size();
  rexcl = cexcl * PI/180.0;

  nu    = nuptree   .getValue<GFTYPE>("nu");
  snut  = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // Convert initial locations from degrees to radians,
  // & compute initial positions of lumps in Cart coords:
  isig.resize(sig0.size());
  for ( auto k=0; k<nlumps; k++ ) {
    lat0[k] *= PI/180.0;
    lon0[k] *= PI/180.0;
    isig[k]  = 1.0/sqrt(sig0[k]*sig0[k] + 4.0*nu*time);
  }
  alpha *= PI/180.0;


  // We use (lat,lon) advecting components from 
  // Williamson JCP 102:211 (1992):
  //   v_lat = -c0 sin(lon)sin(alpha)
  //   v_lon =  c0 ( cos(lat) cos(alpha) + sin(lat) sin(lon) sin(alpha)).
  // We compute the Cartesian components required by the solver:
  for ( auto j=0; j<nxy; j++ ) { 
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    r = sqrt(x*x + y*y + z*z); lat = asin(z/r); lon = atan2(y,x);
    (*utmp[0])[j] = -c0 * sin(lon) * sin(alpha);
    (*utmp[1])[j] =  c0 * ( cos(lat) * cos(alpha) - sin(lat)*cos(lon)*sin(alpha) );
  }
  GMTK::vsphere2cart(grid, utmp, GVECTYPE_PHYS, c);

  *u[0] = 0.0;
  for ( auto k=0; k<nlumps; k++ ) {

        for ( auto j=0; j<nxy; j++ ) { 
          x     = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
          r     = sqrt(x*x + y*y + z*z); lat = asin(z/r); lon = atan2(y,x);

          if ( (0.5*PI-fabs(lat)) < rexcl ) continue;
          // Define following arclengths:
          // s0c:  path from (lat0,lon0) to where solid body rotation
          //       takes peak
          // s0p:  path from (lat0,lon0) to arbitrary point
          // spc:  path from arbitrary point to where solid body rotation
          //       takes peak
         
          // Compute length from initial point to arb point on grid:
//        s0p   = r*acos( sin(lat0[k])*sin(lat) + cos(lat0[k])*cos(lat)*cos(lon-lon0[k]) );
          // Compute where lump should be (convert C to contravar. components):
          latc  = lat0[k] - c0 * irad * sin(lon0[k]) * sin(alpha) * time;
          lonc  = lon0[k] + c0 * irad * ( cos(alpha) + tan(lat0[k]) * sin(lon0[k]) * sin(alpha) ) * time;
          spc   = r*acos( sin(lat)*sin(latc) + cos(lat)*cos(latc)*cos(lonc-lon) );

          // Compute solution along trajectory:
          (*u[0])[j] += u0[k]*pow(sig0[k]*isig[k],GDIM) 
                      * exp( -pow(spc*isig[k], 2.0) ); 
        }

  } // end, lump loop

  return TRUE;

} // end of method impl_icosgauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxdrybubble
// DESC   : Initialize state for GMConv solver with cold/dry bubble
//          on box grids. Taken from Straka et al. 1993, Int. J. Num. 
//          Meth. Fluids 17:1, and from Bryan & Fritsch 2002 MWR
//          130:2817.
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxdrybubble(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{

  GString             serr = "impl_boxdrybubble: ";
  GSIZET              nxy;
  GFTYPE              x, y, z, r;
  GFTYPE              delT, dj, exnerb, exner, L, P0, pj, thetab, T0, Ts;
  GTVector<GFTYPE>   *db, *d, *e, *pb, *Tb;
  std::vector<GFTYPE> xc, xr;  
  GString             sblock;

  PropertyTree inittree   = ptree.getPropertyTree(sconfig);
  sblock                   = ptree.getValue<GString>("pde_name");
  PropertyTree convptree   = ptree.getPropertyTree(sblock);

  GGridBox  *box   = dynamic_cast <GGridBox*>(&grid);
  assert(box && "Must use a box grid");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(u.size() == GDIM+4);

  Tb    = utmp[0];  // background temp
  e     = u  [GDIM];// int. energy density
  d     = u[GDIM+1];// total density fluctuation
  db    = u[GDIM+2];// background density fluct, from solver
  pb    = u[GDIM+3];// background pressure  , from solver
  nxy   = (*xnodes)[0].size(); // same size for x, y, z

  T0    = inittree.getValue<GFTYPE>("T_pert", 15.0);    // temp. perturb. magnitude (K)
  xc    = inittree.getArray<GFTYPE>("x_center");        // center location
  xr    = inittree.getArray<GFTYPE>("x_width");         // bubble width
  P0    = convptree.getValue<GFTYPE>("P0");              // ref pressure (mb or hPa)
  P0   *= 100.0;                                         // convert to Pa
  Ts    = convptree.getValue<GFTYPE>("T_surf");          // surf temp

  assert(xc.size() >= GDIM && xr.size() >= GDIM);

 *u[0]  = 0.0; // sx
 *u[1]  = 0.0; // sy
 if ( GDIM == 3 ) *u[2]  = 0.0; // sz
 *d     = 0.0;

  for ( auto j=0; j<nxy; j++ ) { 
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; 
    if ( GDIM == 3 ) z = (*xnodes)[2][j];
    r = GDIM == 3 ? z : y;
    L         = pow((x-xc[0])/xr[0],2) + pow((y-xc[1])/xr[1],2);
    L        += GDIM == 3 ? pow((z-xc[2])/xr[2],2) : 0.0;
    L         = sqrt(L);
    exnerb    = pow((*pb)[j]/P0, RD/CPD);
    delT      = L <= 1.0 ? 2.0*T0*pow(cos(0.5*PI*L),2.0) : 0.0;
#if 1
    // Ts, delT are pot'l temp, 
    (*Tb)[j]  = (Ts + delT)*exnerb; // T = (theta + dtheta)*exner
    pj        = (*pb)[j]; 
    (*d)[j]   = pj / ( RD * (*Tb)[j]  ) - (*db)[j];
    dj        = (*d)[j] + (*db)[j];
    (*e)[j]   = CVD * dj * ( (*Tb)[j] ); // e = Cv d (T+delT);
//  (*Tb)[j]  = (Ts + delT)*exnerb; // T = (theta + dtheta)*exner
//  (*Tb)[j]  = (thetab + delT)*exnerb;
//  (*d)[j]   = pj / ( RD * (*Tb)[j] )  - (*db)[j];
//  (*e)[j]   = CVD * dj * (thetab+delT)*(exnerb); // e = Cv d (theta+dtheta) * exner;

#else
    // Check that hydrostatic state is maintained:
    (*d)[j]   = 0.0;
    (*Tb)[j]  = Ts*exnerb; // T = (theta + dtheta)*exner
    (*e)[j]   = CVD * (*db)[j] * ( (*Tb)[j] ); // e = Cv d (T);
#endif

  }
//cout << "boxdrybubble: db=" << *db << endl;

  return TRUE;

} // end of method impl_boxdrybubble


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosabcconv
// DESC   : Initialize state for GMConv solver with ABC initial conditions
//          on velocity, and sinusoidal function in theta
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosabcconv(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{

  GString             serr = "impl_icosabcconv: ";
  GSIZET              kdn, kup, nc, nxy;
  GFTYPE              alpha, A, B, C, fact,hphase, poly;
  GFTYPE              x, y, z, r, ri, ro, lat, lon;
  GFTYPE              exner, p, pi2, P0, T0;
  GTVector<GFTYPE>   *d, *e;
  GTVector<GTVector<GFTYPE>*>
                      vh(2);
  GString             sblock;
  std::default_random_engine        generator;
  std::normal_distribution<GFTYPE> *distribution=NULLPTR;


  PropertyTree inittree    = ptree.getPropertyTree(sconfig);
  sblock                   = ptree.getValue<GString>("pde_name");
  PropertyTree convptree   = ptree.getPropertyTree(sblock);
  sblock                   = ptree.getValue<GString>("grid_type");
  PropertyTree gridptree   = ptree.getPropertyTree(sblock);


  GGridIcos  *icos   = dynamic_cast <GGridIcos*>(&grid);
  assert(icos && "Must use ICOS grid");

  nc = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(u.size() >= nc+1);

  distribution = new normal_distribution<GFTYPE>(0,2.0*PI);
  e     = u  [nc];// int. energy density
  d     = u[nc+1];// total density 
  nxy   = (*xnodes)[0].size(); // same size for x, y, z

  #if defined(_G_IS2D)
  ro     = gridptree.getValue <GFTYPE>("radius");
  ri     = 0.0;
  #elif defined(_G_IS3D)
  ri     = gridptree.getValue <GFTYPE>("radiusi");
  ro     = gridptree.getValue <GFTYPE>("radiuso");
  #endif
  hphase = 2.0*PI/(ro - ri);
  T0    = inittree.getValue<GFTYPE>("T0", 15.0);  
  P0    = inittree.getValue<GFTYPE>("P0", 15.0); // in mb 
        P0 *= 100.0; // convert to Pa
  A     = inittree.getValue<GFTYPE>("A", 0.9);  
  B     = inittree.getValue<GFTYPE>("B", 1.0);
  C     = inittree.getValue<GFTYPE>("C", 1.1);
  poly  = inittree.getValue<GFTYPE>("poly", 0.0); 
  kdn   = inittree.getValue  <GINT>("kdn", 1); 
  kup   = inittree.getValue  <GINT>("kup", 10); 

  vh[0] = utmp[0];
  vh[1] = utmp[1];

 *u[0]  = 0.0; // sx
 *u[1]  = 0.0; // sy
 *u[2]  = 0.0; // sz
 *d     = 0.0;

  // Initialize each variable:
  for ( auto k=kdn; k<kup; k++ ) { // sum over wavemodes
    alpha =  hphase*(r - ri)*(*distribution)(generator);
    for ( auto j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r   = sqrt(x*x + y*y + z*z);
      lat = asin(z/r);
      lon = atan2(y,x);
      (*d)   [j] +=  fabs( B*cos(pi2*y+alpha) + C*sin(pi2*z+alpha) ) / pow(k,poly) + 0.001;
      pi2         = 2.0*PI*k;
      fact         = exp(-32.0*lat*lat/(PI*PI));
#if 0
      (*u[0])[j] +=  ( B*cos(k*lon) + C*sin(k*lat) ) / pow(k,poly);
      (*u[1])[j] +=  ( A*sin(k*lon) + C*cos(k*lat) ) / pow(k,poly);
      (*u[2])[j] +=  ( A*cos(k*lon) + B*sin(k*lat) ) / pow(k,poly);
#else
      (*vh[0])[j] +=  fact*A*cos    (k*lat+alpha) / pow(k,poly); // lat
      (*vh[1])[j] +=  fact*B*sin(2.0*k*lat+alpha) / pow(k,poly); // long
#endif
      p           = (*d)[j] * RD * T0;
      exner       = pow(p/P0, RD/CPD);
      (*e)[j]     = CVD * (*d)[j] * T0*exner; 
    }
  } // end, coord j-loop 

  // Convert from 2d surface to 3d Cartesian 
  // momentum densities:
  GMTK::vsphere2cart(grid, vh, GVECTYPE_PHYS, u);
  for ( auto j=0; j<3; j++ ) *u[j] *=  *d;

  if ( distribution != NULLPTR ) delete distribution;

  return TRUE;

} // end of method impl_icosabcconv


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxsod
// DESC   : Initialize state for GMConv solver with Sod shock tube.
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          stinfo : StateInfo
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxsod(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u)
{

  GString             serr = "impl_boxsod: ";
  GSIZET              nxy;
  GFTYPE              a, b;
  GFTYPE              x, y, z;
  GFTYPE              Pfact, P0, pj, T0, width, xc;
  GTVector<GFTYPE>   *d, *e;

  PropertyTree sodptree   = ptree.getPropertyTree(sconfig);

  GGridBox  *box   = dynamic_cast <GGridBox*>(&grid);
  assert(box && "Must use a box grid");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(u.size() == 4);

  e     = u  [GDIM];// int. energy density
  d     = u[GDIM+1];// total density 
  nxy   = (*xnodes)[0].size(); // same size for x, y, z

  T0    = sodptree.getValue<GFTYPE>("T0",300.0);   // ref. temp
  P0    = sodptree.getValue<GFTYPE>("P0",1000.0);  // ref. pressure 
  Pfact = sodptree.getValue<GFTYPE>("Pfact",10.0); // Pleft/Pright
  xc    = sodptree.getValue<GFTYPE>("x_center");   // center location
  width = sodptree.getValue<GFTYPE>("width");      // shock width
  P0   *= 1.0e2;  // convert P0 from mb to Pa

 *u[0]  = 0.0; // sx
 *u[1]  = 0.0; // sy
 if ( GDIM == 3 ) *u[2]  = 0.0; // sz
 *d     = 0.0;

  a = P0*(1.0-Pfact)/PI;
  b = P0*(1.0+Pfact)/2.0;
  for ( auto j=0; j<nxy; j++ ) { 
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; 
    if ( GDIM == 3 ) z = (*xnodes)[2][j];
    pj = a*atan((x-xc)/width) + b;
       
//  (*d)[j]   = (*pb)[j] / ( RD * ( (*Tb)[j] + delT ) ) - (*db)[j];
    (*d)[j]   = pj / ( RD * T0 );
    (*e) [j]  = CVD * (*d)[j]  * T0;
cout << "boxsod: p=" << pj << " d=" << (*d)[j] <<  endl;

  }

  return TRUE;

} // end of method impl_boxsod




} // end, ginitstate namespace
