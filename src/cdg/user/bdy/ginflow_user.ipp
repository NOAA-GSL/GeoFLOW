//==================================================================================
// Module       : ginitbdy_user.cpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations provided
//                by user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================





//**********************************************************************************
//**********************************************************************************
// METHOD : myinflow
// DESC   : Place holder template for inflow update methods
// ARGS   : 
//          ptree  : main PropertyTree
//          sconfig: configuration block name in ptre
//          eqn    : equation implementation
//          grid   : grid
//          t      : time
//          id     : canonical bdy id
//          utmp   : tmp arrays
//          u      : current state, overwritten here
//          ub     : bdy vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
template<typename Types>
GBOOL GInflowUser<Types>::myinflow(const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, const GINT id, State &utmp, State &u, State &ub)
{
  GSIZET              j, nxy;
  GFTYPE              x, y, z, r;
  GFTYPE              gam, igam;
  GFTYPE              dj, ds, N, P0, pj, pt, th, T0, Tb, Ts, U0;
  GFTYPE              eps=100.0*std::numeric_limits<GFTYPE>::epsilon();
  GTVector<GFTYPE>   *db, *d, *e, *pb, *T;
  GTVector<GSIZET>   *igbdy = &(grid.igbdy_bdyface()[id]);

  GString             sblock;
  typename Types::State
                     *ubase;
 
  GTVector<GTVector<GFTYPE>> 
                     *xnodes = &grid.xNodes();
  GTVector<GTVector<GFTYPE>> 
                     *xb     = &grid.xb();
  GMConv<Types>      *ceqn;

  PropertyTree inittree    = ptree.getPropertyTree(sconfig);
  sblock                   = ptree.getValue<GString>("pde_name");
  PropertyTree convptree   = ptree.getPropertyTree(sblock);


  // Check solver type 
  // Remember: eqn is a shared_ptr, so must check 
  //           against its contents
  
  ceqn = dynamic_cast<GMConv<Types>*>(eqn.get());
  assert(ceqn != NULLPTR && "Must initialize for Equation GMConv");

  // Check grid type:
  GridBox  *box   = dynamic_cast <GridBox*>(&grid);
  assert(box && "Must use a box grid");

  // Check state size:
  assert(u.size() == ceqn->state_size());

  // Check tmp size:
  assert(utmp.size() >= 1 );

  // Get base state:
  typename GMConv<Types>::Traits traits = ceqn->get_traits();
  ubase = &ceqn->get_base_state();
//assert( traits.usebase && ubase->size() == 2 );
  assert( ubase->size() == 2 );
   

  T     = utmp[0];  // background temp
  d     = ub [ceqn->DENSITY]; // density
  e     = ub [ceqn->ENERGY]; // int. energy density
  db    = (*ubase)[0];// background density 
  pb    = (*ubase)[1];// background pressure

  N     = inittree.getValue<GFTYPE>("N");          // Brunt-Vaisalla freq
  U0    = inittree.getValue<GFTYPE>("U0");         // Inflow velocity
  P0    = convptree.getValue<GFTYPE>("P0");        // ref pressure (mb or hPa)
  P0   *= 100.0;                                   // convert to Pa
  Ts    = convptree.getValue<GFTYPE>("T_surf");    // surf temp

  // Initialize momentum:
  for ( auto j=0; j<ceqn->ENERGY; j++ ) *ub[j] = 0.0;
//ds = P0 / (RD * Ts); // surf. density
  gam  = RD / CPD;
  igam = 1.0/gam;
  for ( auto jj=0; jj<igbdy->size(); jj++ ) {
    j = (*igbdy)[jj];
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; 
    if ( GDIM == 3 ) z = (*xnodes)[2][j];
    r = GDIM == 3 ? z : y;

    // Compute den from constant Brunt-Vaisalla freq,
    //    N^2 = -g/th_0 d theta/dz = const:
    th = Ts * ( N*N*r/GG + 1.0);
    pj = P0 * pow(th/Ts, -igam);
    dj = pj / ( RD * Ts );
    if ( traits.usebase ) { // There is a base-state
      (*d) [jj]  = dj - (*db)[j];
    }
    else {                  // No base-state
      (*d) [jj]  = dj;
   }
   (*ub[0])[jj] = dj * U0;
   (*e)[jj]    = CVD * pj / RD; // e = Cv * p / R
   
  }

  return TRUE;

} // end of method myinflow


