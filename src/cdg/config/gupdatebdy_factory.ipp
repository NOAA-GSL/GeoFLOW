//==================================================================================
// Module       : gupdatebdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW bdy update object initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : build
// DESC   : Build bdy update object
// ARGS   : ptree  : main property tree
//          grid   : GGrid operator
//          stinfo : StateInfo
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBasePtr
GUpdateBdyFactory<TypePack>::build(const PropertyTree& ptree, Grid &grid, StateInfo &stinfo)
{
  GBOOL         bret = FALSE, use_inits;
  State         uu(u.size());
  GString       sblock, sgrid, supdate;
  PropertyTree  gtree, uptree;
//std::function<void(const geoflow::tbox::PropertyTree& ptree,GString &supdate, Grid &grid,
//                   StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)>
//              mycallback;
  sgrid     = ptree.getValue<GString>("grid_type");
  gtree     = ptree.getPropertyTree(sgrid);
  sblock    = gtree.getValue<GString>("bdy_update_scheme","none");

  UpdateBdyBasePtr base_ptr;
  if      ( "none"         == sblock
    || ""             == sblock ) {
    using UpdateImpl = NullUpdateBdy<TypesPack>

    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = obs_impl;
    return base_ptr;
  }
  else if ( "use_state_init" == sblock ) {
    using UpdateImpl = GFromInitBdy<TypesPack>
    UpdateImpl::Traits traits;
    traits.ptree = ptree;
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = obs_impl;
    return base_ptr;
  }

  
  uptree    = ptree.getPropertyTree(sblock);
  supdate   = uptree.getValue<GString>("update_method");
  if      ( "simple_outflow" == supdate ) {
    using UpdateImpl = GSimpleOutflowBdy<TypesPack>
    UpdateImpl::Traits traits;
    traits.istate.resize(stinfo.nevolve);
    for ( auto j=0; j<traits.istate.size(); j++ ) traits.istate[j] = j;
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( "sponge" == supdate ) {
    using UpdateImpl = GSpongeBdy<TypesPack>
    UpdateImpl::Traits traits;
    traits.istate.resize(stinfo.nevolve);
    for ( auto j=0; j<traits.istate.size(); j++ ) traits.istate[j] = j;
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return base_ptr;

} // end, init method build


//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_from_state
// DESC   : use state var, u, to set bdy, ub
// ARGS   : ptree  : main property tree
//          sconfig: config block name
//          grid   : GGrid operator
//          stinfo : StateInfo
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GUpdateBdyFactory<TypePack>::set_bdy_from_state(const geoflow::tbox::PropertyTree& ptree, GString &sconfig, Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GString       sgrid, supdate;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  bret = GInitStateFactory<TypePack>::init(ptree, grid, stinfo, time, utmp, ub, u);

  // Set from State vector, u and others that we _can_ set:
  for ( auto k=0; k<u.size(); k++ ) {
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_INFLOW].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_INFLOW][j]];
    }
  }

} // end, set_bdy_from_state

