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
// ARGS   : ptree   : main property tree
//          supdate : string naming bdy update prop tree block 
//          grid    : grid object
//          id      : canonical bdy id
//          bdytype : bdy type
//          istate  : state array the method applies to
//          value   : vector of Dirichlet values
//          ibdy    : bdy indirection indices into computational volume
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::UpdateBdyBasePtr
GUpdateBdyFactory<Types>::build(const PropertyTree& ptree, GString &supdate, Grid &grid, const GINT id, GBdyType bdytype, GTVector<GINT> &istate, GTVector<GFTYPE> &value, GTVector<GSIZET> &ibdy, GSIZET igbdy_start)
{
  GBOOL              bret = FALSE;
  UpdateBdyBasePtr   base_ptr;


  if ( "none"         == supdate
    || ""             == supdate ) {
    using UpdateImpl = NullUpdateBdy<Types>;

    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }

  if ( !ptree.isPropertyTree(supdate) ) {
    cout << "GUpdateBdyFactory<Types>::build: PropertyTree " << supdate << " not found" << endl;
    assert(FALSE);
  }

  base_ptr = GUpdateBdyFactory<Types>::get_bdy_class(ptree, supdate, grid, id, bdytype, istate, value, ibdy, igbdy_start);

  return base_ptr;

} // end, init method build


//**********************************************************************************
//**********************************************************************************
// METHOD : get_inflow_callback
// DESC   : Gets CallbackPtr corresponding to sname arg for inflow conditions.
//          the function are gottne from the gns_inflow_user.* namespace
//          collection of methods, which the user may modify.
// ARGS   : sname : inflow function name
//          id    : canonical bdy id
// RETURNS: CallbackPtr for callback function
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::CallbackPtr
GUpdateBdyFactory<Types>::get_inflow_callback(const GString& sname, const GINT id)
{
  GBOOL              bret = FALSE;
  CallbackPtr        callback;


  if      ( ""     == sname 
   ||       "none" == sname ) {
    assert(FALSE); // Must have a named method
  }
  if      ( "myinflow"     == sname ) {
    callback = 

         [](Grid      &grid,
          StateInfo &stinfo,
          Time      &time,
          const GINT id,
          State     &utmp,
          State     &u,
          State     &ub)->GBOOL{return GInflowBdyMethods::myinflow(grid, stinfo, time, id, utmp, u, ub);}; 

  }
  else {
    assert(FALSE && "Specified inflow bdy update method unknown");
  }

  return callback;

} // end, init method get_inflow_callback


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_class
// DESC   : get bdy class for specified bdy type
// ARGS   : ptree   : main prop tree
//          supdate : update block name
//          grid    : grid object
//          id      : canonical bdy id
//          bdytype : bdy type triggering construction of class
//          istate  : vector of state ids
//          value   : vector of Dirichlet values
//          ibdy    : indirection indices into computational volume, 
//                    representing the bdy nodes this method applies to
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::UpdateBdyBasePtr
GUpdateBdyFactory<Types>::get_bdy_class(const PropertyTree& ptree, GString &supdate, Grid &grid, const GINT id, const GBdyType bdytype, GTVector<GINT> &istate, GTVector<GFTYPE> &value, GTVector<GSIZET> &ibdy, GSIZET igbdy_start)
{
  GINT               nstate;
  GLONG              iloc;
  GString            sblock ;
  UpdateBdyBasePtr   base_ptr;
  GTVector<GSIZET>  *igbdy = &grid.igbdy();
  PropertyTree       sptree; // update block tree
  

  nstate = istate.size();
  
  sptree = ptree.getPropertyTree(supdate);

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<Types>;
    typename GDirichletBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.value .resize(nstate); traits.value  = value;
    traits.bdyid = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }

    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_INFLOW == bdytype ) {
    using UpdateImpl = GInflowBdy<Types>;
    typename GInflowBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    assert( sptree.isValue<GBOOL>("use_init") && "use_init boolean missing");
    traits.use_init = sptree.getValue<GBOOL>("use_init");
    if ( !traits.use_init ) {
      assert( sptree.isValue<GString>("inflow_method") 
           && "inflow_method required if use_init==FALSE" ); 
      sblock = sptree.getValue<GString>("inflow_method");
      traits.callback = GUpdateBdyFactory<Types>::get_inflow_callback(sblock, id);
    }
    traits.ptree = ptree;
    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<Types>;
    typename GNoSlipBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<Types>;
    typename G0FluxBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol.resize(ibdy.size()); traits.ibdyvol = ibdy;
    traits.ibdyloc.resize(traits.ibdyvol.size());
    // Find index of ibdyvol in global bdy vector:
    for ( auto j=0; j<traits.ibdyloc.size(); j++ ) {
      iloc = j + igbdy_start;
      assert(iloc >= 0);
      traits.ibdyloc[j] = iloc;
    }
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_OUTFLOW == bdytype ) {
    assert(FALSE); // not available yet
    using UpdateImpl = GOutflowBdy<Types>;
    typename GOutflowBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<Types>;
    typename GSpongeBdy<Types>::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }

    assert(sptree.isArray<Ftype>("farfield") && "farfield array missing");
    traits.idir = sptree.getValue<Ftype>("farfield");
    assert(sptree.isArray<Ftype>("exponent") && "exponent array missing");
    traits.idir = sptree.getValue<Ftype>("exponent");
    assert(sptree.isArray<Ftype>("sigma0") && "sigma0 array missing");
    traits.idir = sptree.getValue<Ftype>("sigma0");
    assert(sptree.isArray<Ftype>("xstart") && "xstart array missing");
    traits.idir = sptree.getValue<Ftype>("xstart");
    
    GTVector<Ftype> xmin, xmax;
    geoflow::coord_dims(ptree, xmin, xmax);// get coord min/max from ptree

    assert(sptree.isValue<GINT>("idir") && "idir value missing");
    traits.idir = sptree.getValue<GINT>("idir");
    assert(traits.idir >= 1 && traits.idir <= xmax.size()); // validate idir

    if ( traits.idir < 0 ) traits.ro = xmin[abs(traits.idir)-1];
    else                   traits.ro = xmax[abs(traits.idir)-1];

    // Allocate update Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return base_ptr;

} // end, init method get_bdy_class


