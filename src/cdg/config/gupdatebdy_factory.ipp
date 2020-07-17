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
//          ibdy    : bdy indirection indices into computational volume
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBdyBaseList
GUpdateBdyFactory<TypePack>::build(const PropertyTree& ptree, GString &supdate, Grid &grid, const GINT id, GBdyType bdytype, GTVector<GINT> &istate, GTVector<GSIZET> &ibdy)
{
  GBOOL              bret = FALSE;
  PropertyTree       sptree; // update block tree
  UpdateBdyBasePtr   base_ptr;


  if ( "none"         == supdate
    || ""             == supdate ) {
    using UpdateImpl = NullUpdateBdy<TypesPack>

    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = obs_impl;
    list.push_back(base_ptr);
    return list;
  }

  if ( !ptree.isPropertyTree(supdate) ) {
    cout << "GUpdateBdyFactory<TypePack>::build: PropertyTree " << supdate << " not found" << endl;
    assert(FALSE);
  }
  sptree = ptree.getPropertyTree(supdate);

  base_ptr = get_bdy_class(sptree, grid, id, bdytype, ibdy);

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
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::CallbackPtr
GUpdateBdyFactory<TypePack>::get_inflow_callback(const GString& sname, const GINT id)
{
  GBOOL              bret = FALSE;
  CallbackPtr        callback;


  if      ( ""     == sname 
   ||       "none" == sname ) {
    assert(FALSE); // Must have a named method
  }
  if      ( "myinflow"     == sname ) {
    callback = 
      [GInflowBdyMethods](Grid      &grid,
                          StateInfo &stinfo,
                          Time      &time,
                          State     &utmp,
                          State     &u,
                          State     &ub){myinflow(grid, stinfo, time, id, utmp, u, ub);}; 
  else {
    assert(FALSE && "Specified inflow bdy update method unknown");
  }

  return callback;

} // end, init method get_inflow_callback


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_class
// DESC   : get bdy class for specified bdy type
// ARGS   : sptree  : bdy spec subtree
//          grid    : grid object
//          id      : canonical bdy id
//          bdytype : bdy type triggering construction of class
//          istate  : vector of state ids
//          ibdy    : indirection indices into computational volume, 
//                    representing the bdy nodes this method applies to
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBasePtr
GUpdateBdyFactory<TypePack>::get_bdy_class(const PropertyTree& sptree, Grid &grid, const GINT id, const GBdyType bdytype, GTVector<GINT> &istate, GTVector<GSIZET> &ibdy)
{
  GINT               nstate;
  GString            sblock ;
  UpdateBdyBasePtr   base_ptr;
  std::vector<Ftype> fvec;
  

  sstr = sptree.getValue<GString>("update_method");
  nstate = istate.size();
  

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    assert(sptree.isArray<Ftype>("value") && "value array missing");
    traits.value .resize(nstate); traits.value = fvec;
    traits.bdyid = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }

    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = obs_impl;
    return base_ptr;
  }
  else if ( GBDY_INFLOW == bdytype ) {
    using UpdateImpl = GInflowBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    assert( sptree.isValue<GBOOL>("use_init") && "use_init boolean missing") {
    traits.use_init = sptree.getValue<GBOOL>("use_init");
    if ( !traits.use_init ) {
      assert( sptree.isValue<GString>("inflow_method") 
           && "inflow_method required if use_init==FALSE" ) 
      sblock = sptree.getValue<GString>("inflow_method");
      traits.callback = get_inflow_callback(sblock, id);
    }
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

  }
  else if ( GBDY_OUTFLOW == bdytype ) {
    assert(FALSE); // not available yet
    using UpdateImpl = GOutflowBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate.resize(nstate); traits.istate = istate;
    traits.bdyid  = id;
    traits.ibdyvol = ibdy;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<TypesPack>
    UpdateImpl::Traits traits;

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
    
    GINT            ndim;
    GTVector<Ftype> xmin, xmax;
    geoflow::coord_dims(ptree, xmin, xmax, ndim);// get coord min/max from ptree

    assert(sptree.isValue<GINT>("idir") && "idir value missing");
    traits.idir = sptree.getValue<GINT>("idir");
    assert(traits.idir >= 1 && traits.idir <= xmax.size()); // validate idir

    if ( traits.idir < 0 ) traits.ro = xmin[abs(traits.idir)-1];
    else                   traits.ro = xmax[abs(traits.idir)-1];

    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

  }
  else {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return base_ptr;

} // end, init method get_bdy_class


