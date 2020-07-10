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
//          supdate: string naming bdy update prop tree
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBasePtr
GUpdateBdyFactory<TypePack>::build(const PropertyTree& ptree, GString &supdate, GBdyType bdytype)
{
  GBOOL              bret = FALSE;
  GString            sblock;
  PropertyTree       bptree;
  UpdateBdyBasePtr   base_ptr;
  std::vector<GINT>  ivec;
  std::vector<Ftype> fvec;

  bptree     = ptree.getPropertyTree(supdate);

  if      ( "none"         == sblock
    || ""             == sblock ) {
    using UpdateImpl = NullUpdateBdy<TypesPack>

    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = obs_impl;
    return base_ptr;
  }

  uptree    = ptree.getPropertyTree(sblock);
  supdate   = uptree.getValue<GString>("update_method");

  assert(bptree.isArray<GINT>("istate") && "istate vector missing");
  ivec = bptree.getArray<GINT>("istate");

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    assert(bptree.isArray<Ftype>("value") && "value array missing");
    traits.value = fvec;
    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
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

    traits.istate = ivec;
    assert(bptree.isArray<Ftype>("value") && "value array missing");
    traits.value = fvec;

    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
    }
    assert( bptree.isValue<GBOOL>("use_init") && "use_init boolean missing") {
    traits.use_init = bptree.getValue<GBOOL>("use_init");
    if ( !traits.use_init ) {
      assert( bptree.isValue<GString>("inflow_method") 
           && "inflow_method required if use_init==FALSE" ) 
      sblock = bptree.getValue<GString>("inflow_method");
      traits.callback = get_inflow_callback(sblock);
    }
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

  }
  else if ( GBDY_OUTFLOW == bdytype ) {
    assert(FALSE); // not available yet
    using UpdateImpl = GOutflowBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    if ( bptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = bptree.getValue<GBOOL>("compute_once");
    }

    assert(bptree.isArray<Ftype>("farfield") && "farfield array missing");
    traits.idir = bptree.getValue<Ftype>("farfield");
    assert(bptree.isArray<Ftype>("exponent") && "exponent array missing");
    traits.idir = bptree.getValue<Ftype>("exponent");
    assert(bptree.isArray<Ftype>("sigma0") && "sigma0 array missing");
    traits.idir = bptree.getValue<Ftype>("sigma0");
    assert(bptree.isArray<Ftype>("xstart") && "xstart array missing");
    traits.idir = bptree.getValue<Ftype>("xstart");
    
    GINT            ndim;
    GTVector<Ftype> xmin, xmax;
    geoflow::coord_dims(ptree, xmin, xmax, ndim);// get coord min/max from ptree

    assert(bptree.isValue<GINT>("idir") && "idir value missing");
    traits.idir = bptree.getValue<GINT>("idir");
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

} // end, init method build


//**********************************************************************************
//**********************************************************************************
// METHOD : get_inflow_callback
// DESC   : Gets CallbackPtr corresponding to sname arg for inflow conditions.
//          the function are gottne from the gns_inflow_user.* namespace
//          collection of methods, which the user may modify.
// ARGS   : sname : inflow function name
// RETURNS: CallbackPtr for callback function
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::CallbackPtr
GUpdateBdyFactory<TypePack>::get_inflow_callback(const GString& sname)
{
  GBOOL              bret = FALSE;
  CallbackPtr        callback;


  if      ( "myinflow"     == sname ) {
    callback = 
      [GInflowBdyMethods](Grid      &grid,
                          StateInfo &stinfo,
                          Time      &time,
                          State     &utmp,
                          State     &u,
                          State     &ub){myinflow(grid, stinfo, time, utmp, u, ub);}; 
  else {
    assert(FALSE && "Specified inflow bdy update method unknown");
  }

  return callback;

} // end, init method get_inflow_callback


