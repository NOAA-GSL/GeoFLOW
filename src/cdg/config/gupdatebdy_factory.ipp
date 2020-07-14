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
//          ibdy    : bdy indirection indices into computational volume
//          tbdy    : bdy types for each ibdy index
//                    
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBasePtr
GUpdateBdyFactory<TypePack>::build(const PropertyTree& sptree, GString &supdate, Grid &grid, const GINT id, BdyIndices &ibdy, BdyTypes &tbdy)
{
  GBOOL              bret = FALSE;
  GString            sclass;
  PropertyTree       sptree; // update block tree
  UpdateBdyBasePtr   base_ptr;
  UpdateBdyBaseList  list;


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
  sclass = sptree.getValue<GString>("bdy_class");


  // If bdy_class is uniform, use standard config method:
  if      ( "uniform"    == sclass) {
    list = handle_uniform(sptree, supdate, grid, id, ibdy, tbdy);   
  }
  else if ( "mixed"      == sclass) {
    list = handle_mixed(sptree, supdate, grid, id, ibdy, tbdy);   
  }
  else {
    assert(FALSE && "Invalid bdy_class option");
  }

} // end, init method build


//**********************************************************************************
//**********************************************************************************
// METHOD : handle_uniform
// DESC   : Build bdy update object for uniform bdys
// ARGS   : ptree   : main property tree
//          supdate : string naming bdy update prop tree block 
//          grid    : grid object
//          id      : canonical bdy id
//          ibdy    : bdy indirection indices into computational volume
//          tbdy    : bdy types for each ibdy index
//                    
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBaseList
GUpdateBdyFactory<TypePack>::handle_uniform(const PropertyTree& sptree, GString &supdate, Grid &grid, const GINT id, BdyIndices &ibdy, BdyTypes &tbdy)
{
  GString            sblock, stype, sclass;
  GBdyType           bdytype;
  PropertyTree       sptree; // update block tree
  UpdateBdyBasePtr   base_ptr;
  UpdateBdyBaseList  list;
  
  
 
  if ( !sptree.isValue<GString>("base_type") ) {
    cout << "GUpdateBdyFactory<TypePack>::handle_uniform: base_type not specified in " << supdate << " block" << endl;
    assert(FALSE);
  };
  stype   = sptree.getValue<GString>("base_type");
  bdytype = geoflow::str2bdytype(stype);

  bret = gspecbdy::impl_uniform    (sptree, grid, id, ibdy, tbdy);
  assert(bret && "gspecbdy::impl_uniform failed"");

  // If canonical bdy is uniform, then bdytype 
  // specifies which update method to call. ibdy, and tbdy
  // are filled outside of this method:
  base_ptr = get_bdy_class (sptree, supdate, grid, id, bdytype); 
  list.push_back(base_ptr);

  return list;

} // end, init method handle_uniform


//**********************************************************************************
//**********************************************************************************
// METHOD : handle_mixed
// DESC   : Build bdy update object for mixed bdys
// ARGS   : ptree   : main property tree
//          supdate : string naming bdy update prop tree block 
//          grid    : grid object
//          id      : canonical bdy id
//          ibdy    : bdy indirection indices into computational volume
//          tbdy    : bdy types for each ibdy index
//                    
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBaseList
GUpdateBdyFactory<TypePack>::handle_mixed(const PropertyTree& sptree, GString &supdate, Grid &grid, const GINT id, BdyIndices &ibdy, BdyTypes &tbdy)
{
  GBOOL              bret;
  GString            sspec;
  GBdyType           bdytype;
  PropertyTree       sptree; // update block tree
  UpdateBdyBasePtr   base_ptr;
  UpdateBdyBaseList  list;
  std::vector<GINT>  ivec;
  std::vector<Ftype> fvec;
  

  sspec = sptree.getValue<GString>("config_method","");
 
  if ( "none"     == sspec
    || ""         == sspec ) {
    cout << "GUpdateBdyFactory<TypePack>::handle_mixed: config_method not specified in " << supdate << " block" << endl;
    assert(FALSE);
  };

  // First, set mixed bdy conditions on bdy:
  if (   "my_mixed_bdy" == sspec ) {
    bret = gspecbdy::impl_my_mixed_bdy(sptree, grid, id, ibdy, tbdy);
    assert(bret && "gspecbdy::impl_my_mixed_bdyfailed"");
  }
  else {
    cout << "GUpdateBdyFactory<TypePack>::handle_mixed: Invalid config_method, " << sspec << " in " << supdate << " block" << endl;
    assert(FALSE);
  {

  // Search tbdy types returned from sspec call, and 
  // create a base_ptr for for each unique type found:
  GTVector<GBdyType> iunique;
  tbdy.unique(0, tbdy.size()-1, iunique);
  for ( auto j=0; j<iunique.size(); j++ ) {
    base_ptr = get_bdy_class (sptree, supdate, grid, id, tbdy[iunique[j]] ); 
    list.push_back(base_ptr);
  }

  return list;

} // end, init method handle_mixed


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


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_class
// DESC   : get bdy class for specified bdy type
// ARGS   : ptree   : main property tree
//          supdate : string naming bdy update prop tree block 
//          grid    : grid object
//          id      : canonical bdy id
//          bdytype : bdy type triggering construction of class
//                    
//                   
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
typename UpdateBdyFactory<TypePack>::UpdateBasePtr
GUpdateBdyFactory<TypePack>::get_bdy_class(const PropertyTree& sptree, GString &supdate, Grid &grid, const GINT id, const GBdyType bdytype)
{
  GString            sblock, stype, sclass;
  GBdyType           bdytype;
  PropertyTree       sptree; // update block tree
  UpdateBdyBasePtr   base_ptr;
  std::vector<GINT>  ivec;
  std::vector<Ftype> fvec;
  

  sstr = sptree.getValue<GString>("update_method");
  assert(sptree.isArray<GINT>("istate") && "istate vector missing");
  ivec = sptree.getArray<GINT>("istate");

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    assert(sptree.isArray<Ftype>("value") && "value array missing");
    traits.value = fvec;
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

    traits.istate = ivec;
    assert(sptree.isArray<Ftype>("value") && "value array missing");
    traits.value = fvec;

    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    assert( sptree.isValue<GBOOL>("use_init") && "use_init boolean missing") {
    traits.use_init = sptree.getValue<GBOOL>("use_init");
    if ( !traits.use_init ) {
      assert( sptree.isValue<GString>("inflow_method") 
           && "inflow_method required if use_init==FALSE" ) 
      sblock = sptree.getValue<GString>("inflow_method");
      traits.callback = get_inflow_callback(sblock);
    }
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
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

    traits.istate = ivec;
    if ( sptree.isValue<GBOOL>("compute_once") ) {
      traits.compute_once = sptree.getValue<GBOOL>("compute_once");
    }
    
    // Allocate observer Implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<TypesPack>
    UpdateImpl::Traits traits;

    traits.istate = ivec;
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


