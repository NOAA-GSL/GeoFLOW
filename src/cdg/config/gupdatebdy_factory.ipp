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
//          sbdy    : bdy condition name
//          grid    : grid object
//          bcblock : stBdyBlock structure
//          ibdy    : bdy indirection indices into computational volume
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::UpdateBdyBasePtr
GUpdateBdyFactory<Types>::build(const PropertyTree& ptree, const GString &sbdy, Grid &grid, stBdyBlock &bcblock, GTVector<GSIZET> &ibdy, GSIZET igbdy_start)
{
  GBOOL              bret = FALSE;
  UpdateBdyBasePtr   base_ptr;


  if ( "none"         == sbdy
    || ""             == sbdy ) {
    using UpdateImpl = NullUpdateBdy<Types>;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }

  if ( !ptree.isPropertyTree(sbdy) ) {
    cout << "GUpdateBdyFactory<Types>::build: PropertyTree " << sbdy << " not found" << endl;
    assert(FALSE);
  }

  base_ptr = GUpdateBdyFactory<Types>::get_bdy_class(ptree, grid, bcblock, ibdy, igbdy_start);

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
//          grid    : grid object
//          bcblock : stBdyBlock structure
//          ibdy    : indirection indices into computational volume, 
//                    representing the bdy nodes this method applies to
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::UpdateBdyBasePtr
GUpdateBdyFactory<Types>::get_bdy_class(const PropertyTree &ptree, Grid &grid, stBdyBlock &bcblock, GTVector<GSIZET> &ibdy, GSIZET igbdy_start)
{
  GINT               nstate;
  GBdyType           bdytype = bcblock.tbdy;
  GLONG              iloc;
  UpdateBdyBasePtr   base_ptr;
  

  nstate = bcblock.istate.size();
  

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<Types>;
    typename GDirichletBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    traits.value .resize(bcblock.value.size())  ; 
    traits.value       = bcblock.value;;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_INFLOW == bdytype ) {
    using UpdateImpl = GInflowBdy<Types>;
    typename GInflowBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    if ( !bcblock.use_init ) {
      traits.callback = GUpdateBdyFactory<Types>::get_inflow_callback(bcblock.smethod, bcblock.bdyid);
    }
    traits.ptree = ptree;
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
#if 0
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<Types>;
    typename GNoSlipBdy<Types>::Traits traits;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
#endif
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<Types>;
    typename G0FluxBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;

    traits.ibdyloc.resize(traits.ibdyvol.size());
    // Find index of ibdyvol in global bdy vector:
    for ( auto j=0; j<traits.ibdyloc.size(); j++ ) {
      iloc = j + igbdy_start;
      assert(iloc >= 0);
      traits.ibdyloc[j] = iloc;
    }
    
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_OUTFLOW == bdytype ) {
    assert(FALSE); // not available yet
    using UpdateImpl = GOutflowBdy<Types>;
    typename GOutflowBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<Types>;
    typename GSpongeBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    traits.idir     = bcblock.idir;
    traits.rs       = bcblock.xstart;
    traits.ro       = bcblock.xmax;

    traits.farfield.resize(bcblock.farfield.size());
    traits.exponent.resize(bcblock.falloff .size());
    traits.sigma   .resize(bcblock.diffusion.size());

    traits.farfield = bcblock.farfield;
    traits.exponent = bcblock.falloff;
    traits.sigma    = bcblock.diffusion;
    
    // Allocate update implementation
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


