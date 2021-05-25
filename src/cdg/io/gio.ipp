//==================================================================================
// Module       : gio.hpp
// Date         : 1/20/20 (DLR)
// Description  : GIO object encapsulating methods for POSIX and collective IO
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : IOBase.
//==================================================================================

#include "tbox/tracer.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : grid  : grid object
//          traits: Traits sturcture
//          comm  : communicator
//**********************************************************************************
template<typename Types>
GIO<Types>::GIO(Grid &grid,  Traits &traits, GC_COMM comm):
IOBase<Types>(grid, traits),
bInit_                     (FALSE),
myrank_   (GComm::WorldRank(comm)),
comm_                       (comm),
cfname_                  (NULLPTR),
nfname_                        (0)
{ 
  GEOFLOW_TRACE();
#if !defined(GEOFLOW_USE_MPI)
  assert(this->traits_->io_type == IOBase<Types>::GIO_COLL && "Collective IO only allowed if MPI is used");
#endif
  init();

} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   :
// ARGS   : none.
//**********************************************************************************
template<typename Types>
GIO<Types>::~GIO()
{ 
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_MPI)
  MPI_Type_free(&mpi_state_type_);
#endif
} // end of destructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : update_type
// DESCRIPTION: Update datatype to handle writes/reads, for each IO call,
//              using StateInfo data
// ARGUMENTS  : info: StateInfo structure
// RETURNS    : none.
//**********************************************************************************
template<typename Types>
void GIO<Types>::update_type(StateInfo &info)
{
  GEOFLOW_TRACE();

  // Update formats based on info. This must be called before 
  // test on traits_->io_type:
  scformat_.str(""); spformat_.clear();
  spformat_.str(""); spformat_.clear();
  if ( info.sttype ==0 ) { // is a physical state
    spformat_   << "%s/%s.%0" << this->traits_.wtime << "d.%0" << this->traits_.wtask << "d.out";
    scformat_   << "%s/%s.%0" << this->traits_.wtime << "d.out";
  }
  else {                   // is a grid file
    spformat_   << "%s/%s.%0" << this->traits_.wtask << "d.out";
    scformat_   << "%s/%s.%0" << this->traits_.wtime << "d.out";
  }

} // end, update_type


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: Initialize object
// ARGUMENTS  : info: StateInfo structure
// RETURNS    : none.
//**********************************************************************************
template<typename Types>
void GIO<Types>::init()
{
  GEOFLOW_TRACE();

  if ( this->traits_.io_type != IOBase<Types>::GIO_COLL ) {
    bInit_ = TRUE;
    return; // nothing more to do
  }


  GINT             iret, state_disp, state_extent, sztot;
  GSIZET           ndof  =this->grid_->ndof();
  GINT             nelems=this->grid_->nelems();
  GTVector<GSIZET> extent;


  // Get extents of a single state component on each task:
  extent.resize(GComm::WorldSize(comm_));
  iret = GComm::Allgather(&ndof, 1, T2GCDatatype<GSIZET>(), extent.data(), 1, T2GCDatatype<GSIZET>(), comm_);

  state_disp   = extent.sum(0,myrank_-1); // cumulative sum
  state_extent = extent[myrank_]; // count
  nbgdof_ = extent.sum()*sizeof(Ftype); // no. bytes of single state comp

  // Find displacements for each set of elem ids:
  GTVector<GINT>   iddisp;;
  GTVector<GINT>   elemcount;;
  elemcount.resize(GComm::WorldSize(comm_));
  iddisp   .resize(GComm::WorldSize(comm_));
  iret = GComm::Allgather(&nelems, 1, T2GCDatatype<GINT>(), elemcount.data(), 1, T2GCDatatype<GINT>(), comm_);
  for ( auto i=0; i<iddisp.size(); i++ ) {
    iddisp[i] = elemcount.sum(0,i-1);
  }
  

#if defined(GEOFLOW_USE_MPI)
//MPI_Type_free(&mpi_state_type_);


  sztot = extent.sum();
  iret = MPI_Type_create_subarray(1, &sztot, &state_extent, &state_disp, MPI_ORDER_C, T2GCDatatype<Ftype>(), &mpi_state_type_);
  assert(iret == MPI_SUCCESS);
  iret = MPI_Type_commit(&mpi_state_type_);
  assert(iret == MPI_SUCCESS);

  // Use id displacements to get element keys/ids 
  // so that task 0 can write them in header:
  GTVector<GKEY> *keys = &this->grid_->elemids();
  assert( keys != NULLPTR && keys->size() == this->grid_->nelems() ); // verify there are the right #
  hkeys_.resize(this->grid_->ngelems());
  iret = GComm::Gatherv(keys->data(), keys->size(), T2GCDatatype<GKEY>(), 
                        hkeys_.data(), elemcount.data(), iddisp.data(), T2GCDatatype<GKEY>(),            
                        0, comm_);
#endif
  bInit_ = TRUE;

} // end of method init



//**********************************************************************************
//**********************************************************************************
// METHOD : write_state_impl
// DESC   : Do GIO POSIX or collective output of state. 
//          If there aren't enough svar specified, or if filename isn't specified,
//          when needed, then it's an error.
//          Note: if info.sttype > 0, then we assume we are printing a 
//          grid, and the filename is created without the time index tag.
// ARGS   : filepref: used if traits.multivar > 0 to specify file name for
//                    all state variables. This works only if we GTypes is GIO_COLL.
//                    If traits.multivar ==0, individual filename refixes are provided 
//                    in info.svar
//          info    : StateInfo structure
//          u       : state
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GIO<Types>::write_state_impl(std::string filepref, StateInfo &info, const State &u)
{
  GEOFLOW_TRACE();
  GString        serr = "write_state_impl: ";
  GSIZET         nb, nc, nd;
  GTVector<GTVector<Ftype>>
                *xnodes = &(this->grid_->xNodes());
  State         ostate(1);
  typename Grid::GElemList     *elems  = &(this->grid_->elems());

  assert(bInit_ && "Object uninitialized");

  update_type(info);

  // Required number of coord vectors:
  nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  if ( info.sttype !=0 ) { // is a grid file 
    assert(nc ==  xnodes->size()
        && "Incompatible grid or coord  dimensions");
  }
    
  this->traits_.dim  = GDIM;

  if ( this->traits_.io_type == IOBase<Types>::GIO_COLL ) {
    info.nelems = this->grid_->ngelems(); // total no. elems among all tasks
  }
  else {
    info.nelems = this->grid_->nelems(); // local nelems for this task
  }
  info.gtype    = this->grid_->gtype();


  // Build format strings:
  resize(this->traits_.wfile);

  // Set porder vector depending on version:
  info.porder.resize(this->traits_.ivers == 0 ? 1 : info.nelems, GDIM);
  for ( auto i=0; i<info.porder.size(1); i++ )  { // for each element
    for ( auto j=0; j<info.porder.size(2); j++ ) info.porder(i,j) = (*elems)[i]->order(j);
  }

  if ( !this->traits_.multivar ) { // one state comp per file
    assert(info.svars.size() >= u.size());
    // Cycle over all fields, and write:
    for ( auto j=0; j<u.size(); j++ ) {
      assert(u[j]->size() > 0 && "Invalid state component");
      svarname_.str(""); svarname_.clear();
      assert(info.svars[j].length() > 0);
      svarname_ << info.svars[j];
      if ( this->traits_.io_type == IOBase<Types>::GIO_POSIX ) {
        sprintf(cfname_, spformat_.str().c_str(), info.odir.c_str(),
                svarname_.str().c_str(), info.index, myrank_);
        fname_.assign(cfname_);
        nb = write_posix(fname_, info, *u[j]); // only writes one variable per file
      }
      else { // GIO_COLL, collective write
        sprintf(cfname_, scformat_.str().c_str(), info.odir.c_str(),
                svarname_.str().c_str(), info.index);
        fname_.assign(cfname_);
        ostate[0] = u[j];
        nb = write_coll(fname_, info, ostate);
      }
      nd = sz_header(info,this->traits_) + u[j]->size()*sizeof(Ftype);
      assert(nb == nd && "Incorrect number of bytes written");
    }
  }
  else {                      // multiple components per file
    assert(this->traits_.io_type == IOBase<Types>::GIO_COLL && "Invalid io_type");
    svarname_.str(""); svarname_.clear();
    assert(filepref.length() > 0);
    svarname_ << filepref;
    sprintf(cfname_, scformat_.str().c_str(), info.odir.c_str(),
            svarname_.str().c_str(), info.index);
    nb = write_coll(fname_, info, u);
    nd = sz_header(info,this->traits_) + u.size()*u[0]->size()*sizeof(Ftype);
    assert(nb == nd && "Incorrect number of bytes written");
  }

} // end, write_state_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : read_state_impl
// DESC   : Read files for entire state, based on StateInfo. 
//          
// ARGS   : filepref: filename prefix, if traits.multivar > 0; else
//                    info.svars are used to form filenames
//          info    : StateInfo structure
//          u       : state object. Method will attempt to fill all components of u.
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GIO<Types>::read_state_impl(std::string filepref, StateInfo &info, State  &u, bool bstate)
{
  GEOFLOW_TRACE();
  GString              serr ="read_state_impl: ";
  GINT                 ivers, nc, nr, nt;
  GElemType            igtype;
  GString              sgtype;
  State                istate(1);
  Traits               ttraits;

  assert(bInit_ && "Object uninitialized");
  assert(info.icomptype.size() >= u.size() && "Stateinfo structure invalid");

  update_type(info);

  nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM; 

  if ( bstate && info.sttype != 0 ) { // Check if reading grid components:
    assert(u.size() == nc && "Incorrect number of grid components");
  }

  resize(this->traits_.wfile);
  
  if ( !this->traits_.multivar ) { // one state comp per file

    assert(info.svars.size() >= 1);
    nt = bstate ? u.size() : 1;
    for ( GSIZET j=0; j<nt; j++ ) { // Retrieve all state/grid components
      if ( info.icomptype[j] == GSC_PRESCRIBED ) continue;
      svarname_.str(""); svarname_.clear();
      assert(info.svars[j].length() > 0);
      svarname_ << info.svars[j];
      if ( this->traits_.io_type == IOBase<Types>::GIO_POSIX ) { // POSIX
        sprintf(cfname_, spformat_.str().c_str(), info.idir.c_str(),
                svarname_.str().c_str(), info.index, myrank_);
        fname_.assign(cfname_);
        nr = read_posix(fname_, info, *u[j], bstate);
      }
      else {  // collective
        sprintf(cfname_, scformat_.str().c_str(), info.idir.c_str(),
                svarname_.str().c_str(), info.index);
        fname_.assign(cfname_);
        istate[0] = u[j];
        nr = read_coll (fname_, info, istate, bstate);
      }
    }
  }
  else {                      // multiple state components in file

    assert(this->traits_.io_type == IOBase<Types>::GIO_COLL && "Invalid io_type");
    svarname_.str(""); svarname_.clear();
    assert(filepref.length() > 0);
    svarname_ << filepref;
    sprintf(cfname_, scformat_.str().c_str(), info.idir.c_str(),
            svarname_.str().c_str(), info.index);
    read_coll(fname_, info, u, bstate);

  }

} // end, read_state_impl method


//**********************************************************************************
//**********************************************************************************
// METHOD : read_state_info_impl
// DESC   : Read files for StateInfo.
//          
// ARGS   : filepref: filename prefix 
//          info    : StateInfo structure, returned
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GIO<Types>::read_state_info_impl(std::string filename, StateInfo &info)
{
  GEOFLOW_TRACE();
  GString              serr ="read_state_info_impl: ";
  GSIZET               nh;
  Traits               ttraits;


  nh  = read_header(filename, info, ttraits);

  assert( nh == sz_header(info,this->traits_) );

} // end, read_state_info_impl method


//**********************************************************************************
//**********************************************************************************
// METHOD : write_posix
// DESC   : Do simple (lowes-level) GIO POSIX output of specified field
// ARGS   : 
//          filename: filename, fully resolved
//          info    : StateInfo structure
//          u       : field to output
// RETURNS: number bytes written
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::write_posix(GString filename, StateInfo &info, const GTVector<Ftype> &u) 
{
  GEOFLOW_TRACE();

    GString  serr ="write_posix: ";
    FILE     *fp;
    GSIZET    i, j, nb, nd;

    // Write header; remember that file is closed on exit:
    nb = write_header_posix(filename, info, this->traits_);
    assert(nb == sz_header(info,this->traits_));
   
    fp = fopen(filename.c_str(),"wb");
    assert(fp != NULL && "Error opening file");

    if ( this->traits_.ivers > 0 && info.porder.size(1) <= info.nelems ) {
      cout << serr << " porder of insufficient size for version: " << this->traits_.ivers
                   << ". Error writing file " << filename << endl;
      exit(1);
    }


    // Write field data:
    nd = fwrite(u.data(), sizeof(Ftype), u.size(), fp);
    fclose(fp);

    nb +=  nd*sizeof(Ftype);

    return nb;

} // end, write_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : read_posix
// DESC   : Do simple GIO POSIX read
// ARGS   : 
//          filename: filename, fully resolved
//          info    : StateInfo structure, filled here
//          u       : field to output; resized if required
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: no. bytes read
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::read_posix(GString filename, StateInfo &info, GTVector<Ftype> &u, bool bstate)
{
  GEOFLOW_TRACE();

    GString serr = "read_posix: ";
    FILE     *fp;
    GSIZET    nb, nd, nh, nt;
    Traits    ttraits;
    
    nh  = read_header(filename, info, ttraits);

    assert(ttraits.ivers == this->traits_.ivers
                                         && "Incompatible file version number");
    assert(ttraits.dim   == GDIM         && "File dimension incompatible with GDIM");
    assert(info   .gtype == this->grid_->gtype() 
                                           && "File grid type incompatible with grid");
    if ( !bstate ) return nh;

    fp = fopen(filename.c_str(),"rb");
    assert(fp != NULL && "Error opening file");

    // Seek to first byte after header:
    fseek(fp, nh, SEEK_SET); 

    // Compute field data size from header data:
    nd = 0;
    if ( ttraits.ivers == 0 ) { // expansion order is constant
      nt = 1; 
      for ( GSIZET j=0; j<GDIM; j++ ) nt *= (info.porder(0,j) + 1);
      nd += nt * info.nelems;
    }
    else {                     // expansion order varies
      for ( GSIZET i=0; i<info.porder.size(1); i++ ) {
        nt = 1; 
        for ( GSIZET j=0; j<GDIM; j++ ) nt *= (info.porder(i,j) + 1);
        nd += nt;
      }
    }

    u.resize(nd);
    nb = fread(u.data(), sizeof(Ftype), nd, fp);
    
    fclose(fp);

    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << filename << endl;
      exit(1);
    }
    nh += nb*sizeof(Ftype);

   return nh;

} // end, read_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : write_header_posix
// DESC   : Write GIO POSIX file header. Note file pointer position 
//          may be changed on exit.
// ARGS   : 
//          filename : file name
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes written
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::write_header_posix(GString filename, StateInfo &info, Traits &traits)
{
  GEOFLOW_TRACE();
    GString serr ="write_header_posix: ";
    GINT   imulti = static_cast<GINT>(traits.multivar);
    GSIZET nb, nd, nh, numr;
    FILE  *fp;
  
    fp = fopen(filename.c_str(),"wb");
    assert(fp != NULLPTR && "error opening file");

    nb = 0;
  
    fseek(fp, 0, SEEK_SET);
    // Write header: dim, numelems, poly_order:
    nh=fwrite(&traits.ivers     , sizeof(GINT)  ,    1, fp);   // GIO version number
      nb += nh*sizeof(GINT);
    nh=fwrite(&traits.dim       , sizeof(GINT)  ,    1, fp);   // dimension
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.nelems      , sizeof(GSIZET),    1, fp);   // num elements
      nb += nh*sizeof(GSIZET);
    nd = info.porder.size(1) * info.porder.size(2);
    nh=fwrite(info.porder.data().data()
                               , sizeof(GINT)     ,   nd, fp); // exp order
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.gtype         , sizeof(GINT)  ,    1, fp); // grid type
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.cycle         , sizeof(GSIZET),    1, fp); // time cycle stamp
      nb += nh*sizeof(GSIZET);
    nh=fwrite(&info.time          , sizeof(Ftype) ,    1, fp); // time stamp
      nb += nh*sizeof(Ftype);
    nh=fwrite(&imulti             , sizeof(GINT)  ,    1, fp); // time stamp
      nb += nh*sizeof(GINT);

    // Add list of element ids (keys) always at the end of header:
    GTVector<GKEY> *keys = &this->grid_->elemids();
    assert( keys != NULLPTR && keys->size() == info.nelems ); // verify there are the right number
    nh=fwrite(keys->data()        , sizeof(GKEY)  ,    info.nelems, fp);// keys/ids
      nb += nh*sizeof(GKEY);

    fclose(fp);
  
    return nb;

} // end, write_header_posix


#if defined(GEOFLOW_USE_MPI)
//**********************************************************************************
//**********************************************************************************
// METHOD : write_header_coll
// DESC   : Write GIO MPI file header. Note file pointer position 
//          may be changed on exit.
// ARGS   : 
//          filename : filename
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes written
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::write_header_coll(GString filename, StateInfo &info, Traits &traits)
{
  GEOFLOW_TRACE();
  GString serr ="write_header: ";
  GINT       nh, imulti, iret;
  GSIZET     nb, numr;
  MPI_File   fp;
  MPI_Status status;

  
  iret = MPI_File_open(comm_, filename.c_str(), MPI::MODE_CREATE|MPI::MODE_WRONLY, MPI::INFO_NULL, &fp);
  assert(iret == MPI_SUCCESS && "MPI_File_open failure");

  nb = sz_header(info,traits);
  if ( myrank_ == 0 ) {
    imulti = static_cast<GINT>(traits.multivar);
    nb = 0;
    MPI_File_seek(fp, 0, MPI_SEEK_SET); // set to 0-displacement
    
    MPI_File_write(fp, &traits.ivers   , 1   , T2GCDatatype  <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &traits.dim     , 1   , T2GCDatatype  <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &info.nelems    , 1   , T2GCDatatype<GSIZET>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    numr = traits.ivers == 0 ? 1 : info.nelems;
    info.porder.resize(numr,traits.dim);
    numr = info.porder.size(1)*info.porder.size(2);
    MPI_File_write(fp, info.porder.data().data()
                                               , numr, T2GCDatatype   <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    MPI_File_write(fp, &info.gtype       , 1   , T2GCDatatype  <GINT>(), &status);
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status);
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    MPI_File_write(fp, &info.time        , 1   , T2GCDatatype<Ftype>(),  &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    MPI_File_write(fp, &imulti           , 1   , T2GCDatatype  <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;

    // Element ids/keys:
    GTVector<GKEY> *keys = &this->grid_->elemids();
    assert( hkeys_.size() == info.nelems ); // verify there are the right number
    MPI_File_write(fp, hkeys_.data()     , info.nelems, T2GCDatatype  <GKEY>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
  }


  MPI_File_close(&fp); 
     
  GComm::Synch(comm_);

  return nb;

} // end, write_header_coll
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : read_header
// DESC   : Read GIO file header
// ARGS   : 
//          filename : file name (fully resolved)
//          info     : StateInfo structure, filled with what header provides
//          traits   : this object's traits
// RETURNS: no. header bytes read
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::read_header(GString filename, StateInfo &info, Traits &traits)
{
  GEOFLOW_TRACE();

    GString serr ="read_header: ";
    GINT imulti ;
    GSIZET nb, nd, nh, numr;
  
    nb = 0;
//  if ( traits.io_type == IOBase<Types>::GIO_POSIX ) {
      // Read field data:
      FILE *fp;
      fp = fopen(filename.c_str(),"rb");
      assert(fp != NULLPTR && "gio.cpp: error opening file");
    
      // Read header: 
      nh = fread(&traits.ivers     , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&traits.dim       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.nelems      , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = fread(info.porder.data().data()
                                   , sizeof(GINT)  , numr, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.gtype       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.cycle       , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
      nh = fread(&info.time        , sizeof(Ftype) ,    1, fp); nb += nh*sizeof(Ftype);
      nh = fread(&imulti           , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      traits.multivar = static_cast<GBOOL>(imulti);

      info.elemids.resize(info.nelems);
      nh = fread( info.elemids.data(), sizeof(GKEY),    info.nelems, fp); nb += nh*sizeof(GKEY);
    
      fclose(fp);
  
#if 0
    }
  
    else {
#if defined(GEOFLOW_USE_MPI)
      MPI_File        fh;
      MPI_Status      status;
      MPI_File_open(comm_, filename.c_str(), MPI::MODE_RDWR, MPI::INFO_NULL, &fh);
//    MPI_File_seek(fh, 0, MPI_SEEK_SET); // set to 0-displacement
      
      // Read header: 
      nh = MPI_File_read(fh, &traits.ivers     , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &traits.dim       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.nelems      , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = MPI_File_read(fh, info.porder.data().data()
                                               , numr, T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.gtype       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      nh = MPI_File_read(fh, &info.time        , 1   , T2GCDatatype<Ftype>(), &status); nb += nh*sizeof(Ftype);

      info.elemids.resize(info.nelems);
      nh = MPI_File_read(fh, info.elemids.data() , info.nelems  , T2GCDatatype<GKEY>(), &status); nb += nh*sizeof(GKEY);

      MPI_File_close(&fh); 
#else
# error "MPI not defined; cannot read GIO_COLL"
#endif
    } 
#endif

    // Check number read vs expected value:
    nd = sz_header(info, traits);
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << filename << endl;
      exit(1);
    }

    return nb;

} // end, read_header


//**********************************************************************************
//**********************************************************************************
// METHOD : sz_header
// DESC   : Get header size in bytes
// ARGS   : 
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes 
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::sz_header(const StateInfo &info, const Traits &traits)
{
  GEOFLOW_TRACE();

    GString serr ="sz_header: ";
    GSIZET nd, numr;
  
    // Get no. bytes in header (should agree with read_posx, write_posix):
    numr = traits.ivers == 0 ? 1 : info.nelems;
    numr *= traits.dim;
    nd = (numr+4)*sizeof(GINT) + 2*sizeof(GSIZET) + sizeof(Ftype);

    if ( this->traits_.io_type == IOBase<Types>::GIO_COLL ) {
      nd += this->grid_->ngelems() * sizeof(GKEY);
    }
    else if ( this->traits_.io_type == IOBase<Types>::GIO_POSIX ) {
      nd += this->grid_->nelems () * sizeof(GKEY);
    }

    return nd;

} // end, sz_header


//**********************************************************************************
//**********************************************************************************
// METHOD : resize
// DESC   : Resize global data
// ARGS   : n  : new number bytes
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GIO<Types>::resize(GINT n)
{
  GEOFLOW_TRACE();

  if ( n > nfname_ ) {
    if ( cfname_ != NULLPTR ) delete [] cfname_;
    cfname_ = new char [n];
    fname_.resize(n);
    nfname_ = n;
  }

} // end, resize


//**********************************************************************************
//**********************************************************************************
// METHOD : write_coll
// DESC   : Collective write of state components
// ARGS   : filename: filename
//          info    : StateInfo structure
//          u       : state
// RETURNS: number bytes written
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::write_coll(GString filename, StateInfo &info, const State &u)
{
  GEOFLOW_TRACE();
#if !defined(GEOFLOW_USE_MPI)
  #error "Illegal entry into GIO<Types>::write_coll: MPI not defined"
#endif

#if defined(GEOFLOW_USE_MPI)

    GString        serr = "write_coll: ";
    GINT           iret, nbheader, nc, nh;
    GSIZET         nb;
    GSIZET         ntot;

    MPI_Offset     disp;
    MPI_File       fh;
    MPI_Status     status;

    // Required number of coord vectors:
    nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

    nbheader = sz_header(info, this->traits_);

    // Write header; remember that file is closed on exit:
    nh = write_header_coll(filename, info, this->traits_);
    assert(nh == nbheader && "Expected header size not written");
    ntot = nh;

    iret = MPI_File_open(comm_, filename.c_str(), MPI::MODE_CREATE|MPI::MODE_WRONLY, MPI::INFO_NULL, &fh);
    assert(iret == MPI_SUCCESS && "MPI_File_open failure");


    // Cycle over all fields, and write:
    if ( !this->traits_.multivar ) { // print each comp to sep. file
        assert(u[0]->size() > 0 && "Invalid state component");
        disp = nbheader ;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<Ftype>(), mpi_state_type_, "native", MPI::INFO_NULL);
        assert(iret == MPI_SUCCESS);
        iret = MPI_File_write_all(fh, u[0]->data(), u[0]->size(), T2GCDatatype<Ftype>(), &status);
        assert(iret == MPI_SUCCESS);
        MPI_Get_count(&status, MPI_BYTE, &nh);  
        ntot += nh;
    }
    else {                       // print each comp to same file
      for ( auto j=0; j<u.size(); j++ ) {
        disp = nbheader + j*nbgdof_ ; 
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<Ftype>(), mpi_state_type_, "native", MPI::INFO_NULL);
        assert(iret == MPI_SUCCESS);
        iret = MPI_File_write_all(fh, u[j]->data(), u[j]->size(), T2GCDatatype<Ftype>(), &status);
        assert(iret == MPI_SUCCESS);
        MPI_Get_count(&status, MPI_BYTE, &nh);  
        ntot += nh;
      }
    } 

    MPI_File_close(&fh);
#endif

    return ntot;

} // end, write_coll


//**********************************************************************************
//**********************************************************************************
// METHOD : read_coll
// DESC   : Collective read of state components
// ARGS   : filename: filename
//          info    : StateInfo structure
//          u       : state
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: none
//**********************************************************************************
template<typename Types>
GSIZET GIO<Types>::read_coll(GString filename, StateInfo &info, State &u, bool bstate)
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_MPI)

    GString        serr = "read_coll: ";
    GINT           iret;
    GSIZET         nb, nbheader, nh, ntot;
    Traits         ttraits=this->traits_;
    MPI_Offset     disp;
    MPI_File       fh;
    MPI_Status     status;


    // Read header and do some checks:
    nh = read_header(filename, info, ttraits);
    assert(ttraits.ivers == this->traits_.ivers
                                         && "Incompatible file version number");
    assert(ttraits.dim   == GDIM         && "File dimension incompatible with GDIM");
    assert(info   .gtype == this->grid_->gtype() 
                                           && "File grid type incompatible with grid");

    ntot = nh;

    if ( !bstate ) return ntot;

    // NOTE: read_header closes its file handle, so, when reopened, filepointer
    //       starts at displacement of 0"

    // Re-open file, and read state data:
    iret = MPI_File_open(comm_, filename.c_str(), MPI::MODE_RDONLY, MPI::INFO_NULL, &fh);
    assert(iret == MPI_SUCCESS && "MPI_File_open failure");
    nbheader = sz_header(info, this->traits_);


    // Cycle over all fields, and read:
    //   Note: any variable polynomial order element-by-element
    //         should be handled in ::update_type method via 
    //         state_disp_ & mpi_state_type_. Currently,
    //         variable order is not fully supported on read:
    if ( !this->traits_.multivar ) { // print each comp to sep. file
        assert(u[0]->size() > 0 && "Invalid state component");
        disp = nbheader ;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<Ftype>(), mpi_state_type_, "native", MPI::INFO_NULL);
        assert(iret == MPI_SUCCESS);
        iret = MPI_File_read_all(fh, u[0]->data(), u[0]->size(), T2GCDatatype<Ftype>(), &status);
        assert(iret == MPI_SUCCESS);
        ntot += iret == MPI_SUCCESS ? u[0]->size() : 0;
    }
    else {                       // read each comp from same file
      for ( auto j=0; j<u.size(); j++ ) {
        disp = nbheader + j*nbgdof_ ; 
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<Ftype>(), mpi_state_type_, "native", MPI::INFO_NULL);
        assert(iret == MPI_SUCCESS);
        iret = MPI_File_read_all(fh, u[j]->data(), u[j]->size(), T2GCDatatype<Ftype>(), &status);
        assert(iret == MPI_SUCCESS);
        ntot += iret == MPI_SUCCESS ? u[j]->size() : 0;
      } 
    } 

    MPI_File_close(&fh);
#endif

    return ntot;

} // end, read_coll

