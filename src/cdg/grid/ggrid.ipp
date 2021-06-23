//==================================================================================
// Module       : ggrid
// Date         : 8/31/18 (DLR)
// Description  : GeoFLOW grid object. Is primarily a container for
//                a list of elements, as an indirection buffer
//                that breaks down elements by type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : ptree: property tree
//          b    : GNBasis
//          comm : communicator
// RETURNS: none
//**********************************************************************************
template<typename Types>
GGrid<Types>::GGrid(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,Ftype>*> &b, GC_COMM &comm)
:
bInitialized_                   (FALSE),
bapplybc_                       (FALSE),
do_face_normals_                 (TRUE),
bpconst_                         (TRUE),
gderivtype_                (GDV_CONSTP),
do_gbdy_test_                   (FALSE),
doQDealias_                     (FALSE),
bInitQDealias_                  (FALSE),
nprocs_        (GComm::WorldSize(comm)),
ngelems_                            (0),
irank_         (GComm::WorldRank(comm)),
minnodedist_                      (0.0),
eps_                              (0.0),
volume_                           (0.0),
ivolume_                          (0.0),
comm_                            (comm),
mass_                         (NULLPTR),
imass_                        (NULLPTR),
ggfx_                         (NULLPTR),
ptree_                          (ptree),
bdy_apply_callback_           (NULLPTR)
{
       GEOFLOW_TRACE();

  GString snorm;
  GString gname          = ptree.getValue<GString>("grid_type");
  GString tname          = ptree.getValue<GString>("terrain_type","");
  PropertyTree gridptree = ptree.getPropertyTree(gname);

  // Get solver traits for terrain:
  cgtraits_.maxit = gridptree.getValue<GDOUBLE>("maxit", 128);
  cgtraits_.tol   = gridptree.getValue<GDOUBLE>("tol", 1.0e-8);
  snorm           = gridptree.getValue<GString>("norm_type", "GCG_NORM_INF");
  cgtraits_.normtype = LinSolverBase<CGTypePack>::str2normtype(snorm);

  cudat_.nstreams = ptree.getValue<GINT>("nstreams",1);
  cudat_.nstreams = MAX(cudat_.nstreams,1);
  
  
  doQDealias_ = ptree.getValue<GBOOL>("do_dealiasing", false);
  if ( doQDealias_ &&  ptree.keyExists("qdealias_order") ) {
     pqdealias_      = ptree.getArray<GINT>("qdealias_order");
     assert( pqdealias_.size() >= GDIM );
  }
  if ( !doQDealias_ ) EH::displayWarning("Quadratic dealiasing will not be done!");

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
GGrid<Types>::~GGrid()
{
       GEOFLOW_TRACE();

  if ( mass_ != NULLPTR ) delete mass_;
  if ( imass_ != NULLPTR ) delete imass_;
  for ( auto j=0; j<gelems_.size(); j++ ) {
    if ( gelems_[j] != NULLPTR ) delete gelems_[j];
  }

  GCBLAS::handle_destroy(cudat_.hcublas );
  GCBLAS::handle_destroy(cudat_.hbatch_cublas);
  for ( auto j=0; j<cudat_.pStream.size(); j++ ) {
    GCBLAS::stream_destroy(cudat_.pStream[j]);
  }

} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : do_typing
// DESC   : Do classification of element list into its
//          types, and store in member data structure. Each
//          type element in structure is an index into gelems_
//          corresp to that type
//          array.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::do_typing()
{
       GEOFLOW_TRACE();

  GTVector<GElemType> itmp(gelems_.size());

  GSIZET *ind=NULLPTR;
  GSIZET  nd=0;
  GSIZET  nfound;
  for ( auto j=0; j<gelems_.size(); j++ ) itmp[j] = gelems_[j]->elemtype();

  itype_.resize(GE_MAX);
  ntype_.resize(GE_MAX);
  ntype_ = 0;
  for ( auto j=0; j<GE_MAX; j++ ) {
    nfound = itmp.contains(static_cast<GElemType>(j),ind,nd);
    for ( auto i=0; i<nfound; i++ ) itype_[j].push_back(ind[i]);
    ntype_[j] = nfound;
  }
  if ( ind != NULLPTR ) delete [] ind;

  // Do sanity check:
  assert(ntype_.sum() == gelems_.size() && "Element typing failed");
} // end of method do_typing


//**********************************************************************************
//**********************************************************************************
// METHOD : print
// DESC   : Write basic grid info to specified file for plotting
// ARGS   : filename : string containing file name 
//          bdof     : print internal dof too?
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::print(GString filename, GBOOL bdof)
{
       GEOFLOW_TRACE();
  GString serr = "GridIcos::print: ";
  std::ofstream ios;

  GSIZET n;
  GTVector<GINT> N;
  GTPoint<Ftype> *vert;
  GTVector<GTVector<Ftype>> *xnodes;

  ios.open(filename);

  if ( !bdof ) { // print only wire mesh (vertices)
    for ( auto i=0; i<gelems_.size(); i++ ) { // for each element
       
       for ( auto j=0; j<gelems_[i]->nvertices(); j++ ) { // for each vertex of element
        vert = &gelems_[i]->xVertices(j);
        for ( auto k=0; k<vert->size()-1; k++ ) {
            ios << (*vert)[k] << " ";
        }
        ios << (*vert)[vert->size()-1] << std::endl;
         
      }
    }
    return;
  }

  // Print internal dofs too:
  for ( auto i=0; i<gelems_.size(); i++ ) { // for each element
    xnodes  = &gelems_[i]->xNodes();  
    N       = gelems_[i]->size();
    n       = (*xnodes)[0].size();
    if ( GDIM == 1 ) {
      for ( auto k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " " << std::endl;
      }
    }
    else if ( GDIM == 2 && gelems_[i]->elemtype() != GE_2DEMBEDDED ) {
      xnodes  = &gelems_[i]->xNodes();  
      for ( auto k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " "
            << (*xnodes)[1][k] << std::endl;
      }
    }
    else if (GDIM==2 && gelems_[i]->elemtype() == GE_2DEMBEDDED ) {
      // Lay down in separate 'sub-quads':
      for ( auto k=0; k<N[1]-1; k++ ) {
        for ( auto j=0; j<N[0]-1; j++ ) {
          // For each sub-quad, print its vertices:
          ios << (*xnodes)[0][j+k*N[0]] << " "
              << (*xnodes)[1][j+k*N[0]] << " "
              << (*xnodes)[2][j+k*N[0]] << std::endl;
          ios << (*xnodes)[0][j+1+k*N[0]] << " "
              << (*xnodes)[1][j+1+k*N[0]] << " "
              << (*xnodes)[2][j+1+k*N[0]] << std::endl;
          ios << (*xnodes)[0][j+1+(k+1)*N[0]] << " "
              << (*xnodes)[1][j+1+(k+1)*N[0]] << " "
              << (*xnodes)[2][j+1+(k+1)*N[0]] << std::endl;
          ios << (*xnodes)[0][j+(k+1)*N[0]] << " "
              << (*xnodes)[1][j+(k+1)*N[0]] << " "
              << (*xnodes)[2][j+(k+1)*N[0]] << std::endl;
        }
      }
    }
    else if ( GDIM == 3 ) {
      for ( auto k=0; k<n-1; k++ ) {
        ios << (*xnodes)[0][k] << " "
            << (*xnodes)[1][k] << " "
            << (*xnodes)[2][k] << std::endl;
      }
    }
  }


  ios.close();

} // end of method print


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD :  << operator method (1)
// DESC   : output stream operator
// ARGS   :
// RETURNS: ostream &
//**********************************************************************************
std::ostream &operator<<(std::ostream &str, GGrid &e)
{
       GEOFLOW_TRACE();
  return str;
} // end of operator <<
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : ndof
// DESC   : Find number of dof in grid
// ARGS   : none
// RETURNS: GSIZET number local dof
//**********************************************************************************
template<typename Types>
GSIZET GGrid<Types>::ndof()
{
       GEOFLOW_TRACE();
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET Ntot=0;
   for ( auto i=0; i<gelems_.size(); i++ ) Ntot += gelems_[i]->nnodes();

   return Ntot;
} // end of method ndof


//**********************************************************************************
//**********************************************************************************
// METHOD : nfacedof
// DESC   : Find number of elem face (not bdy) dof in grid. 
// ARGS   : none
// RETURNS: GSIZET number surface/face  dof
//**********************************************************************************
template<typename Types>
GSIZET GGrid<Types>::nfacedof()
{
       GEOFLOW_TRACE();
   return gieface_.size();
} // end of method nfacedof


//**********************************************************************************
//**********************************************************************************
// METHOD : nbdydof
// DESC   : Find number of bdy dof in grid. These will include
//          global domain boundary _and_ may also include embedded
//          boundary surface nodes.
// ARGS   : none
// RETURNS: GSIZET number surface dof
//**********************************************************************************
template<typename Types>
GSIZET GGrid<Types>::nbdydof()
{
       GEOFLOW_TRACE();
   return igbdy_.size();;
} // end of method nbdydof


//**********************************************************************************
//**********************************************************************************
// METHOD : minlength
// DESC   : Find elem lengths for each element
// ARGS   : pmin: poly. order corresponding to min element length
//          dx  : vector that gives the min length for each element;
//                filled if non-NULL
// RETURNS: Ftype separation
//**********************************************************************************
template<typename Types>
typename Types::Ftype GGrid<Types>::minlength(GTVector<GINT> *pmin, GTVector<Ftype> *dx)
{
       GEOFLOW_TRACE();
   assert(gelems_.size() > 0 && "Elements not set");

   GINT           pp;
   Ftype          dn, emin, lmin, gmin;
   GTPoint<Ftype> dr;
   GTVector<GTPoint<Ftype>> *xverts;
  
   if ( dx   != NULLPTR ) dx  ->resizem(gelems_.size());
   if ( pmin != NULLPTR ) pmin->resizem(gelems_.size());

   emin = lmin = std::numeric_limits<Ftype>::max();
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       dn = dr.norm();
       if ( dn < emin ) {
         emin = dn;
         pp = gelems_[i]->order(j % 2);
       }
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       dn = dr.norm();
       if ( dn < emin ) {
         emin = dn;
         pp = gelems_[i]->order(j % 2);
       }
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       dn = dr.norm();
       if ( dn < emin ) {
         emin = dn;
         pp = gelems_[i]->order(j % 2);
       }
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       dn = dr.norm();
       if ( dn < emin ) {
         emin = dn;
         pp = gelems_[i]->order(2);
       }
     }
     #endif
     lmin = MIN(lmin,emin);
     if ( dx   != NULLPTR ) (*dx)  [i] = emin; 
     if ( pmin != NULLPTR ) (*pmin)[i] = pp; 
   }

   GComm::Allreduce(&lmin, &gmin, 1, T2GCDatatype<Ftype>() , GC_OP_MIN, comm_);

   return gmin;
} // end of method minlength


//**********************************************************************************
//**********************************************************************************
// METHOD : maxlength
// DESC   : Find max elem length 
// ARGS   : 
// RETURNS: Ftype length
//**********************************************************************************
template<typename Types>
typename Types::Ftype GGrid<Types>::maxlength()
{
       GEOFLOW_TRACE();
   assert(gelems_.size() > 0 && "Elements not set");

   Ftype emax, lmax, gmax;
   GTPoint<Ftype> dr;
   GTVector<GTPoint<Ftype>> *xverts;

   lmax = 0.0;
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       emax = MAX(emax,dr.norm());
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       emax = MAX(emax,dr.norm());
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       emax = MAX(emax,dr.norm());
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       emax = MAX(emax,dr.norm());
     }
     #endif
     lmax = MAX(lmax,emax);
   }

   GComm::Allreduce(&lmax, &gmax, 1, T2GCDatatype<Ftype>() , GC_OP_MAX, comm_);

   return gmax;
} // end of method maxlength


//**********************************************************************************
//**********************************************************************************
// METHOD : avglength
// DESC   : Find avg elem length 
// ARGS   : none
// RETURNS: Ftype length
//**********************************************************************************
template<typename Types>
typename Types::Ftype GGrid<Types>::avglength()
{
       GEOFLOW_TRACE();
   assert(gelems_.size() > 0 && "Elements not set");

   Ftype gavg, lavg, navg, lv[2], gv[2];;
   GTPoint<Ftype> dr;
   GTVector<GTPoint<Ftype>> *xverts;

   navg = 0.0;
   lavg = 0.0;
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     #endif
   }

   lv[0] = lavg;
   lv[1] = navg;
   GComm::Allreduce(&lv, &gv, 2, T2GCDatatype<Ftype>() , GC_OP_SUM, comm_);
   gavg = gv[0] / gv[1];

   return gavg;
} // end of method avglength


//**********************************************************************************
//**********************************************************************************
// METHOD : grid_init (1)
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::grid_init()
{

       GEOFLOW_TRACE();
  do_elems(); // generate element list from derived class

  bpconst_ = ispconst();

  GComm::Synch(comm_);

  do_typing(); // do element-typing check


  // Have elements been set yet?
  assert(gelems_.size() > 0 && "Elements not set");

  // Restrict grid to a single element type:
  assert(ntype_.multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");

  if      ( itype_[GE_2DEMBEDDED].size() > 0 ) gtype_ = GE_2DEMBEDDED;
  else if ( itype_[GE_DEFORMED]  .size() > 0 ) gtype_ = GE_DEFORMED;
  else if ( itype_[GE_REGULAR]   .size() > 0 ) gtype_ = GE_REGULAR;

  globalize_coords    (); // set glob vec of node coords

  minnodedist_ = find_min_dist();
  find_min_dist(dxmin_);
  eps_         = 0.125*minnodedist_;


  // Create metric data, Jacobians:
  if ( gtype_ == GE_2DEMBEDDED || gtype_  == GE_DEFORMED ) {
    def_geom_init();
  }

  if ( gtype_ == GE_REGULAR ) {
    reg_geom_init();
  }

  // Configure and initialize domain bdy data:
  init_bc_info();


  // Get global number of elements:
  GSIZET nelems = gelems_.size();
  GComm::Allreduce(&nelems, &ngelems_, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);


  // Set element id vector:
  set_elemids();

  bInitialized_ = TRUE;

  mass_ = new Mass(*this);
  

  // Compute (global) grid volume:
  GTVector<Ftype> tmp0(ndof()), tmp1(ndof());
  tmp0 = 1.0;
  volume_  = integrate(tmp0, tmp1);
  ivolume_ = 1.0 / volume_;

  // Fill cuMatBlockDat structure:
  GCBLAS::handle_create(cudat_.hcublas );
  GCBLAS::handle_create(cudat_.hbatch_cublas);
  cudat_.ibblk.resize(cudat_.nstreams);
  cudat_.ieblk.resize(cudat_.nstreams);
  cudat_.pStream.resize(cudat_.nstreams);
  cudat_.nbatch = nelems;
  GSIZET idel = nelems/cudat_.nstreams;
  for ( auto j=0; j<cudat_.nstreams; j++ ) {
    GCBLAS::stream_create(cudat_.pStream[j]);
    cudat_.ibblk[j] = j*idel;
    cudat_.ieblk[j] = j<cudat_.nstreams-1 ? cudat_.ibblk[j] + idel
                    : nelems-1;
  }
  
} // end of method grid_init (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : grid_init (2)
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type. Called for restart.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::grid_init(GTMatrix<GINT> &p,
                      GTVector<GTVector<Ftype>> &xnodes)
{
       GEOFLOW_TRACE();

  do_elems(p, xnodes); // generate element list from derived class

  GComm::Synch(comm_);

  do_typing(); // do element-typing check

  // Have elements been set yet?
  assert(gelems_.size() > 0 && "Elements not set");


  // Restrict grid to a single element type:
  assert(ntype_.multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");


  if      ( itype_[GE_2DEMBEDDED].size() > 0 ) gtype_ = GE_2DEMBEDDED;
  else if ( itype_[GE_DEFORMED]  .size() > 0 ) gtype_ = GE_DEFORMED;
  else if ( itype_[GE_REGULAR]   .size() > 0 ) gtype_ = GE_REGULAR;

  globalize_coords    (); // set glob vec of node coords

  minnodedist_ = find_min_dist();
  find_min_dist(dxmin_);
  eps_         = 0.125*minnodedist_;

  if ( itype_[GE_2DEMBEDDED].size() > 0
    || itype_  [GE_DEFORMED].size() > 0 ) {
    def_geom_init();
  }

  if ( itype_[GE_REGULAR].size() > 0 ) {
    reg_geom_init();
  }

  // Configure and initialize domain bdy data:
  init_bc_info();


  bInitialized_ = TRUE;

  // Compute grid volume:
  GTVector<Ftype> tmp0(ndof()), tmp1(ndof());
  tmp0 = 1.0;
  volume_  = integrate(tmp0, tmp1);
  ivolume_ = 1.0 / volume_;

} // end of method grid_init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_geom_init
// DESC   : Initialize global (metric) variables for deformed elems. 
//          All elements are assumed to be of the same type. Note:
//          do_elems and globalize_coords must be called prior to entry.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::def_geom_init()
{
       GEOFLOW_TRACE();

   assert(gelems_.size() > 0 && "Elements not set");
   assert(gtype_ == GE_2DEMBEDDED
       || gtype_ == GE_DEFORMED && "Invalid element type");


   GString serr = "GGrid<Types>::def_geom_init: ";
   GSIZET nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
   GTVector<GTVector<Ftype>> *xe;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,nxy);
   dXdXi_.resize(nxy,nxy);
   for ( auto j=0; j<nxy; j++ ) {
     for ( auto i=0; i<nxy; i++ )  {
       dXidX_(i,j).resize(ndof());
       dXdXi_(i,j).resize(ndof());
     }
   }
   Jac_.resize(ndof());
   faceJac_.resize(nfacedof());

   // Resize surface-point-wise normals, tangents:
   faceNormals_.resize(nxy); // no. coords for each normal at each face point
   bdyNormals_.resize(nxy); // no. coords for each normal at each domain bdy point
   for ( auto i=0; i<bdyNormals_.size(); i++ ) {
     faceNormals_[i].resize(nfacedof());
     bdyNormals_ [i].resize(nbdydof());
   }

   bdyTangents_.resize(GDIM-1);
   for ( auto j=0; j<bdyTangents_.size(); j++ ) {
     bdyTangents_[j].resize(nxy);
     for ( auto i=0; i<bdyTangents_.size(); i++ ) {
       bdyTangents_[j][i].resize(nbdydof());
     }
   }


   // Now, set the geometry/metric quanties from the elements:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

     xe    = &gelems_[e]->xNodes();

     // Restrict global arrays to local scope:
     Jac_.range(ibeg, iend);
     for ( auto j=0; j<nxy; j++ ) {
       for ( auto i=0; i<nxy; i++ )  {
         dXidX_(i,j).range(ibeg, iend);
         dXdXi_(i,j).range(ibeg, iend);
       }
     }

     // Set the geom/metric quantities using element data:
     // (Note: faceJac_ is not modified in this call)
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(dXdXi_, dXidX_, Jac_, faceJac_);
     }
     else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(dXdXi_, dXidX_, Jac_, faceJac_);
     }

     // Zero-out local xe; only global allowed now:
//   for ( auto j=0; j<nxy; j++ ) (*xe)[j].clear(); 
     
   } // end, element loop

   Jac_.range_reset();
   for ( auto j=0; j<nxy; j++ )  {
     for ( auto i=0; i<nxy; i++ )  {
       dXidX_(i,j).range_reset();
       dXdXi_(i,j).range_reset();
     }
   }

   do_face_data();

   GComm::Synch(comm_);
   
} // end of method def_geom_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_geom_init
// DESC   : Initialize global (metric) variables for regular elemes. 
//          All elements are assumed to be of the same type. Note:
//          do_elems and globalize_coords must be called prior to entry.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::reg_geom_init()
{
       GEOFLOW_TRACE();
   assert(gelems_.size() > 0 && "Elements not set");
   assert( gtype_ == GE_REGULAR && "Invalid element type");

   GString serr = "GridIcos::reg_geom_init: ";
   GSIZET nxy = GDIM;
   GTVector<GTVector<Ftype>> *xe;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,1);
   dXdXi_.resize(nxy,1);
   for ( auto i=0; i<nxy; i++ ) {
     dXidX_(i,0).resize(ndof());
     dXdXi_(i,0).resize(ndof());
   }
   Jac_.resize(ndof());
   faceJac_.resize(nfacedof());


   // Resize surface-point-wise normals, tangents:
   faceNormals_.resize(nxy); // no. coords for each normal at each face point
   bdyNormals_ .resize(nxy); // no. coords for each normal at each bdy point
   for ( auto i=0; i<nxy; i++ ) {
     faceNormals_[i].resize(nfacedof());
     bdyNormals_ [i].resize(nbdydof());
   }

   bdyTangents_.resize(GDIM-1);
   for ( auto j=0; j<bdyTangents_.size(); j++ ) {
     bdyTangents_[j].resize(nxy);
     for ( auto i=0; i<bdyTangents_[j].size(); i++ ) {
       bdyTangents_[j][i].resize(nbdydof());
     }
   }


   // Now, set the geometry/metric quanties from the elements:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     xe    = &gelems_[e]->xNodes();
     Jac_.range(ibeg, iend);
     for ( auto j=0; j<dXidX_.size(2); j++ ) {
       for ( auto i=0; i<dXidX_.size(1); i++ )  {
         dXidX_(i,j).range(ibeg, iend);
         dXdXi_(i,j).range(ibeg, iend);
       }
     }

     // Set the geom/metric quantities using element data:
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(dXdXi_, dXidX_, Jac_, faceJac_); 
     } 
     else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(dXdXi_, dXidX_, Jac_, faceJac_);
     }
      
     // Zero-out local xe; only global allowed now:
//   for ( auto j=0; j<nxy; j++ ) (*xe)[j].clear(); 

   } // end, element loop

   Jac_.range_reset();
   for ( auto j=0; j<dXidX_.size(2); j++ )  {
     for ( auto i=0; i<dXidX_.size(1); i++ )  {
       dXidX_(i,j).range_reset();
       dXdXi_(i,j).range_reset();
     }
   }

   do_face_data();

   GComm::Synch(comm_);
   
} // end of method reg_geom_init



//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_data
// DESC   : Compute normals to elem faces, and to domain boundary 
//          nodes. Note: init_bc_info must be called prior to entry.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::do_face_data()
{
       GEOFLOW_TRACE();

  assert(gelems_.size() > 0 && "Elements not set");

  GString         serr = "GridIcos::do_face_data: ";


  // Set element face normals. Note: arrays for 
  // normals are allocated in these calls:
  if ( do_face_normals_ ) {
    elem_face_data(dXdXi_, gieface_, faceMass_, faceNormals_);
  }

} // end of method do_face_data


//**********************************************************************************
//**********************************************************************************
// METHOD : globalize_coords
// DESC   : Create 'global' vectors from element-based coordinates
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::globalize_coords()
{
       GEOFLOW_TRACE();
   GString serr = "GridIcos::globalize_coords: ";
   GSIZET  nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
   GTVector<GTVector<Ftype>> *xe;

   xNodes_.resize(nxy);
   for ( auto j=0; j<nxy; j++ ) xNodes_[j].resize(ndof());

   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     xe    = &gelems_[e]->xNodes();

     // Set global nodal Cart coords from element coords:
     for ( auto j=0; j<nxy; j++ ) {
       xNodes_[j].range(ibeg, iend);
       xNodes_[j] = (*xe)[j];
     }

   } // end, element loop

   // Reset global scope:
   for ( auto j=0; j<nxy; j++ ) xNodes_[j].range_reset();


} // end, method globalize_coords


//**********************************************************************************
//**********************************************************************************
// METHOD : dXidX (1)
// DESC   : return global metric terms
// ARGS   : none
// RETURNS: GTMatrix<GTVector<Ftype>> &
//**********************************************************************************
template<typename Types>
GTMatrix<GTVector<typename Types::Ftype>> &GGrid<Types>::dXidX()
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   return dXidX_;

} // end of method dXidX


//**********************************************************************************
//**********************************************************************************
// METHOD : dXidX (2)
// DESC   : return global metric element
// ARGS   : i,j : matrix element indices
// RETURNS: GTVector<Ftype> &
//**********************************************************************************
template<typename Types>
GTVector<typename Types::Ftype> &GGrid<Types>::dXidX(GSIZET i, GSIZET j)
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   return dXidX_(i,j);

} // end of method dXidX


//**********************************************************************************
//**********************************************************************************
// METHOD : massop
// DESC   : return global diagonal mass operator
// ARGS   : none
// RETURNS: GTVector<Ftype> &
//**********************************************************************************
template<typename Types>
typename Types::Mass &GGrid<Types>::massop()
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   return *mass_;

} // end of method massop


//**********************************************************************************
//**********************************************************************************
// METHOD : imassop
// DESC   : return inverse of global diagonal mass operator
// ARGS   : none
// RETURNS: GTVector<Ftype> &
//**********************************************************************************
template<typename Types>
typename Types::Mass &GGrid<Types>::imassop()
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   if ( imass_ == NULLPTR ) imass_ = new Mass(*this, TRUE);
   return *imass_;

} // end of method imassop


//**********************************************************************************
//**********************************************************************************
// METHOD : Jac
// DESC   : return global coord transform Jacobian
// ARGS   : none
// RETURNS: GTVector<Ftype> &
//**********************************************************************************
template<typename Types>
GTVector<typename Types::Ftype> &GGrid<Types>::Jac()
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   return Jac_;

} // end of method Jac


//**********************************************************************************
//**********************************************************************************
// METHOD : faceJac
// DESC   : Return global coord transform Jacobian for faces
// ARGS   : none
// RETURNS: GTVector<Ftype> &
//**********************************************************************************
template<typename Types>
GTVector<typename Types::Ftype> &GGrid<Types>::faceJac()
{
       GEOFLOW_TRACE();
   assert(bInitialized_ && "Object not inititaized");
   return faceJac_;

} // end of method faceJac


//**********************************************************************************
//**********************************************************************************
// METHOD : find_min_dist (1)
// DESC   : Compute min inter-node distance for grid
//          NOTE: this method isn't recommended for regular
//                use, as it's slow.
// ARGS   : none
// RETURNS: minimum nodal dist over entire grid
//**********************************************************************************
template<typename Types>
typename Types::Ftype GGrid<Types>::find_min_dist()
{
       GEOFLOW_TRACE();
  assert(gelems_.size() > 0 && "Elements not set");

  Ftype            xmin, xgmin;
  GTVector<GINT>   pmin;
  GTVector<Ftype>  dx;
 
#if 0
  Ftype           tiny = 1000.0*std::numeric_limits<Ftype>::epsilon();
  GTPoint<Ftype>  dx(3), p0(3);

  // Find min node distance on grid:
  dx = 0.0;
  p0 = 0.0;

  Ftype xmin = std::numeric_limits<Ftype>::max();
  Ftype xgmin, xn;

  GComm::Synch(comm_);

  // Get grid dimension:
  GINT nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  // Cycle over all node points and find min distance:
  // NOTE: this isn't fast....
  for ( auto j=0; j<xNodes_[0].size()-1; j++ ) {
    for ( auto i=0; i<nxy; i++ ) {
      p0[i] = xNodes_[i][j];
      dx[i] = xNodes_[i][j+1];
    }
    dx -= p0; 
    xn = dx.norm();
    if ( xn > tiny ) xmin = MIN(xmin, xn); 
  }
#else
  minlength(&pmin, &dx);
  for ( auto j=0; j<dx.size(); j++ ) {
    dx[j] *= 1.0 / (pmin[j]*pmin[j]);
  }
  xmin = dx.min();
#endif
  GComm::Allreduce(&xmin, &xgmin, 1, T2GCDatatype<Ftype>() , GC_OP_MIN, comm_);

  return xgmin;

} // end of method find_min_dist (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : find_min_dist (2)
// DESC   : Compute min inter-node distance for each element on grid
//          NOTE: this method isn't recommended for regular
//                use, as it's slow.
// ARGS   : dx  : array of smallest nodal distance for each element
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::find_min_dist(GTVector<Ftype> &dx)
{
       GEOFLOW_TRACE();
  assert(gelems_.size() > 0 && "Elements not set");

  GTVector<GINT>   pmin;
 
#if 0
  GSIZET           ibeg, iend;
  GTPoint<Ftype>   dd(3), p0(3);

  // Find min node distance on grid:
  dd = 0.0;
  p0 = 0.0;

  Ftype tiny = 1000.0*std::numeric_limits<Ftype>::epsilon();
  Ftype xmin; 
  Ftype xn;

  dx.resizem(gelems_.size());

  // Cycle over all node points and find min distance:
  // NOTE: this isn't fast....
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    xmin  = std::numeric_limits<Ftype>::max();
    for ( auto j=ibeg; j<iend; j++ ) {
      for ( auto i=0; i<xNodes_.size(); i++ ) {
        p0[i] = xNodes_[i][j];
        dd[i] = xNodes_[i][j+1];
      }
      dd -= p0; 
      xn = dd.norm();
      if ( xn > tiny ) xmin  = MIN(xmin, xn); 
    }
    dx[e] = xmin;
  }
#else
  minlength(&pmin, &dx);
  for ( auto j=0; j<dx.size(); j++ ) {
    dx[j] *= 1.0 / (pmin[j]*pmin[j]);
  }
#endif

} // end of method find_min_dist (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : integrate
// DESC   : Compute spatial integral of global input vector. Result 
//          is a sum over all MPI tasks. NOTE: This method extracts weights
//          from the elements, so don't use this method very often.
// ARGS   : u     : 'global' integral argument
//          tmp   :  tmp vector, same size as u
//          blocal:  do global reduction?
// RETURNS: Ftype integral
//**********************************************************************************
template<typename Types>
typename Types::Ftype GGrid<Types>::integrate(GTVector<Ftype> &u, GTVector<Ftype> &tmp, GBOOL bglobal)
{
       GEOFLOW_TRACE();
  assert(bInitialized_ && "Object not inititaized");

  GSIZET                       ibeg, iend; // beg, end indices for global array
  GSIZET                       n;
  Ftype                        xint, xgint;
  GTVector<GINT>               N(GDIM);    // coord node sizes
  GTVector<GTVector<Ftype>*>   W(GDIM);    // element weights

  // NOTE: We could instantiate a mass matrix, but this would
  //       require computing a new GDIM-operator each time
  //       integration is required or setting it.
#if defined(_G_IS2D)
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
    u.range(ibeg, iend);

    for ( auto k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( auto k=0; k<N[1]; k++ ) {
      for ( auto j=0; j<N[0]; j++ ) {
        tmp[n] = (*W[1])[k]*(*W[0])[j] * u[n];
        n++;
      }
    }
  } // end, element loop
#elif defined(_G_IS3D)
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
    u.range(ibeg, iend);

    for ( auto k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( auto k=0; k<N[2]; k++ ) {
      for ( auto j=0; j<N[1]; j++ ) {
        for ( auto i=0; i<N[0]; i++ ) {
          tmp[n] = (*W[2])[k]*(*W[1])[j]*(*W[0])[i] * u[n];
          n++;
        }
      }
    }
  } // end, element loop
#endif
  tmp.range_reset(); 
  u.range_reset(); 

  // Multiply by Jacobian:
  tmp.pointProd(Jac_);
  xint = tmp.sum();  

  xgint = xint;
  if ( bglobal ) {
    GComm::Allreduce(&xint, &xgint, 1, T2GCDatatype<Ftype>() , GC_OP_SUM, comm_);
  }

  return xgint;

} // end of method integrate


//**********************************************************************************
//**********************************************************************************
// METHOD : deriv (1)
// DESC   : Compute (collocated) spatial derivative of u in direction idir, 
//          and return in du.
// ARGS   : u   : 'global' integral argument
//          idir: coord wrt which to take derivative (1, 2, or 3)
//          utmp: tmp vector of same size as u
//          du  : derivtive, returned
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::deriv(GTVector<Ftype> &u, GINT idir, GTVector<Ftype> &utmp, 
                  GTVector<Ftype> &du)
{
       GEOFLOW_TRACE();
  assert(bInitialized_ && "Object not inititialized");


  GINT nc =  gtype_ == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GTMatrix<GTVector<Ftype>> *dXidX = &this->dXidX();

  // du/dx_idir = Sum_j=[1:N] dXi_j/dX_idir D_j u:
  if ( this->gtype() == GE_REGULAR ) {
    assert(idir > 0 && idir <= GDIM && "Invalid spatial derivative");
    compute_grefderiv(u, etmp_, idir, FALSE, du); // D_idir u

    du.pointProd((*dXidX)(idir-1, 0));
  }
  else {  // compute dXi_j/dX_idir D^j u:
    assert(idir > 0 && idir <= nc && "Invalid spatial derivative");
    compute_grefderiv(u, etmp_, 1, FALSE, du); // D_xi u
    du.pointProd((*dXidX)(0,idir-1));
    for ( auto j=1; j<nc; j++ ) {
      compute_grefderiv(u, etmp_, j+1, FALSE, utmp); // D_xi^j u
//    utmp.pointProd((*dXidX)(j,idir-1));
      utmp.pointProd((*dXidX)(idir-1,j));
      du += utmp; 

    }
  }
    
} // end of method deriv (1)

//**********************************************************************************
//**********************************************************************************
// METHOD : wderiv 
// DESC   : Compute weak spatial derivative of u in direction idir, and
//          return in du. Allow user to select transpose of reference 
//          derivative matrix.
// ARGS   : u      : 'global' integral argument
//          idir   : coord wrt which to take derivative (1, 2, or 3)
//          dotrans: select transpose (TRUE) or not (FALSE)
//          utmp   : tmp vector of same size as u
//          du     : derivtive, returned
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::wderiv(GTVector<Ftype> &u, GINT idir, GBOOL dotrans, 
                   GTVector<Ftype> &utmp, GTVector<Ftype> &du)
{
  GEOFLOW_TRACE();
  assert(bInitialized_ && "Object not inititialized");


  GINT       nxy = this->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GTMatrix<GTVector<Ftype>> *dXidX = &this->dXidX();
  GTVector<Ftype> *mass = this->massop().data();

GTVector<Ftype> t1(ndof());

  assert(idir > 0 && idir <= nxy && "Invalid spatial derivative");

  // du/dx_idir = Sum_j=[1:N] dxi_j/dx_idir D_j u:
  if ( this->gtype() == GE_REGULAR ) {
    if ( dotrans ) {
#if 1
      u.pointProd((*dXidX)(idir-1,0), utmp);
      utmp.pointProd(*mass);
      compute_grefderiv(utmp, etmp_, idir, dotrans, du); // D_idir u
#else
      compute_grefderiv(u, etmp_, idir, dotrans, du); // D_idir u
      du.pointProd((*dXidX)(idir-1,0));
      du.pointProd(*mass);
#endif

    }
    else {
      compute_grefderiv(u, etmp_, idir, dotrans, du); // D_idir u
      du.pointProd((*dXidX)(idir-1, 0));
      du.pointProd(*mass);
    }
  }
  else {  // compute dXi_j/dX_idir D^j u:
    if ( dotrans ) {
      u.pointProd((*dXidX)(0,idir-1), utmp);
      utmp.pointProd(*mass);
      compute_grefderiv(utmp, etmp_, 1, dotrans, du); // D_xi u
      for ( auto j=1; j<GDIM; j++ ) {
        u.pointProd((*dXidX)(idir-1,j),t1);
        t1.pointProd(*mass);
        compute_grefderiv(t1, etmp_, j+1, dotrans, utmp); // D_xi^j u
        du += utmp; 
      }
    }
    else {
      compute_grefderiv(u, etmp_, 1, dotrans, du); // D_xi u
      du.pointProd((*dXidX)(0,idir-1));
      for ( auto j=1; j<GDIM; j++ ) {
        compute_grefderiv(u, etmp_, j+1, dotrans, utmp); // D_xi^j u
        utmp.pointProd((*dXidX)(idir-1,j));
        utmp.pointProd(*mass);
        du += utmp; 
      }
    }
  }
    
} // end of method wderiv 


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : wderiv
// DESC   : Compute _weak_ spatial derivative of a scalar, q,  in direction 
//          idir, and return in du. This is computed as:
//               -D^T_idir p
// ARGS   : u     : scalar field
//          idir  : coord wrt which to take derivative (1, 2, or 3)
//          bwghts: include Gaussian weither and Jacobian
//          utmp  : tmp vector of same size as u
//          du    : derivtive, returned
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::wderiv(GTVector<Ftype> &u, GINT idir, GBOOL bwghts, GTVector<Ftype> &utmp, 
                   GTVector<Ftype> &du)
{
  GEOFLOW_TRACE();
  assert(bInitialized_ && "Object not inititialized");


  GTMatrix<GTVector<Ftype>> *dXidX = &this->dXidX();
  GTVector<Ftype> *Jac = &this->Jac();
  GTVector<Ftype> *mass=  this->massop().data();


  // du/dx_idir = Sum_j=[1:N] dxi_j/dx_idir D_j u:
  if ( this->gtype() == GE_REGULAR ) {
    assert(idir > 0 && idir <= GDIM && "Invalid spatial derivative");
    compute_grefderivW(u, etmp_, idir, FALSE, du); 
    du.pointProd((*dXidX)(idir-1, 0));
  }
  else {  // compute dXi_j/dX_idir D^j u:
    assert(idir > 0 && idir <= GDIM+1 && "Invalid spatial derivative");
    compute_grefderivW(u, etmp_, 1, FALSE, du); 
    du.pointProd((*dXidX)(0,idir-1));
    for ( auto j=1; j<GDIM; j++ ) {
      compute_grefderivW(u, etmp_, j+1, FALSE, utmp); 
      utmp.pointProd((*dXidX)(idir-1,j ));
      du += utmp; 
    }
  }
  if ( bwghts ) {
    for ( auto i=0; i<du.size(); i++ ) {
      du[i] *= (*Jac)[i] ; // * (*mass)[i];
    }
  }
    
} // end of method wderiv
#endif


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : init_local_face_info
// DESC   : Set local face info from element face data,
//          to be called after elements have been set.
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::init_local_face_info()
{
  GEOFLOW_TRACE();
  GBOOL                        bret;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GINT>               itmp; 
  GTVector<GUINT>              utmp; 
  GTVector<GTVector<GINT>>    *ieface;     // element face indices
  GTVector<GTVector<GUINT>>   *deface;     // element face node description
  GTVector<Ftype>            *efacemass;  // element face weights
  GTVector<Ftype>             ftmp; 

  GSIZET  m, n, nn; 
  GSIZET        ig; // index into global array

  if ( !do_face_normals_ ) return;

  n = 0;
  for ( auto e=0; e<gelems_.size(); e++ ) { // get global # face nodes
    n += gelems_[e]->nfnodes();
  }
  gieface_.resize(n);
  faceMass_.resize(n);

  itmp.resize(n);
  utmp.resize(n);
  ftmp.resize(n);

  // Find list over all elemes of element face nodes:
  nn = 0; // global reference index
  m  = 0; // current num of global face indices
  for ( auto e=0; e<gelems_.size(); e++ ) { // cycle over all elems
    ibeg      = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    ieface    = &gelems_[e]->face_indices(); // set in child class
    deface    = &gelems_[e]->face_desc(); // set in child class
    efacemass = &gelems_[e]->face_mass();
    for ( auto j=0; j<ieface->size(); j++ ) { // cycle over all elem faces
      for ( auto k=0; k<(*ieface)[j].size(); k++ ) {
        ig = nn + (*ieface)[j][k];
#if 0
        if ( !gieface_.containsn(ig, m) ) { // don't include repeated face ind
          itmp  [m] = ig;
          utmp  [m] = (*deface)[j][k];
          ftmp  [m] = (*efacemass)[(*ieface)[j][k]];
          m++;
        }
#else
        itmp  [m] = ig;
        utmp  [m] = (*deface)[j][k];
        ftmp  [m] = (*efacemass)[(*ieface)[j][k]];
        m++;
#endif
      } // end, face node loop
    } // end, face loop
    nn += gelems_[e]->nnodes();
  } // end, element loop

  gieface_ .resize(m); 
  faceMass_.resize(m);
  for ( auto j=0; j<m; j++ ) {
    gieface_ [j] = itmp[j];
    faceMass_[j] = ftmp[j];
  }


} // end, init_local_face_info
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : init_bc_info
// DESC   : Set global bdy condition data from the element bdy data,
//          to be called after elements have been set.
// ARGS   : bterrain: tells us if we're doing terrain (TRUE) or not (FALSE)
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::init_bc_info(GBOOL bterrain)
{
       GEOFLOW_TRACE();
  GSIZET   nind;
  GBdyType btype;


  // Find boundary indices & types from config file 
  // specification, for _each_ natural/canonical domain face.
  // Bdy normals etc. are also computed internally:
  config_gbdy(ptree_, bterrain, igbdy_bdyface_, igbdyt_bdyface_, igbdy_, degbdy_);


  // Flatten bdy index indirection array; this is
  // done in child classes, and stored in igbdy_.

  // Create bdy type bins for each domain bdy:
  //   [Dom bdy][bdytype][volume index]:
  GBdyType         itype;
  GSIZET          *ind=NULLPTR;
  GSIZET           n, nbdy;
  GTVector<GSIZET> iunique;

//cout << "init_bc_info: igbdy=" << igbdy_ << endl;
   
  nind = 0; 
  // Compute 'binned' structures for global indices
  // defining bdys, and to which index in bdy normal vectors
  // a given bdy node corresponds:
  n = 0; // cycle over all bdy nodes
  igbdy_binned_.resize(igbdy_bdyface_.size());
  for ( auto k=0; k<igbdy_bdyface_.size(); k++ ) { // cycle over canonical bdy face
//cout << "init_bc_info: igbdy[" << k << "]=" << igbdy_bdyface_[k] << endl;
    nbdy = igbdyt_bdyface_[k].size();
    igbdy_binned_ [k].resize(GBDY_MAX);
    for ( auto j=0; j<GBDY_MAX; j++ ) { // cycle over each bc type
      itype = static_cast<GBdyType>(j);
      nbdy = igbdyt_bdyface_[k].multiplicity(itype, ind, nind);
      igbdy_binned_[k][j].resize(nbdy);
      for ( auto i=0; i<nbdy; i++ ) { // assign comp. volume index
        igbdy_binned_[k][j][i] = igbdy_bdyface_[k][ind[i]];
        n++;
      }
    } // end, bdy cond type loop
  } // end, can. bdy loop
  if ( ind != NULLPTR ) delete [] ind;

  // Compute mask matrix from bdy vector:
  mask_.resize(ndof()); 
  mask_ = 1.0;
  for ( auto k=0; k<igbdy_binned_.size(); k++ ) {
    for ( auto j=0; j<GBDY_MAX; j++ ) {
      btype = static_cast<GBdyType>(j);
      if ( btype != GBDY_NONE && btype != GBDY_PERIODIC ) {
        for ( auto i=0; i<igbdy_binned_[k][j].size(); i++ ) {
          mask_[igbdy_binned_[k][j][i]] = 0.0;
        }
      }
    }
  }

  // Compute node locations of bdy nodes:
  gbdyNodes_.resize(xNodes_.size());
  for ( auto i=0; i<gbdyNodes_.size(); i++ ) {
    gbdyNodes_[i].resize(igbdy_.size());
    for ( auto j=0; j<igbdy_.size(); j++ ) { 
      gbdyNodes_[i][j] = xNodes_[i][igbdy_[j]];
    }
  }

} // end of method init_bc_info


//**********************************************************************************
//**********************************************************************************
// METHOD : add_terrain
// DESC   : Add bdy terrain by deforming coordinates. 
//          Works with only non-regular grids, and must be called 
//          after 'base' grid is completely constructed.
// ARGS   : xb  : boundary coordinates
//          utmp: tmp array pool
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::add_terrain(const State &xb, State &utmp)
{
       GEOFLOW_TRACE();
   assert(gtype_ == GE_2DEMBEDDED
       || gtype_ == GE_DEFORMED && "Invalid element type");

   assert(utmp.size() >= 6);

  // Set up temp arrays from pool:
  GINT       iret;
  StateComp *x0=utmp[0];
  StateComp *b =utmp[1];
  State      tmp(utmp.size()-2);
  GTVector<GTVector<Ftype>> *xe;


  for ( auto j=0; j<tmp.size(); j++ ) tmp[j] = utmp[j+2];

  assert(ggfx_ != NULLPTR && "GGFX not set");

  // Construct solver and weak Laplacian operator:
  GCG<CGTypePack>   cg(cgtraits_, *this, *ggfx_, tmp);
  TerrainOp         H(*this);

  H.use_metric(TRUE); // do Laplacian in real space
  xb_.resize(xNodes_.size());

  // Solve Nabla^2 (Xnew + Xb ) = 0 
  // for new (homgogeneous) grid solution, Xnew, 
  // given terrain, Xb, and // 'base' grid, XNodes:
  GEOFLOW_TRACE_START("GGrid<Types>::add_terrain: Solve");
  for ( auto j=0; j<xNodes_.size(); j++ ) {
   *b  = 0.0;
   *x0 = 0.0; // first guess
    xb_[j].resize(xb[j]->size()); xb_[j] = *xb[j];
//if ( j == 1 )
//cout << "GGrid<Types>::add_terrain: xb[" << j << "]=" << *xb[j] << endl;
    iret = cg.solve(H, *b, *xb[j], *x0);

    assert(iret == GCG<CGTypePack>::GCGERR_NONE);
    xNodes_[j] = *x0;             // Reset XNodes = x0
//GPP(comm_,"GGrid<Types>::add_terrain: new_xNodes[" << j << "]=" << xNodes_[j]);
  }
  GEOFLOW_TRACE_STOP();
 
  // Before computing new metric, Jacobian, etc, must set
  // new coordinates in elements that have already been initialized:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     for ( auto j=0; j<xNodes_.size(); j++ ) xNodes_[j].range(ibeg, iend);
     gelems_[e]->set_nodes(xNodes_);
   }
   for ( auto j=0; j<xNodes_.size(); j++ ) xNodes_[j].range_reset();


  // Now, with new coordinates, recompute metric terms, 
  // Jacobian. Only deformed element grids can have terrain,
  // so only these are updated:
  def_geom_init();

  // Re-compute Mass matrices:
  if ( mass_  != NULLPTR ) delete mass_;
  if ( imass_ != NULLPTR ) delete imass_;
  mass_  = new Mass(*this, FALSE);
  imass_ = new Mass(*this, TRUE);


  // Configure and initialize domain bdy data:
  init_bc_info(TRUE); // TRUE ==> adding terrain

  // Compute element face data:
  do_face_data();

  // Compute new minimum node distance:
  minnodedist_ = find_min_dist();
  find_min_dist(dxmin_);

  // Compute new (global) grid volume:
  *utmp[0] = 1.0;
  volume_  = integrate(*utmp[0], *utmp[1]);
  ivolume_ = 1.0 / volume_;

} // end of method add_terrain


//**********************************************************************************
//**********************************************************************************
// METHOD : smooth
// DESC   :
//
// DESC   : Computes a weighted average.
//              Calculates:
//                u <- DSS(M_L u) / DSS(M_L),
//          where u is the field expressed in local format;
//          M_L is the local mass matrix (unassembled);
//          DSS is the application of Q Q^T, or the direct stiffness operator.
// ARGS   :
//          tmp  : tmp space
//          u    : Locally expressed field to smooth
// RETURNS: none.
//************************************************************************************
template<typename Types>
void GGrid<Types>::smooth(GTVector<Ftype> &tmp, GTVector<Ftype> &u)
{
  GEOFLOW_TRACE();

  tmp = u;

  u.pointProd(*(this->massop().data()));
  tmp = *(this->imassop().data());
  this->ggfx_->doOp(tmp, typename GGFX<Ftype>::Smooth());  // DSS(mass_local)

  u.pointProd(tmp);      // M_assembled u M_L

} // end, smooth method

//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivW
// DESC   : Compute weighted tensor product derivative in specified direction
//          of specified field, u, in ref space, using grid object.
//          Compute
//            du = [ Mz_X_My_X_MxDx, or
//                   Mz_X_MyDy_X_Mx, or
//                   MzDz_X_My_X_Mx].
//     
//          depending on whether idir = 1, 2, or 3, respectively,
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          idir   : coordinate direction (1, 2, or 3)
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          du     : vector of length of u containing the derivative.
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::compute_grefderivW(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                               GINT idir, GBOOL dotrans, GTVector<Ftype> &du)
{
       GEOFLOW_TRACE();
  GSIZET               ibeg, iend; // beg, end indices for global array
  GBOOL                bembedded;
  GTVector<GSIZET>     N(GDIM);
  GTVector<GTVector<Ftype>*>
                       W(GDIM);    // element weights
  GTMatrix<Ftype>     *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &this->elems();

  bembedded = this->gtype() == GE_2DEMBEDDED;

#if defined(_G_IS2D)
  switch (idir) {
  case 1:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      W[1]= (*gelems)[e]->gbasis(1)->getWeights();
      Di  = (*gelems)[e]->gbasis(0)->getDerivMatrixW(dotrans);
      GMTK::Dg2_X_D1(*Di, *W[1], u, etmp, du); 
    }
    break;
  case 2:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      W[0]= (*gelems)[e]->gbasis(0)->getWeights();
      Di  = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans);
      GMTK::D2_X_Dg1(*W[0], *Di, u, etmp, du); 
    }
    break;
  case 3:
    assert( GDIM == 3
         && "Only GDIM reference derivatives");
    du = 0.0; //u;
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset to global range
  du.range_reset();

#elif defined(_G_IS3D)

  switch (idir) {
  case 1:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      W[1]= (*gelems)[e]->gbasis(1)->getWeights();
      W[2]= (*gelems)[e]->gbasis(2)->getWeights();
      Di  = (*gelems)[e]->gbasis(0)->getDerivMatrixW(dotrans); 
      GMTK::Dg3_X_Dg2_X_D1(*Di, *W[1], *W[2], u, etmp, du); 
    }
    break;

  case 2:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      W[0]= (*gelems)[e]->gbasis(0)->getWeights();
      W[2]= (*gelems)[e]->gbasis(2)->getWeights();
      Di  = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans); 
      GMTK::Dg3_X_D2_X_Dg1(*W[0], *Di, *W[2], u, etmp, du); 
    }
    break;

  case 3:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      W[0]= (*gelems)[e]->gbasis(0)->getWeights();
      W[1]= (*gelems)[e]->gbasis(1)->getWeights();
      Di  = (*gelems)[e]->gbasis(2)->getDerivMatrix(!dotrans); 
      GMTK::D3_X_Dg2_X_Dg1(*W[0], *W[1], *Di, u, etmp, du); 
    }
    break;

  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset global vec to globalrange
  du.range_reset();

#endif

} // end of method compute_grefderivW


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivsW
// DESC   : Compute tensor product derivs of specified field, u, in ref space 
//          for grid using grid object to determine which to compute, and 
//          include weights.
//          Compute:
//            du = (Mz_X_My_Mx) [ I_X_I_X_Dx
//                                I_X_Dy_X_I
//                                Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions, and
//          Mi are the (diagonal) 1d-weights (or 'mass functions'). This can be 
//          re-written as
//            du = [ Mz_X_My_X_(Mx Dx)
//                   Mz_X_(My Dy)_X_Mx
//                  (Mz Dz)_X_My_X_Mx],
//           with comparable expressions for 2d.
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need two vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::compute_grefderivsW(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                GBOOL dotrans, GTVector<GTVector<Ftype>*> &du)
{
       GEOFLOW_TRACE();
  assert(du.size() >= GDIM
  && "Insufficient number of derivatives specified");


  GBOOL                         bembedded;
  GINT                          nxy;
  GSIZET                        k;
  GSIZET                        ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>              N(GDIM);
  GTVector<GTVector<Ftype>*>    W(GDIM);    // element weights
  GTVector<GTMatrix<Ftype>*>    Di(GDIM);   // element-based 1d derivative operators
  GElemList                    *gelems = &this->elems();

  bembedded = this->gtype() == GE_2DEMBEDDED;
  assert(( (bembedded && du.size()>=3)
  || (!bembedded&& du.size()>=GDIM) )
  && "Insufficient number of derviatves provided");

  nxy = bembedded ? GDIM+1 : GDIM;

#if defined(_G_IS2D)
  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( k=0; k<nxy ; k++ ) du[k]->range(ibeg, iend);
    for ( k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    for ( k=0; k<GDIM; k++ ) W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixW (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans);
    GMTK::Dg2_X_D1(*Di[0], *W [1], u, etmp, *du[0]); 
    GMTK::D2_X_Dg1(*W [0], *Di[1], u, etmp, *du[1]); 
    #if 0
    if ( bembedded ) { // ref 3-deriv is just W u:
      k = 0;
      for ( auto j=0; j<N[1]; j++ ) {
        for ( auto i=0; i<N[0]; i++ ) {
          (*du[2])[k] = u[k]*(*W[0])[i]*(*W[1])[j];  
          k++;
        }
      }
    }
    #endif
  }
  u.range_reset(); // reset global vec to global range
  for ( k=0; k<nxy; k++ ) du[k]->range_reset();

#elif defined(_G_IS3D)

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( k=0; k<nxy; k++ ) du[k]->range(ibeg, iend);
    for ( k=0; k<nxy; k++ ) {
      N[k]= (*gelems)[e]->size(k);
      W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    }
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixW (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(2)->getDerivMatrixW(!dotrans); 
    GMTK::Dg3_X_Dg2_X_D1(*Di[0], *W [1], *W [2], u, etmp, *du[0]); 
    GMTK::Dg3_X_D2_X_Dg1(*W [0], *Di[1], *W [2], u, etmp, *du[1]); 
    GMTK::D3_X_Dg2_X_Dg1(*W [0], *W [1], *Di[2], u, etmp, *du[2]); 
  }
  u.range_reset(); // reset global vec to globalrange
  for ( k=0; k<nxy; k++ ) du[k]->range_reset();

#endif

} // end of method compute_grefderivsW
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefdiv
// DESC   : Compute tensor product 'divergence' of input fields in ref space
//          for grid:
//             Div u = [I_X_I_X_Dx     |u1|
//                      I_X_Dy_X_I   . |u2|
//                      Dz_X_I_X_I]    |u3|
//          or
//     
//             Div u = [I_X_I_X_DxT)   |u1|
//                      I_X_DyT_X_I .  |u2|
//                      DzT_X_I_X_I]   |u3|
//          where Dx(T), Dy(T), Dz(T) are 1d derivative objects from basis functions
// ARGS   : 
//          u      : input vector field whose divergence we want, allocated globally 
//                   (e.g., for all elements). Must have GDIM components, unless
//                   we're using an embedded grid, when GDIM=2, when it should have
//                   3 components. Will not be altered on exit. If a component is
//                   NULL, we assume it's 0 here.
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          divu   : scalar result
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::compute_grefdiv(GTVector<GTVector<Ftype>*> &u, GTVector<Ftype> &etmp,
                            GBOOL dotrans, GTVector<Ftype> &divu)
{
       GEOFLOW_TRACE();
  GBOOL                        bembedded;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GTVector<Ftype>*>   W(GDIM);    // element 1/weights
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<Ftype>*>   Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &this->elems();

  bembedded = this->gtype() == GE_2DEMBEDDED;
#if 0
  assert(( (bembedded && u.size()==GDIM) 
  || (!bembedded&& u.size()==GDIM) )
  && "Insufficient number of vector field components provided");
#else
  assert(  u.size()==GDIM  
  && "Insufficient number of vector field components provided");
#endif

  divu = 0.0;

#if defined(_G_IS2D)

  for ( auto e=0; e<gelems->size(); e++ ) {
    // restrict global vecs to local range
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( auto k=0; k<u.size(); k++ ) if ( u[k]!=NULLPTR) u[k]->range(ibeg, iend); 
    for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    etmp.resizem((*gelems)[e]->nnodes());
    if ( u[0] != NULLPTR && u[0]->size() > 1 ) {
      GMTK::I2_X_D1(*Di[0], *u[0], N[0], N[1], etmp); // D1 u1
      divu += etmp;
    }
    if ( u[1] != NULLPTR && u[1]->size() > 1 ) {
      GMTK::D2_X_I1(*Di[1], *u[1], N[0], N[1], etmp); // D2 u2
      divu += etmp;
    }
#if 0
    if ( bembedded ) divu += *u[2]; // D3 acts as I
#endif
  }
  divu.range_reset();  // reset range to global scope
  for ( auto k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#elif defined(_G_IS3D)

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( auto k=0; k<u.size(); k++ ) if (u[k]!=NULLPTR) u[k]->range(ibeg, iend); 
    for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    etmp.resizem((*gelems)[e]->nnodes());
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 

    if ( u[0] != NULLPTR && u[0]->size() > 1 ) {
      GMTK::I3_X_I2_X_D1(*Di[0], *u[0], N[0], N[1], N[2], etmp); // D1 u1
      divu += etmp;
    }
    if ( u[1] != NULLPTR && u[1]->size() > 1 ) {
      GMTK::I3_X_D2_X_I1(*Di[1], *u[1], N[0], N[1], N[2], etmp); // D2 u2
      divu += etmp;
    }
    if ( u[2] != NULLPTR && u[2]->size() > 1 ) {
      GMTK::D3_X_I2_X_I1(*Di[2], *u[2], N[0], N[1], N[2], etmp); // D3 u3
      divu += etmp;
    }
  }
  divu.range_reset();  // reset range to global scope
  for ( auto k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#endif


} // end, method compute_grefdiv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivs
// DESC   : Compute tensor product derivs of specified field, u, in ref space
//          for grid, using grid object to determine which to compute. Compute:
//            du = [ I_X_I_X_Dx
//                   I_X_Dy_X_I
//                   Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need to vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::compute_grefderivs(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                               GBOOL dotrans, GTVector<GTVector<Ftype>*> &du)
{
	GEOFLOW_TRACE();
  assert(du.size() >= GDIM
       && "Insufficient number of derivatives specified");
  


  GBOOL                        bembedded;
  GINT                         nxy;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<Ftype>*>   Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &this->elems();

  bembedded = this->gtype() == GE_2DEMBEDDED;
  assert(( (bembedded && du.size()>=3) 
        || (!bembedded&& du.size()>=GDIM) )
       && "Insufficient number of derviatves provided");

  nxy = bembedded ? GDIM+1 : GDIM;

#if defined(_G_IS2D)

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( auto k=0; k<nxy ; k++ ) du[k]->range(ibeg, iend);
    for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::I2_X_D1(*Di[0], u, N[0], N[1], *du[0]); 
    GMTK::D2_X_I1(*Di[1], u, N[0], N[1], *du[1]); 
#if 0
    if ( bembedded ) { // ref 3-deriv is just W u:
      *du[2] = u;  
    }
#endif
  }
  u.range_reset(); // reset to global range
  for ( auto k=0; k<nxy; k++ ) du[k]->range_reset();

#elif defined(_G_IS3D)

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( auto k=0; k<GDIM; k++ ) du[k]->range(ibeg, iend);
    for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(2)->getDerivMatrix(!dotrans); 
    GMTK::I3_X_I2_X_D1(*Di[0], u, N[0], N[1], N[2], *du[0]); 
    GMTK::I3_X_D2_X_I1(*Di[1], u, N[0], N[1], N[2], *du[1]); 
    GMTK::D3_X_I2_X_I1(*Di[2], u, N[0], N[1], N[2], *du[2]); 
  }
  u.range_reset(); // reset global vec to globalrange
  for ( auto k=0; k<nxy; k++ ) du[k]->range_reset();

#endif

} // end of method compute_grefderivs


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderiv
// DESC   : Compute tensor product derivative in specified direction
//          of specified field, u, in ref space, using grid object.
//          Compute
//            du = [ I_X_I_X_Dx, or
//                   I_X_Dy_X_I, or
//                   Dz_X_I_X_I].
//     
//          depending on whether idir = 1, 2, or 3, respectively,
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          idir   : coordinate direction (1, 2,...,GDIM)
//          dotrans: flag telling us to take transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          du     : vector of length of u containing the derivative.
//          
//          Serves as main entry point to deriv computation
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::compute_grefderiv(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                              GINT idir, GBOOL dotrans, GTVector<Ftype> &du)
{
	GEOFLOW_TRACE();
  switch ( gderivtype_ ) {
    case GDV_VARP:
      grefderiv_varp(u, etmp, idir, dotrans, du);
      break;
    case GDV_CONSTP:
      grefderiv_constp (u, etmp, idir, dotrans, du);
      break;
    default:
      assert(false);
  }


} // end of method compute_grefderiv


//**********************************************************************************
//**********************************************************************************
// METHOD : grefderiv_varp
// DESC   : Compute tensor product derivative in specified direction
//          of specified field, u, in ref space, using grid object.
//          Compute
//            du = [ I_X_I_X_Dx, or
//                   I_X_Dy_X_I, or
//                   Dz_X_I_X_I].
//     
//          depending on whether idir = 1, 2, or 3, respectively,
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          idir   : coordinate direction (1, 2,...,GDIM)
//          dotrans: flag telling us to take transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          du     : vector of length of u containing the derivative.
//
//          Computes matrix by formaulating tensor products over individual
//          elements elements. May be used for case when order varies
//          amongh elements.
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::grefderiv_varp(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                             GINT idir, GBOOL dotrans, GTVector<Ftype> &du)
{
	GEOFLOW_TRACE();
  GSIZET               ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>     N(GDIM);
  GTMatrix<Ftype>     *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &this->elems();


#if defined(_G_IS2D)
  switch (idir) {
  case 1:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
      GMTK::I2_X_D1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 2:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
      GMTK::D2_X_I1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 3:
    assert( GDIM == 3
         && "Only GDIM reference derivatives");
    du = 0.0; //u;
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset to global range
  du.range_reset();

#elif defined(_G_IS3D)

  switch (idir) {
  case 1:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
      GMTK::I3_X_I2_X_D1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  case 2:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
      GMTK::I3_X_D2_X_I1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  case 3:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(2)->getDerivMatrix(!dotrans); 
      GMTK::D3_X_I2_X_I1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset global vec to globalrange
  du.range_reset();

#endif

} // end of method grefderiv_varp


//**********************************************************************************
//**********************************************************************************
// METHOD : grefderiv_constp
// DESC   : Compute tensor product derivative in specified direction
//          of specified field, u, in ref space, using grid object.
//          Compute
//            du = [ I_X_I_X_Dx, or
//                   I_X_Dy_X_I, or
//                   Dz_X_I_X_I].
//     
//          depending on whether idir = 1, 2, or 3, respectively,
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          idir   : coordinate direction (1, 2,...,GDIM)
//          dotrans: flag telling us to take transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          du     : vector of length of u containing the derivative.
//
//          Computes matrix by formaulating tensor products over all
//          elements elements. May be used for case when order doesn't vary
//          amongh elements.
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GGrid<Types>::grefderiv_constp(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                            GINT idir, GBOOL dotrans, GTVector<Ftype> &du)
{
	GEOFLOW_TRACE();
  GSIZET               ibeg, iend, Ne,  NN; // beg, end indices for global array
  GTVector<GSIZET>     N(GDIM);
  GTMatrix<Ftype>     *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &this->elems();


  for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[0]->size(k);
  Ne = gelems->size();


#if defined(_G_IS2D)
  switch (idir) {
  case 1:
    Di = (*gelems)[0]->gbasis(0)->getDerivMatrix (dotrans);
    GMTK::I2_X_D1<Ftype>(*Di, u, N[0], N[1], Ne, cudat_, du); 
    break;
  case 2:
    Di = (*gelems)[0]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::D2_X_I1<Ftype>(*Di, u, N[0], N[1], Ne, cudat_, du); 
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }

#elif defined(_G_IS3D)
  switch (idir) {
  case 1:
    Di = (*gelems)[0]->gbasis(0)->getDerivMatrix (dotrans); 
    NN = N[2]*N[1] * Ne;
    GMTK::I3_X_I2_X_D1<Ftype>(*Di, u, N[0], N[1], N[2], Ne, cudat_, du); 
    break;

  case 2:
    Di = (*gelems)[0]->gbasis(1)->getDerivMatrix(!dotrans); 
    GMTK::I3_X_D2_X_I1<Ftype>(*Di, u, N[0], N[1], N[2], Ne, cudat_, du); 
    break;

  case 3:
    Di = (*gelems)[0]->gbasis(2)->getDerivMatrix(!dotrans); 
    GMTK::D3_X_I2_X_I1<Ftype>(*Di, u, N[0], N[1], N[2], Ne, cudat_, du); 
    break;

  default:
    assert(FALSE && "Invalid coordinate direction");
  }

#endif

} // end of method grefderiv_constp


//**********************************************************************************
//**********************************************************************************
// METHOD : ispconst
// DESC   : Check if p is const over all elements
// ARGS   : none. 
// RETURNS: TRUE if p is const; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GGrid<Types>::ispconst()
{
	GEOFLOW_TRACE();
  GBOOL       bconst;
  GSIZET      Ne, N0[GDIM];
  GElemList  *gelems = &this->elems();

  Ne = gelems->size();
  for ( auto k=0; k<GDIM; k++ ) N0[k]= (*gelems)[0]->size(k);

  bconst = TRUE;
  for ( auto i=1; i<Ne && bconst; i++ ) {
    for ( auto k=0; k<GDIM; k++ ) {
      bconst = bconst && N0[k] == (*gelems)[i]->size(k);
    }
  }

  return bconst;

} // end of method ispconst


//**********************************************************************************
//**********************************************************************************
// METHOD : set_derivtype
// DESC   : Set derivative type. Does some checking to ensure
//          that it's valid.
// ARGS   : GDerivType flag
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::set_derivtype(GDerivType gt)
{
	GEOFLOW_TRACE();
  if ( !bpconst_ ) {
    assert( gt == GDV_VARP );
  }

  gderivtype_ = gt;

} // end of method set_derivtype




//**********************************************************************************
//**********************************************************************************
// METHOD : set_elemids
// DESC   : Set element id vector from individiual elements
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::set_elemids()
{
	GEOFLOW_TRACE();

  gelemids_.resize(gelems_.size());
  for ( auto e=0; e<gelems_.size(); e++ ) {
    gelemids_[e] = gelems_[e]->elemid();
  }

} // end of method set_elemids




//**********************************************************************************
//**********************************************************************************
// METHOD : init_qdealias
// DESC   : Do initialization for quadratic dealiasing 
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::init_qdealias()
{
	GEOFLOW_TRACE();
  if ( !bpconst_ ) {
    cout << "GGrid<Types>::dealias: dealiasing requires constant p" << endl;
    assert( FALSE ); 
  }

  GSIZET                       n;
  GSIZET                       nelems = gelems_.size();
  GTVector<GINT>               N(GDIM);
  GTVector<Ftype>              qxi; 
  GTVector<GTVector<Ftype>*>   W(GDIM), qW1d(GDIM);    // element weights
  GTVector<GLLBasis<GCTYPE,Ftype>>
                               dbasis(GDIM);


  // Compute dealias bases and interp matrices:
  IQPdealias_ .resize(GDIM);
  IQPdealiasT_.resize(GDIM);
  for ( auto k=0; k<GDIM; k++ ) {
    W[k] = gelems_[0]->gbasis(k)->getWeights();
    N[k] = gelems_[0]->size(k);
    assert( pqdealias_[k] >= N[k]-1 );
    dbasis[k].resize(pqdealias_[k]);
    dbasis[k].getXiNodes(qxi);
    qW1d        [k] = dbasis[k].getWeights();
    IQPdealias_ [k].resize(pqdealias_[k]+1, N[k]);
    IQPdealiasT_[k].resize(N[k],pqdealias_[k]+1);
    gelems_[0]->gbasis(k)->evalBasis(qxi,IQPdealias_[k]);
    IQPdealias_[k].transpose(IQPdealiasT_[k]);
  }


  qdN_.resize(1);
  // Compute dealias tensor prod weight matrix (diagnonal:
#if defined(_G_IS2D)
  qdN_[0] = (pqdealias_[0]+1)*(pqdealias_[1]+1); // # nodes 
  qW_.resize(nelems*qdN_[0]);
  for ( auto e=0, n=0; e<gelems_.size(); e++ ) {
    for ( auto j=0; j<pqdealias_[1]+1; j++ ) {
      for ( auto i=0; i<pqdealias_[0]+1; i++ ) {
        qW_[n] = (*qW1d[0])[i] * (*qW1d[1])[j];
        n++;
      }
    }
  }
#elif defined(_G_IS3D)
  qdN_[0] = (pqdealias_[0]+1)*(pqdealias_[1]+1)*(pqdealias_[2]+1); // # nodes
  qW_.resize(nelems*qdN_[0]);
  for ( auto e=0, n=0; e<gelems_.size(); e++ ) {
    for ( auto k=0; k<pqdealias_[2]+1; k++ ) {
      for ( auto j=0; j<pqdealias_[1]+1; j++ ) {
        for ( auto i=0; i<pqdealias_[0]+1; i++ ) {
          qW_[n] = (*qW1d[0])[i] * (*qW1d[1])[j] * (*qW1d[2])[k];
          n++;
        }
      }
    }
  }
#endif

  // Compute inverse tensor prod weights (not mass matrix)
  // for original explansion basis:

  StateComp *pjac   = &this->Jac();
  StateComp *pimass =  this->imassop().data();
  iWp_.resize(pimass->size());
  for ( auto j=0; j<pimass->size(); j++ ) {
    iWp_[j] = (*pjac)[j] * (*pimass)[j]; // 1/p-Weights: divide out Jac
  }

  // Allocate quadratic dealising tmp space:
  qdtmp_.resize(2);
  for ( auto j=0; j<qdtmp_.size(); j++ ) qdtmp_[j].resize(qW_.size());

  bInitQDealias_ = TRUE;

} // end of method init_qdealias


//**********************************************************************************
//**********************************************************************************
// METHOD : dealias
// DESC   : Do dealiasing for quadratic nonlinearity, v1 * v2
// ARGS   : v1      : first variable in product
//          v2      : second vairable in product
//          prod    : dealiased product of v1 * v2
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGrid<Types>::dealias(StateComp &v1, StateComp &v2, StateComp &prod)
{
	GEOFLOW_TRACE();

  if ( !doQDealias_ ) {  // multiply w/o dealiasing
    prod = v1; prod.pointProd(v2);
    return;
  }

  if ( !bInitQDealias_ ) init_qdealias();

  if ( !bpconst_ ) {
    cout << "GGrid<Types>::dealias: dealiasing requires constant p" << endl;
    assert( FALSE ); 
  }

  GSIZET ibeg, iend, iqbeg, iqend;
  
  // First, compute IQPdealias * (v1, v2):
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    iqbeg = e*qdN_[0]; iqend = iqbeg + qdN_[0] - 1;
    // Restrict global data to local scope:
    v1       .range  (ibeg,iend); v2       .range  (ibeg,iend);
    qdtmp_[0].range(iqbeg,iqend); qdtmp_[1].range(iqbeg,iqend);
#if defined(_G_IS2D)
    GMTK::D2_X_D1<Ftype>(IQPdealias_[0], IQPdealiasT_[1], v1, tptmp_, qdtmp_[0]);
    GMTK::D2_X_D1<Ftype>(IQPdealias_[0], IQPdealiasT_[1], v2, tptmp_, qdtmp_[1]);
#elif defined(_G_IS3D)
    GMTK::D3_X_D2_X_D1<Ftype>(IQPdealias_[0], IQPdealiasT_[1], IQPdealiasT_[2], v1, tptmp_, qdtmp_[0]);
    GMTK::D3_X_D2_X_D1<Ftype>(IQPdealias_[0], IQPdealiasT_[1], IQPdealiasT_[2], v2, tptmp_, qdtmp_[1]);
#endif
  } // end, element loop
  v1.range_reset(); v2.range_reset();
  qdtmp_[0].range_reset(); qdtmp_[1].range_reset();

  // Compute pointProd of over-sampled terms:
  qdtmp_[0].pointProd(qdtmp_[1]);

  // Do a "Galerkin projection" back to original space:
  //    p_prod =  pW^-1 IQPdealiasT qW  qProd:
  qdtmp_[0].pointProd(qW_); // qW * qProd

  // Do tensor product application of IQPdealiasT:
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    iqbeg = e*qdN_[0]; iqend = iqbeg + qdN_[0] - 1;
    // Restrict global data to local scope:
    prod.range       (ibeg,iend); 
    qdtmp_[0].range(iqbeg,iqend); 
#if defined(_G_IS2D)
    GMTK::D2_X_D1<Ftype>(IQPdealiasT_[0], IQPdealias_[1], qdtmp_[0], tptmp_, prod);
#elif defined(_G_IS3D)
    GMTK::D3_X_D2_X_D1<Ftype>(IQPdealiasT_[0], IQPdealias_[1], IQPdealias_[2], qdtmp_[0], tptmp_, prod);
#endif
  } // end, element loop
  prod.range_reset();
  qdtmp_[0].range_reset();
  
  prod.pointProd(iWp_); // apply pW^-1

} // end of method dealias


