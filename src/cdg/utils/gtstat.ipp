//==================================================================================
// Module       : gtstat.ipp
// Date         : 10/1/19 (DLR)
// Description  : Encapsulates the methods and data associated with
//                GeoFLOW statistics methods.
// Copyright    : Copyright 2019, Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantitate with no. bins
// ARGS   : GSIZET n: no. bins
// RETURNS: none
//**********************************************************************************
template<typename T>
GTStat<T>::GTStat(GSIZET n, GC_COMM comm):
bfixedwidth_   (FALSE),
nbins_         (n),
nkeep_         (0),
fixedwidth_    (0.0),
comm_       (comm)
{
  assert(nbins_ > 0 && "Invalid bin count");

  myrank_  = GComm::WorldRank(comm_);

  lpdf_ .resize(nbins_);

} // end, constuctor (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Instantitate with constant bin widsth 
// ARGS   : T w: bin width
// RETURNS: none
//**********************************************************************************
template<typename T>
GTStat<T>::GTStat(T w, GC_COMM comm):
bfixedwidth_   (TRUE),
nbins_         (0),
nkeep_         (0),
fixedwidth_    (w),
comm_       (comm)
{
  assert(nbins_ > 0 && "Invalid bin count");

  myrank_  = GComm::WorldRank(comm_);


} // end, constuctor (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method
// DESC   : Override default copy constructor
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<typename T>
GTStat<T>::GTStat(const GTStat<T> &obj):
nbins_      (obj.get_nbins())
{
  lpdf_ .resize(nbins_);
}


//**********************************************************************************
//**********************************************************************************
// METHOD : dopdf1d (1)
// DESC   : Do 1d pdf of input scalar field, and output to specified file.
// ARGS   : 
//          u      : scalar field on which to operate
//          ifixmin: if TRUE, then use specified fmin, to set lower dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          ifixmax: if TRUE, then use specified fmax, to set lower dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf, used for binning if ifixmin/max=TRUE,
//                   returned with values found in u if ifixmin/max=FALSE
//          iside  : if 1, considers only u>0 data; if -1, considers only
//                   u<0 data; if 0, considers all data. 
//          dolog  : take log of |u| when creating bins?
//          utmp   : tmp vector of same size as u
//          pdf    : final pdf
// RETURNS: none.
//**********************************************************************************
template<typename T> 
void GTStat<T>::dopdf1d(GTVector<T> &u, GBOOL ifixmin, GBOOL ifixmax, T &fmin, T &fmax, GINT iside, GBOOL dolog, GTVector<T> &utmp, GTVector<T> &pdf)
{
  GSIZET ibin, iend, j, lkeep;
  T      bmin, bmax, del, test;
  T      fbin, ftmp;
  T      sumr, xnorm;
  T      tiny;


  tiny = fabs(100.0*std::numeric_limits<T>::epsilon());

  utmp = u;
  iend = utmp.size()-1;

  // Check if doing 1-sided pdf:
  if      ( iside ==  1 ) { // keep positive entries
    for ( auto j=0; j< utmp.size(); j++ ) {
      utmp[j] = MAX(0.0, utmp[j]);
    }   
    utmp.sortdecreasing();
    utmp.contains(0.0, iend); iend--;
  }
  else if ( iside == -1 ) { // keep negative entries
    for ( auto j=0; j< utmp.size(); j++ ) {
      utmp[j] = MIN(0.0, utmp[j]);
    }   
    utmp.sortincreasing();
    utmp.contains(0.0, iend); iend--;
  }

  if ( dolog ) {
    utmp.abs();
  } 

  // Compute bin dynamic range, if not using specified range:
  if ( !ifixmin ) {
    utmp.range(0,iend);
    ftmp = utmp.min();
    GComm::Allreduce(&ftmp, &fmin, 1, T2GCDatatype<T>() , GC_OP_MIN, comm_);
    if ( dolog ) {
      fmin += tiny;
    }
    utmp.range_reset();
  }

  if ( !ifixmax ) {
    utmp.range(0,iend);
    ftmp = utmp.max();
    GComm::Allreduce(&ftmp, &fmax, 1, T2GCDatatype<T>() , GC_OP_MAX, comm_);
    if ( dolog ) {
      fmax += tiny;
    }
    utmp.range_reset();
  }

  bmin = fmin; // set bin min/max
  bmax = fmax;
  if ( dolog ) {
    bmin = log10(fabs(fmin)); 
    bmax = log10(fabs(fmax));
  }

  if ( bfixedwidth_ ) {
    if ( dolog ) {
      nbins_ = static_cast<GSIZET>( (bmax - bmin) / log10(fixedwidth_) );
    }
    else {
      nbins_ = static_cast<GSIZET>( (bmax - bmin) / fixedwidth_ );
    }
  }

  gpdf_.resize(nbins_);
  ikeep_.resize(u.size());
  lpdf_.resizem(nbins_);
  pdf  .resizem(nbins_);
  lpdf_ = 0.0;
  pdf   = 0.0;

  // Find indirection indices that meet dyn. range criterion:
  lkeep = 0;
  #pragma omp for
  for ( j=0; j<utmp.size(); j++ ) {
    if ( utmp[j] >= fmin && utmp[j] <= fmax ) {
      ikeep_[lkeep++] = j;
    }
  }
  GComm::Allreduce(&lkeep, &nkeep_, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);
  assert(nkeep_ > 0  && "No samples within dynamic range");

  xnorm = 1.0 / static_cast<T>(nkeep_);

  // Compute average:
  sumr = 0.0;
  #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
  for ( j=0; j<lkeep; j++ ) {
    sumr += utmp[ikeep_[j]];
  }
  GComm::Allreduce(&sumr, &gavg_, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  gavg_ *= xnorm;


  // Compute std deviation:
  sumr = 0.0;
  #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
  for ( j=0; j<lkeep; j++ ) {
    sumr += (utmp[ikeep_[j]]-gavg_)*(utmp[ikeep_[j]]-gavg_);
  }
  GComm::Allreduce(&sumr, &sig_, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  sig_ = sqrt(sig_*xnorm);

  // Note: We _may_ want to compute higher order quantities like
  //       skewness, flatness. If so, do that here.

  // Compute local PDF:
  del = fabs(bmax - bmin) / nbins_;
  if ( dolog ) {
    #pragma omp parallel for  private(ibin,test)
    for ( j=0; j<lkeep; j++ ) {
      test = log10(utmp[ikeep_[j]]+tiny);
      ibin = static_cast<GSIZET> ( ( test - bmin )/del );
      ibin = MIN(MAX(ibin,0),nbins_-1);
      #pragma omp atomic
      lpdf_[ibin] += 1.0;
    }
  }
  else {
    #pragma omp parallel for  private(ibin,test)
    for ( j=0; j<lkeep; j++ ) {
      test = utmp[ikeep_[j]];
      ibin = static_cast<GSIZET> ( ( test - bmin )/del );
      ibin = MIN(MAX(ibin,0),nbins_-1);
      #pragma omp atomic
      lpdf_[ibin] += 1.0;
    }
  }
  
  // Compute global reduction between MPI tasks to find final (global) pdf:
  GComm::Allreduce(lpdf_.data(), pdf.data(), nbins_, T2GCDatatype<T>() , GC_OP_SUM, comm_);


  // Do sanity check:
  fbin = 0.0;
  #pragma omp parallel for  default(shared) reduction(+:fbin)
  for ( j=0; j<nbins_; j++ ) {
    #pragma omp atomic
    fbin += pdf[j];
  }
  assert( fbin == nkeep_ && "Inconsistent binning");


} // end, dopdf1d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : dopdf1d (2)
// DESC   : Do 1d pdf of input scalar field, and return pdf. Additional member
//          is allocated in this call (gpdf_).
// ARGS   : 
//          u      : scalar field on which to operate
//          ifixmin: if TRUE, then use specified fmin, to set lower dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          ifixmax: if TRUE, then use specified fmax, to set lower dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf bins
//          iside  : if 1, considers only u>0 data; if -1, considers only
//                   u<0 data; if 0, considers all data. 
//          dolog  : take log of |u| when creating bins?
//          utmp   : tmp vector of same size as u
//          fname  : filename to which to output pdf
// RETURNS: none.
//**********************************************************************************
template<typename T> 
void GTStat<T>::dopdf1d(GTVector<T> &u, GBOOL ifixmin, GBOOL ifixmax, T &fmin, T &fmax, GINT iside, GBOOL dolog, GTVector<T> &utmp, const GString &fname)
{
  GSIZET            j;
  std::ofstream     ios;
  std::stringstream header;


  dopdf1d(u, ifixmin, ifixmax, fmin, fmax, iside, dolog, utmp, gpdf_);
  if ( myrank_ == 0 )  {

    header << std::scientific << std::setprecision(8);
    header << "# range=[" 
           << fmin   << ","
           << fmax   << "]; avg="
           << gavg_  << "; sig="
           << sig_   << "; nbins="
           << nbins_ << "; blog ="
           << dolog  << "; nkeep="
           << nkeep_ ;

     ios.open(fname,std::ios_base::trunc);
     ios << std::scientific << std::setprecision(15);
     ios << header.str() << std::endl;
     for ( j=0; j<nbins_-1; j++ ) {
       ios << gpdf_[j] << " " << std::endl;
     }
     ios << gpdf_[nbins_-1] << std::endl;
    
     ios.close();
  }


} // end, dopdf1d (2)


