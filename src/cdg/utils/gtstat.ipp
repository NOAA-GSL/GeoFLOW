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
// METHOD : Constructor method
// DESC   : Instantitate with no. bins
// ARGS   : GSIZET n: no. bins
// RETURNS: none
//**********************************************************************************
template<typename T>
GTStat<T>::GTStat(GSIZET n, GC_COMM comm):
nbins_         (n),
nkeep_         (0),
comm_       (comm)
{
  assert(nbins_ > 0 && "Invalid bin count");

  myrank_  = GComm::WorldRank(comm_);

  lpdf_ .resize(nbins_);
}


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
//          ifixdr : if TRUE, then use specified fmin, fmax to set dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf
//          dolog  : take log of |u| when creating bins?
//          pdf    : final pdf
// RETURNS: none.
//**********************************************************************************
template<typename T> 
void GTStat<T>::dopdf1d(GTVector<T> &u, GBOOL ifixdr, T &fmin, T &fmax, GBOOL dolog, GTVector<T> &pdf)
{
  GSIZET ibin, j, lkeep;
  T      del, test, tmax, tmin;
  T      fbin, fmin1, fmax1;
  T      sumr, xnorm;
  T      tiny;

  gpdf_.resize(nbins_);
  ikeep_.resize(u.size());
  lpdf_ = 0.0;
  pdf   = 0.0;

  tiny = fabs(100.0*std::numeric_limits<T>::epsilon());

  // Compute dynamic range, if not using specified range:
  if ( !ifixdr ) {
    fmin1 = u.min();
    fmax1 = u.max();
    GComm::Allreduce(&fmin1, &fmin, 1, T2GCDatatype<T>() , GC_OP_MIN, comm_);
    GComm::Allreduce(&fmax1, &fmax, 1, T2GCDatatype<T>() , GC_OP_MAX, comm_);
  }
  
  if ( dolog ) {
    fmin = log10(fabs(fmin)+tiny);
    fmax = log10(fabs(fmax)+tiny);
    if ( fmin > fmax ) {
      fmin1 = fmin;
      fmin = fmax;
      fmax = fmin1;
    }
  }

  tmin = fmin;
  tmax = fmax;
  if ( dolog ) {
    tmin = pow(10.0, fmin);
    tmax = pow(10.0, fmax);
  }


  // Find indirection indices that meet dyn. range criterion:
  lkeep = 0;
  if ( dolog ) {
  #pragma omp for
    for ( auto j=0; j<u.size(); j++ ) {
      if ( fabs(u[j]) >= tmin && fabs(u[j]) <= tmax ) {
        ikeep_[lkeep++] = j;
      }
    }
  }
  else {
  #pragma omp for
    for ( auto j=0; j<u.size(); j++ ) {
      if ( u[j] >= tmin && u[j] <= tmax ) {
        ikeep_[lkeep] = j;
        lkeep++;
      }
    }
  }
  GComm::Allreduce(&lkeep, &nkeep_, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);
  assert(nkeep_ > 0  && "No samples are within dynamic range");

  xnorm = 1.0 / static_cast<T>(nkeep_);

  // Compute average:
  sumr = 0.0;
  if ( dolog ) {
    #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
    for ( j=0; j<lkeep; j++ ) {
      sumr += fabs(u[ikeep_[j]]);
    }
  }
  else {
    #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
    for ( j=0; j<lkeep; j++ ) {
      sumr += u[ikeep_[j]];
    }
  }
  GComm::Allreduce(&sumr, &gavg_, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  gavg_ *= xnorm;


  // Compute std deviation:
  sumr = 0.0;
  #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
  for ( j=0; j<lkeep; j++ ) {
    sumr += pow(u[ikeep_[j]]-gavg_,2);
  }
  GComm::Allreduce(&sumr, &sig_, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  sig_ = sqrt(sig_/xnorm);

  // Compute local PDF:
  del = fabs(fmax - fmin) / nbins_;
  if ( dolog ) {
    #pragma omp parallel for  private(ibin,test)
    for ( j=0; j<lkeep; j++ ) {
      test = log10(fabs(u[ikeep_[j]])+tiny);
      ibin = static_cast<GSIZET> ( ( test - fmin )/del );
      ibin = MIN(MAX(ibin,0),nbins_-1);
      #pragma omp atomic
      lpdf_[ibin] += 1.0;
    }
  }
  else {
    #pragma omp parallel for  private(ibin,test)
    for ( j=0; j<lkeep; j++ ) {
      test = u[ikeep_[j]];
      ibin = static_cast<GSIZET> ( ( test - fmin )/del );
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
//          ifixdr : if TRUE, then use specified fmin, fmax to set dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf
//          dolog  : take log of |u| when creating bins?
//          fname  : filename to which to output pdf
// RETURNS: none.
//**********************************************************************************
template<typename T> 
void GTStat<T>::dopdf1d(GTVector<T> &u, GBOOL ifixdr, T &fmin, T &fmax, GBOOL dolog, const GString &fname)
{
  GSIZET            j;
  T                 ofmax, ofmin;
  std::ofstream     ios;
  std::stringstream header;

  gpdf_.resize(nbins_);

  dopdf1d(u, ifixdr, fmin, fmax, dolog, gpdf_);
  if ( myrank_ == 0 )  {

    ofmin = fmin;
    ofmax = fmax;
    if ( dolog ) {
      ofmin = pow(10.0, fmin);
      ofmax = pow(10.0, fmax);
    }

    header << std::scientific << std::setprecision(8);
    header << "# range=[" 
           << ofmin  << ","
           << ofmax  << "]; avg="
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


