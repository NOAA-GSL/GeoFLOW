//==================================================================================
// Module       : ggfx.ipp
// Date         : 5/9/18 (DLR)
// Description  : Template functions implementations for GGFX
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <math.h>
#include <type_traits>
#include <assert.h>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"

using namespace std;


//************************************************************************************
//************************************************************************************
// METHOD : doOp (1)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Must
//                 be same size as glob_index in Init call.
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T>  
GBOOL GGFX<T>::doOp(GTVector<T> &u, GGFX_OP op) 
{
#if 0
  GSIZET   nu = u.size();
  T       *uu = u.data();
  return doOp(uu, nu,  op);
#else

  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GINT      irank;
  GINT      i, j;
  GBOOL     bret=TRUE;
  GString   serr = "GGFX<T>::doOp(1): ";

  GPP(comm_, serr << "iOpL2LIndices=" << iOpL2LIndices_);
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "Doing first localGS...");
#endif
  
  // For each global index row in iOpL2LIndices, gather and
  // scatter data locally: do operation for each of the
  // shared nodes; 
  bret = localGS(u, iOpL2LIndices_, nOpL2LIndices_, op);
  if ( !bret ) {
    std::cout << serr << "localGS (1) failed " << std::endl;
    exit(1);
  }
  
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "First localGS done.");
#endif

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "iL2RIndices=" << iOpL2RIndices_);
#endif

  // If no other tasks, return:
  if ( nprocs_ == 1 ) bret = TRUE;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Doing dataExchange...");
#endif

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Doing dataExchange");
#endif
  // Perform a exchange of field data:
  bret = dataExchange(u);
  if ( !bret ) {
    std::cout << serr << "dataExchange.failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "dataExchange done.");
#endif
  GPP(comm_,serr << " recvBuff=" << recvBuff_);
  GPP(comm_,serr << " iOpR2LIndices=" << iOpR2LIndices_);
  GPP(comm_,serr << " nOpR2LIndices=" << nOpR2LIndices_);
//GPP(comm_,serr << " recvBuff.size=" << recvBuff_.data().size());
  bret = localGS(u, iOpL2RIndices_, nOpL2RIndices_, op, &recvBuff_);
  if ( !bret ) {
    std::cout << serr << "localGS (2) failed " << std::endl;
    exit(1);
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "done.");
#endif
  return bret;
#endif

} // end, method doOp (1)


//************************************************************************************
//************************************************************************************
// METHOD : doOp (2)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Can be 
//                 _any_ size, but the number of elems use din iind array 
//                 _must_ be the same size as the number in glob_index in Init call. 
//          iind : indirection array into u. Must contain indices that
//                 corresp to those in glob_index in Init call.
//          nind : number of indirection indices in iind; must be the same
//                 as the number of glob_index in Init call.
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T> GBOOL GGFX<T>::doOp(GTVector<T> &u, GTVector<GSIZET> &iind, GGFX_OP op)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GINT      irank;
  GINT      i, j;
  GBOOL     bret=TRUE;
  GString   serr = "GGFX<T>::doOp(2): ";

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing first localGS..." << std::endl;
#endif

  // For each global index row in iOpL2LIndices, gather and
  // scatter data locally: get one operation for each of the
  // shared nodes; 
  bret = localGS(u, iind, iOpL2LIndices_, nOpL2LIndices_, op);
  if ( !bret ) {
    std::cout << serr << "localGS (1) failed " << std::endl;
    exit(1);
  }
  
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "First localGS done ." << std::endl;
#endif

  // If no other tasks, return:
  if ( nprocs_ == 1 ) bret = TRUE;


#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing dataExchange..." << std::endl;
#endif

  // Perform a exchange of field data:
  bret = dataExchange(u, iind);
  if ( !bret ) {
    std::cout << serr << "dataExchange failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "dataExchange done ." << std::endl;
#endif


#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing second localGS..." << std::endl;
#endif
  bret = localGS(u, iind, iOpL2RIndices_, nOpL2RIndices_, op, &recvBuff_);
  if ( !bret ) {
    std::cout << serr << "localGS (2) failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Second localGS done." << std::endl;
#endif

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing final assignment..." << std::endl;
#endif

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "done." << std::endl;
#endif
  return bret;

} // end of method doOp(2)


//************************************************************************************
//*************************************************************************************
// METHOD : localGS (1)
// DESC   : Performs GGFX_OP (op) on the array of field values, u
//          and scatters them to the same dof
// ARGS   : u      : input field whose elements are to be combined. May
//                   be larger than number of elements in glob_index array
//                   used in call to Init, as long as indirection indices,
//                   iind, are being used
//          ilocal : Matrix of indirection indices into qu array. Each column represents
//                   a shared node, and each row represents the other local indices 
//                   that point to the same shared node, so that the values at the same
//                   global index can be summed, multiplied, max'd or min'd.
//          nlocal : for each shared node (column in ilocal), local indices to it
//          op     : operation to be performed in the combination
//          qop    : operand; used if non_NULL; must have same no elements as columns in
//                   ilocal. qop==NULLPTR is default. When called, will likely be the 
//                   receive buffer.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX<T>::localGS(GTVector<T> &qu, GSZMatrix &ilocal, GSZBuffer &nlocal, GGFX_OP op, GTMatrix<T> *qop)
{

  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString serr = "GGFX<T>::localGS (1): ";
  T       res;
  GLLONG  i, j;

  GPP(comm_,serr << "op=" << op);

  // Perform GGFX_OP on the nodes shared by this proc:
  switch(op) {
    case GGFX_OP_SUM:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(0);
          for ( i=0; i<nlocal[j]; i++ ) { // do gather
             res += qu[ilocal(i,j)];
          }
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] += (*qop)(i,j);
          }
        }
      }
      break;
    case GGFX_OP_PROD:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(1);
          for ( i=0; i<nlocal[j]; i++ ) {
             res *= qu[ilocal(i,j)];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] *= (*qop)(i,j);
          }
        }
      }
      break;
    case GGFX_OP_MAX:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::min();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MAX(res,qu[ilocal(i,j)]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = MIN(qu[ilocal(i,j)],(*qop)(i,j));
          }
        }
      }
      break;
    case GGFX_OP_MIN:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MIN(res,qu[ilocal(i,j)]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = MIN(qu[ilocal(i,j)],(*qop)(i,j));
          }
        }
      }
      break;
    case GGFX_OP_SMOOTH:
      if ( qop == NULLPTR ) {
GPP(comm_,serr<<"imult="<<imult_);
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(0);
          for ( i=0; i<nlocal[j]; i++ ) { // do gather
             res += qu[ilocal(i,j)];
          }
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] += (*qop)(i,j);
          }
        }
        qu.pointProd(imult_); // apply inverse multiplicity to avg:
      }
      break;
    default:
      return FALSE;
  }

  return TRUE;

} // end of method localGS (1)


//************************************************************************************
//*************************************************************************************
// METHOD : localGS (2)
// DESC   : Performs GGFX_OP (op) on the array of field values, u,
//          and scatters them to the same dof
// ARGS   : u      : input field whose elements are to be combined. May
//                   be larger than number of elements in glob_index array
//                   used in call to Init, as long as indirection indices,
//                   iind, are being used
//          iind   : indirection array into u. Must contain indices that
//                   corresp to those in glob_index in Init call.
//          ilocal : Matrix of indirection indices into qu array. Each row represents
//                   a shared node, and each column represents the other local indices 
//                   that point to the same shared node, so that the values at the same
//                   global index can be summed, multiplied, max'd or min'd.
//          nlocal : for each shared node (row in ilocal), local indices to it
//          op     : operation to be performed in the combination
//          qop    : operand; used if non_NULL; must have same no elements as columns in
//                   ilocal.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX<T>::localGS(GTVector<T> &qu, GTVector<GSIZET> &iind, 
                    GSZMatrix &ilocal, GSZBuffer &nlocal,  GGFX_OP op, GTMatrix<T> *qop)
{

  assert( std::is_arithmetic<T>::value && "Illegal template type");


  GString serr = "GGFX<T>::localGS (2): ";
  T       res;
  GLLONG  i, j;


  // Perform GGFX_OP on the nodes shared by this proc:
  switch(op) {
    case GGFX_OP_SUM:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(0);
          for ( i=0; i<nlocal[j]; i++ ) {
             res += qu[iind[ilocal(i,j)]];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] += (*qop)(i,j);
          }
        }
      }
      break;
    case GGFX_OP_PROD:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(1);
          for ( i=0; i<nlocal[j]; i++ ) {
             res *= qu[iind[ilocal(i,j)]];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] *= (*qop)(i,j);
          }
        }
      }
      break;
    case GGFX_OP_MAX:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::min();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MAX(res,qu[iind[ilocal(i,j)]]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = MAX(qu[iind[ilocal(i,j)]],(*qop)(i,j));
          }
        }
      }
      break;
    case GGFX_OP_MIN:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MIN(res,qu[iind[ilocal(i,j)]]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = MAX(qu[iind[ilocal(i,j)]],(*qop)(i,j));
          }
        }
      }
      break;
    case GGFX_OP_SMOOTH:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = static_cast<T>(0);
          for ( i=0; i<nlocal[j]; i++ ) { // do gather
             res += qu[iind[ilocal(i,j)]];
          }
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[iind[ilocal(i,j)]] += (*qop)(i,j);
          }
        }
        for ( i=0; i<iind.size(); i++ ) { // do H1-smoothing
          qu[iind[i]] *= imult_[iind[i]];
        }
      }
      break;
    default:
      return FALSE;
  }

  return TRUE;

} // end of method localGS (2)


//************************************************************************************
//************************************************************************************
// METHOD : dataExchange (1)
// DESC   : Perform data exchange between procs by using asynchronous
//              recvs...
// ARGS   :
//          u    : field data 
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX<T>::dataExchange(GTVector<T> &u)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString       serr = "GGFX<T>::dataExchange(1): ";
  GCommDatatype dtype=T2GCDatatype<T>();
  GINT          i, j;
  GBOOL         bret=TRUE;

  sendBuff_.set(-1);
  recvBuff_.set(-1);

  // Fill send buffer with data:
  for ( j=0; j<iOpL2RTasks_.size() && bret; j++ ) {
 // if ( iOpL2RTasks_[j] == rank_ ) continue;
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) {
      sendBuff_(i,j) = u[iOpL2RIndices_(i,j)];
    }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GIBuffer rtasks(iOpL2RTasks_);
  GIBuffer isend(iOpL2RTasks_.size());
  GPP(comm_,serr << "iOpL2RTasks=" << iOpL2RTasks_);
  GPP(comm_,serr << "sendBuff=" << sendBuff_);
#endif

  bret = GComm::ASendRecv(recvBuff_.data().data(),iOpL2RTasks_.size(),NULLPTR,recvBuff_.size(1),dtype,iOpL2RTasks_.data(), TRUE, 
                          sendBuff_.data().data(),iOpL2RTasks_.size(),NULLPTR,sendBuff_.size(1),dtype,iOpL2RTasks_.data(), comm_);

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "(2)sendBuff=" << sendBuff_);
  GPP(comm_,serr << "(2)recvBuff=" << recvBuff_);
#endif

  return bret;
} // end of dataExchange (1)


//************************************************************************************
//************************************************************************************
// METHOD : dataExchange (2)
// DESC   : Perform data exchange between procs by using asynchronous
//              recvs...
// ARGS   :
//          u    : field data 
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX<T>::dataExchange(GTVector<T> &u, GTVector<GSIZET> &iind)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString       serr = "GGFX<T>::dataExchange(2): ";
  GCommDatatype dtype=T2GCDatatype<T>();
  GINT          i, j;
  GBOOL         bret=TRUE;


  // Fill send buffer with data:
  for ( j=0; j<iOpL2RTasks_.size() && bret; j++ ) {
//  if ( iOpL2RTasks_[j] == rank_ ) continue;
    for ( i=0; j<nOpL2RIndices_[j]; i++ ) {
      sendBuff_(i,j) = u[iind[iOpL2RIndices_(i,j)]];
    }
  }

recvBuff_.set(-1);

  bret = GComm::ASendRecv(recvBuff_.data().data(),(GINT)recvBuff_.size(2),NULLPTR,(GINT)recvBuff_.size(1),dtype,iOpL2RTasks_.data(), TRUE, 
                          sendBuff_.data().data(),(GINT)sendBuff_.size(2),NULLPTR,(GINT)sendBuff_.size(1),dtype,iOpL2RTasks_.data(), comm_);


  return bret;
} // end of dataExchange (2)

//************************************************************************************
//************************************************************************************
// METHOD     : initMult
// DESCRIPTION: Initialize data for H1-smoothing
// ARGUMENTS  : 
// RETURNS    : 
//************************************************************************************
template<typename T> 
void GGFX<T>::initMult()
{

  assert(bInit_ && "Operator not initialized");

  GDVector mult(nglob_index_);

  mult = 1.0;

  // Do DSS sum to find multiplicity:
  doOp(mult, GGFX_OP_SUM);


  // Compute 1/mult:
  imult_.resize(mult.size());
  for ( GSIZET j=0; j<mult.size(); j++ ) {
    imult_[j] = 1.0/mult[j];
  }

  GPP(comm_,"GGFX<T>::initMult: mult=" << mult);
  GPP(comm_,"GGFX<T>::initMult: imult=" << imult_);

} // end of method initMult


