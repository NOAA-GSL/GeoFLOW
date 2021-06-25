//==================================================================================
// Module       : gcomm.cpp
// Date         : 3/1/18 (DLR)
// Description  : Namespace to 'wrap' MPI communication calls
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gcomm.hpp"

#if !defined(GCOMM_GLOBAL_DATA)
#define GCOMM_GLOBAL_DATA
#if defined(GEOFLOW_USE_MPI)
GINT         ngrstatus_ = 0;
GINT         ngsstatus_ = 0;
MPI_Request *mpi_recv_req_=NULLPTR;
MPI_Request *mpi_send_req_=NULLPTR;
MPI_Status  *grstatus_    =NULLPTR;
MPI_Status  *gsstatus_    =NULLPTR;
#else
GINT ngrstatus_ = 0;
GINT ngsstatus_ = 0;
void *grstatus_ = NULLPTR;
void *gsstatus_ = NULLPTR;
#endif

#endif  

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD     : InitComm
// DESC       : initializes communication layer
// ARGS       :
// RETURNS    : 
//**********************************************************************************
void GComm::InitComm(int *argc, char **argv[])
{
#if defined(GEOFLOW_USE_MPI)
  int          ierrMPI, iMPILen;
  char         sMPIErr[GMAX_ERROR_STRING];

  ierrMPI = MPI_Init(argc, argv);
  if ( ierrMPI != MPI_SUCCESS ) {
    MPI_Error_string(ierrMPI, sMPIErr, &iMPILen);
    std::cout << "GComm::InitComm: MPI_Init failed: MPI_Err: " << sMPIErr << std::endl;
    exit(1);
  }

  GINT    nprocs = GComm::WorldSize(MPI_COMM_WORLD);
  if ( mpi_recv_req_ == NULLPTR ) delete [] mpi_recv_req_; 
  if ( grstatus_     == NULLPTR ) delete [] grstatus_    ;
  ngrstatus_    = nprocs;
  mpi_recv_req_ = new MPI_Request[ngrstatus_];
  grstatus_     = new MPI_Status [ngrstatus_];
#endif

} // end of InitComm


//************************************************************************************
//************************************************************************************
// METHOD     : WorldRank
// DESC       : Gets task rank from all available tasks
// ARGS       :
// RETURNS    : rank
//************************************************************************************

GINT GComm::WorldRank(GC_COMM comm)
{
  GINT  rank=0;

#if defined(GEOFLOW_USE_MPI)
  MPI_Comm_rank(comm, &rank);
#endif

  return rank;
} // end of WorldRank


//**********************************************************************************
//**********************************************************************************
// METHOD     : WorldSize
// DESC       : Gets number of tasks
// ARGS       :
// RETURNS    : number of tasks
//**********************************************************************************

GINT  GComm::WorldSize(GC_COMM comm)
{
  GINT  ntasks=1;

#if defined(GEOFLOW_USE_MPI)
  MPI_Comm_size(comm,&ntasks);
#endif

  return ntasks;
} // end of WorldSize


//**********************************************************************************
//**********************************************************************************
// METHOD     : TermComm
// DESC       : Terminates communication layer
// ARGS       :
// RETURNS    : 
//**********************************************************************************
void GComm::TermComm()
{

#if defined(GEOFLOW_USE_MPI)
  MPI_Finalize();
  if ( mpi_recv_req_ != NULLPTR ) delete [] mpi_recv_req_; 
  if ( mpi_send_req_ != NULLPTR ) delete [] mpi_send_req_; 
  if ( grstatus_     != NULLPTR ) delete [] (MPI_Status *)grstatus_;
  if ( gsstatus_     != NULLPTR ) delete [] (MPI_Status *)gsstatus_;
#endif

} // end of TermComm


//**********************************************************************************
//**********************************************************************************
// METHOD     : ASendRecv (1)
// DESC       : Performs asynchronous send/recv from any source. RecvBuff 
//              and SendBuff assumed to be 'matrices', one row for each buffer 
//              to be sent, received.

//              The source task id's are returned in 'source' variable,
//              which must be allocated always.
// ARGS       :
//              RecvBuff  : receive buffers (GTMatrix)
//              nRecvBuff : # recv buffers (rows) in RecvBuff, # irecv & source elems
//              irecv     : indirection array giving for each recv buffer, the 
//                          column index into RecvBuff to recv from. Used if non-NULL
//              maxRecvLen: length of each RecvBuff
//              rtype     : GCommDatatype of received data
//              source    : contains sources of recvd data. If bUseSource==TRUE, must 
//                          contain the source task ids from which data is received. If
//                          bUseSource==FALSE, the library will fill it with the 
//                          sources that send the RecvBuff's their data.
//              bUseSource: if FALSE, use array source as source for posting recvs; if not,
//                          then these contain the sources from MPI
//              SendBuff  : send buffer (matrix)
//              nSendBuff : # send buffers (rows) in SendBuff
//              isend     : indirection array giving for each send buffer, the 
//                          column index into SendBuff to send. Used if non-NULL
//              maxSendLen: length of each SendBuff
//              stype     : GCommDatatype of sent data
//              dest      : destination task for sent data
//              comm      : GC_COMM
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************

GBOOL GComm::ASendRecv(void *RecvBuff, GINT  nRecvBuff, GINT  *irecv, GINT maxRecvLen, GCommDatatype rtype, GINT  *source, GBOOL bUseSource,
                       void *SendBuff, GINT  nSendBuff, GINT  *isend, GINT maxSendLen, GCommDatatype stype, GINT  *dest  , GC_COMM comm)
{

  GString serr = "GComm::ASendRecv: ";
  GINT  i, ierr;
  GD_DATATYPE irtype = GCommData2Index(rtype);
  GD_DATATYPE istype = GCommData2Index(stype);


#if defined(GEOFLOW_USE_MPI)

#if 0
  if ( RecvBuff == NULLPTR || nRecvBuff == 0 || maxRecvLen == 0 
    || SendBuff == NULLPTR || nSendBuff == 0 || maxSendLen == 0  ) return TRUE;
  // Note: may want one task to send, another to receive, so this 
  //       test isn't advised!
#endif
 
    
  GBYTE *pbeg;
  if ( mpi_recv_req_ == NULLPTR || ngrstatus_ <  nRecvBuff ) {
    if ( mpi_recv_req_ == NULLPTR ) delete [] mpi_recv_req_; 
    if ( grstatus_     == NULLPTR ) delete [] grstatus_    ;
    ngrstatus_    = nRecvBuff;
    mpi_recv_req_ = new MPI_Request[ngrstatus_];
    grstatus_     = new MPI_Status [ngrstatus_];
  }

  // Post receives:
  if ( !bUseSource ) {
    if ( irecv == NULLPTR ) {
      for ( i=0; i<nRecvBuff; i++ ) {
        pbeg = (GBYTE*)RecvBuff + i*maxRecvLen*GD_DATATYPE_SZ[irtype];
        ierr = MPI_Irecv(pbeg, maxRecvLen, rtype, MPI_ANY_SOURCE, 111, comm, &mpi_recv_req_[i]);
      }
    }
    else {
      for ( i=0; i<nRecvBuff; i++ ) {
        pbeg = (GBYTE*)RecvBuff + irecv[i]*maxRecvLen*GD_DATATYPE_SZ[irtype];
        ierr = MPI_Irecv(pbeg, maxRecvLen, rtype, MPI_ANY_SOURCE, 111, comm, &mpi_recv_req_[i]);
      }
    }
  }
  else {
    if ( irecv == NULLPTR ) {
      for ( i=0; i<nRecvBuff; i++ ) {
        pbeg = (GBYTE*)RecvBuff + i*maxRecvLen*GD_DATATYPE_SZ[irtype];
        ierr = MPI_Irecv(pbeg, maxRecvLen, rtype, source[i], 111, comm, &mpi_recv_req_[i]);
      }
    }
    else {
      for ( i=0; i<nRecvBuff; i++ ) {
        pbeg = (GBYTE*)RecvBuff + irecv[i]*maxRecvLen*GD_DATATYPE_SZ[irtype];
        ierr = MPI_Irecv(pbeg, maxRecvLen, rtype, source[i], 111, comm, &mpi_recv_req_[i]);
      }
    }
  }

 
  // Send data:
  if ( isend == NULLPTR ) {
    for ( i=0; i<nSendBuff; i++ ) {
      pbeg = (GBYTE*)SendBuff + i*maxSendLen*GD_DATATYPE_SZ[istype];
      ierr = MPI_Send(pbeg, maxSendLen, stype, dest[i], 111, comm);
    }
  }
  else {
    for ( i=0; i<nSendBuff; i++ ) {
      pbeg = (GBYTE*)SendBuff + isend[i]*maxSendLen*GD_DATATYPE_SZ[istype];
      ierr = MPI_Send(pbeg, maxSendLen, stype, dest[i], 111, comm);
    }
  }

    
  GINT  slen=256;
  GINT  irerr, iserr;
  GCHAR srerr[slen];
  GCHAR sserr[slen];

  MPI_Barrier(comm);
  MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

  ierr = MPI_Waitall(nRecvBuff, mpi_recv_req_, grstatus_);
  if ( ierr != MPI_SUCCESS ) {
    for ( i=0; i<nRecvBuff; i++ ) {
      irerr = grstatus_[i].MPI_ERROR;
      if ( irerr != MPI_SUCCESS ) {
        MPI_Error_string(irerr,srerr,&slen);
        GPP(comm, serr << " MPI_ERROR(RECV)[" << i << "]=" << srerr);
      }
    }
    MPI_Abort(comm,0);
  }

  // if we didn't use the sources, then we return
  // them:
  if ( !bUseSource ) {
    for ( i=0; i<nRecvBuff; i++ )
      source[i] = grstatus_[i].MPI_SOURCE;
  }
    
  return TRUE;

#else

  GBYTE *psend, *precv;
  if ( isend == NULLPTR ) {
    if ( irecv == NULLPTR ) {
      for ( i=0; i<nSendBuff; i++ ) {
        psend = (GBYTE*)SendBuff + i*maxSendLen*GD_DATATYPE_SZ[istype];
        precv = (GBYTE*)RecvBuff + i*maxRecvLen*GD_DATATYPE_SZ[irtype];
        memcpy(precv, psend, maxSendLen*GD_DATATYPE_SZ[istype]);
      }
    }
    else {
      for ( i=0; i<nSendBuff; i++ ) {
        psend = (GBYTE*)SendBuff + i       *maxSendLen*GD_DATATYPE_SZ[istype];
        precv = (GBYTE*)RecvBuff + irecv[i]*maxRecvLen*GD_DATATYPE_SZ[irtype];
        memcpy(precv, psend, maxSendLen*GD_DATATYPE_SZ[istype]);
      }
    }
  }
  else {
    if ( irecv == NULLPTR ) {
      for ( i=0; i<nSendBuff; i++ ) {
        psend = (GBYTE*)SendBuff + isend[i]*maxSendLen*GD_DATATYPE_SZ[istype];
        precv = (GBYTE*)RecvBuff + i       *maxRecvLen*GD_DATATYPE_SZ[irtype];
        memcpy(precv, psend, maxSendLen*GD_DATATYPE_SZ[istype]);
      }
    }
    else {
      for ( i=0; i<nSendBuff; i++ ) {
        psend = (GBYTE*)SendBuff + isend[i]*maxSendLen*GD_DATATYPE_SZ[istype];
        precv = (GBYTE*)RecvBuff + irecv[i]*maxRecvLen*GD_DATATYPE_SZ[irtype];
        memcpy(precv, psend, maxSendLen*GD_DATATYPE_SZ[istype]);
      }
    }
  }
  return TRUE;
#endif

} // end of method ASendRecv (1)



//**********************************************************************************
//**********************************************************************************
// METHOD     : BRecv
// DESC       : Performs blocking receive
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::BRecv(void *rbuff, GINT count, GCommDatatype rtype, GINT src, GC_COMM comm)
{   
  
#if !defined(GEOFLOW_USE_MPI)
  return FALSE;
#else
  GINT iret;
  MPI_Status mstatus;

  if ( rbuff == NULLPTR ) return FALSE;

  iret = MPI_Recv(rbuff, count, rtype, src, 0, comm, &mstatus);
  
  return iret == MPI_SUCCESS ? TRUE : FALSE;
#endif
    
} // end of method BRecv


//**********************************************************************************
//**********************************************************************************
// METHOD     : IRecv
// DESC       : Performs non-blocking receive
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::IRecv(void *rbuff, GINT count, GCommDatatype rtype, GINT src, void *hreq, GC_COMM comm )
{   
  
#if !defined(GEOFLOW_USE_MPI)
  return FALSE;
#else
  GINT iret;

  if ( rbuff == NULLPTR ) return FALSE;

  iret = MPI_Irecv(rbuff, count, rtype, src, 0, comm, (MPI_Request *)hreq);
  
  return iret == MPI_SUCCESS ? TRUE : FALSE;
#endif
    
} // end of method IRecv


//**********************************************************************************
//**********************************************************************************
// METHOD     : BSend
// DESC       : Performs blocking send 
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::BSend(void *sendbuff, GINT  count, GCommDatatype stype, GINT dest, GC_COMM comm )
{   

#if !defined(GEOFLOW_USE_MPI)
  return FALSE;
#else
  GINT iret; 

  if ( sendbuff == NULLPTR ) return FALSE;

  iret = MPI_Send(sendbuff, count, stype, dest, 0, comm);
  
  return iret == MPI_SUCCESS ? TRUE : FALSE;
#endif
  
    
} // end of method BSend


//**********************************************************************************
//**********************************************************************************
// METHOD     : ISend
// DESC       : Performs non-blocking send 
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::ISend(void  *sendbuff, GINT count, GCommDatatype stype, GINT dest, void *hreq, GC_COMM comm )
{   

#if !defined(GEOFLOW_USE_MPI)
  return FALSE;
#else
  GINT iret;

  if ( sendbuff == NULLPTR ) return FALSE;

  iret = MPI_Isend(sendbuff, count, stype, dest, 0, comm, (MPI_Request *)hreq);
  
  return iret == MPI_SUCCESS ? TRUE : FALSE;
#endif
    
} // end of method ISend


//**********************************************************************************
//**********************************************************************************
// METHOD     : BWaitAll (1)
// DESC       : Performs blocking wait for recv of specified handle.
//              Remember, this handle may represent multiple posts.
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::BWaitAll(GCMHandle &ch)
{   

#if !defined(GEOFLOW_USE_MPI)
  return FALSE;
#else

  if ( ngrstatus_ < ch.nposts_ ) {
    delete [] grstatus_;
    grstatus_ = new MPI_Status[ch.nposts_];
    ngrstatus_ = ch.nposts_;
  }
   
#if 0
  MPI_Wait((MPI_Request*)(ch.mhandle_)+i, grstatus_);
  for ( i=0; i<ch.nposts_; i++ ) { 
     MPI_Wait((MPI_Request*)(ch.mhandle_)+i, grstatus_);
  }
#else
  MPI_Waitall(ch.nposts_, (MPI_Request*)(ch.mhandle_), (MPI_Status *) grstatus_);
#endif

  return TRUE;

#endif

    
} // end of method BWaitAll(1)
 

//**********************************************************************************
//**********************************************************************************
// METHOD     : BWaitAll(2)
// DESC       : Performs blocking wait for recv of specified handle.
//              Remember, this handle may represent multiple posts.
// ARGS       :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GComm::BWaitAll(void *mreq, GINT nreq)
{   
  GINT iret;
#if !defined(GEOFLOW_USE_MPI)
  std::cout << "GComm::BWaitAll: GEOFLOW_USE_MPI undefined" << std::endl;
  return FALSE;
#else

  if ( ngrstatus_ < nreq ) {
    delete [] grstatus_;
    grstatus_ = new MPI_Status[nreq];
    ngrstatus_ = nreq;
  }
   
#if 0
  MPI_Wait((MPI_Request*)(ch.mhandle_)+i, grstatus_);
  for ( i=0; i<ch.nposts_; i++ ) { 
     MPI_Wait((MPI_Request*)(ch.mhandle_)+i, grstatus_);
  }
#else
  iret = MPI_Waitall(nreq, (MPI_Request*)mreq, (MPI_Status *) grstatus_);
#endif

  if ( iret != MPI_SUCCESS )
    std::cout << "GComm::BWaitAll: iret= " << iret << std::endl;

  return iret == MPI_SUCCESS ? TRUE : FALSE;

#endif

    
} // end of method BWaitAll(2)
 

//**********************************************************************************
//**********************************************************************************
// METHOD     : Allreduce
// DESC       : Performs reduction
// ARGS       :
// RETURNS    : 
//**********************************************************************************
GINT  GComm::Allreduce(void  *operand, void *result, const GINT  count, GCommDatatype itype, GC_OP iop, GC_COMM comm)
{

  GINT   ircount=count;
#if defined(GEOFLOW_USE_MPI)

  ircount = MPI_Allreduce(operand, result, count, itype, GC_Optype[iop], comm);
#else
  GD_DATATYPE igtype = GCommData2Index(itype);
  memcpy((GBYTE*)result, (GBYTE*)operand, count*GD_DATATYPE_SZ[igtype]);
#endif

  return ircount;

} // end of method Allreduce


//**********************************************************************************
//**********************************************************************************
// METHOD     : Allgather
// DESC       : Performs all-gather operation
// ARGS       :
// RETURNS    : 
//**********************************************************************************
GINT GComm::Allgather(void *operand, GINT  sendcount, GCommDatatype stype, 
                      void *result , GINT  recvcount, GCommDatatype rtype, GC_COMM comm)
{

  GINT    iret=sendcount, rank=GComm::WorldRank(comm);
  GString serr = "GComm::Allgather: ";

#if defined(GEOFLOW_USE_MPI)
  if ( operand == NULLPTR ) {
    std::cout << serr << "MPI_IN_PLACE requested, but not allowed" << std::endl;
    exit(1);
  }

  iret = MPI_Allgather(operand, sendcount, stype, result, recvcount, rtype, comm);
#else
  if ( recvcount < sendcount ) return 0;
  if ( operand == NULLPTR || result == NULLPTR ) return 0;
  GD_DATATYPE irtype = GCommData2Index(rtype);
  GD_DATATYPE istype = GCommData2Index(stype);
  memcpy((GBYTE*)result+rank*recvcount*GD_DATATYPE_SZ[irtype], 
          (GBYTE*)operand  , sendcount*GD_DATATYPE_SZ[istype]);
  
#endif

  return iret;

} // end of method Allgather


//**********************************************************************************
//**********************************************************************************
// METHOD     : Gatherv
// DESC       : Performs MPI_Gatherv operation
// ARGS       :
// RETURNS    : 
//**********************************************************************************
GINT GComm::Gatherv(void *sbuff, GINT sendcount, GCommDatatype stype, void *rbuff, GINT *recvcount, GINT *disp, GCommDatatype rtype, GINT root, GC_COMM comm)
{

  GINT    iret=sendcount, rank=GComm::WorldRank(comm);
  GString serr = "GComm::Gatherv: ";

#if defined(GEOFLOW_USE_MPI)
  if ( sbuff == NULLPTR ) {
    std::cout << serr << "MPI_IN_PLACE requested, but not allowed" << std::endl;
    exit(1);
  }

  iret = MPI_Gatherv(sbuff, sendcount, stype, rbuff, recvcount, disp, rtype, root, comm);
#else
  if ( recvcount[0] < sendcount ) return 0;
  if ( sbuff == NULLPTR || result == NULLPTR ) return 0;
  GD_DATATYPE irtype = GCommData2Index(rtype);
  GD_DATATYPE istype = GCommData2Index(stype);
  for ( auto i=0; i<sendcount; i++ ) {
    memcpy((GBYTE*)rbuff  + disp[i]*GD_DATATYPE_SZ[irtype], 
           (GBYTE*)sbuff, sendcount*GD_DATATYPE_SZ[istype]);
  } 
#endif

  return iret;

} // end of method Gatherv


//**********************************************************************************
//**********************************************************************************
// METHOD     : DataTypeFromStruct
// DESC       : Builds a derived data type usable by Comm methods, which encapsulates
//              the data associated with a C-type structure, *WHICH CONTAINS
//              CONTIGUOUS DATA*.
// ARGS       :
//              blk_ptr    : array of pointers to C-type structure members, 
//                           of length 'n-blk'. Note, that for platform-independence,
//                           the blk_ptr array should be formed by calls to
//                           GComm::Address.
//              blk_types  : array of GCommDatatype's of length 'n_blk', each of
//                           member of which represents the GCommDatatype's of the
//                           c_struct data members in order. The sz_types MUST
//                           NOT be G_STRUCTTYPE, but MUST be a basic type!
//              n_types    : array of length 'n_blk' that gives the corresp. no. of
//                           elements of type sz_types.
//              n_blk      : number of data members in the c_struct structure
//              return_type: the GCommDatatype which will be passed in subsequent 
//                           GComm calls, is an MPI_Datatype when MPI is used, else
//                           it is the structure size (an GINT ), in bytes.
// RETURNS    : TRUE on success; else FALSE.
//**********************************************************************************
GBOOL GComm::DataTypeFromStruct(AGINT  offset[], GCommDatatype blk_types[], GINT  n_types[], 
                                const GINT n_blk, GCommDatatype *return_type)
{
#if defined(GEOFLOW_USE_MPI)
  GBOOL        bret=TRUE;
  GINT         i, ierrMPI, n  ;
  int          iMPILen;
  char         sMPIErr[GMAX_ERROR_STRING];

  // Build data structure for blocks type sz_types[i]:

  // Build MPI datatype:
#if defined(G_MPI1)
  ierrMPI=MPI_Type_struct       (n_blk, n_types, offset, blk_types, (MPI_Datatype*)return_type);
#else
  ierrMPI=MPI_Type_create_struct(n_blk, n_types, offset, blk_types, (MPI_Datatype*)return_type);
#endif
  if ( ierrMPI != MPI_SUCCESS ) {
    MPI_Error_string(ierrMPI, sMPIErr, &iMPILen);
    std::cout << "GComm::DataTypeFromStruct: MPI_Type_create_struct failed: MPI_Err: " << sMPIErr << std::endl;
    bret = FALSE;
  }

  return bret;
#else
  GINT  i, n, r;

  r = 0;
  n = 0;
  for ( i=0; i<n_blk; i++ ) {
    n += blk_types[i]*n_types [i];
  }
  // Now add the size of padding required to align
  // on integer-byte boundaries:
  r = n % GWORDSIZE_BYTES;
  n = (r == 0 ? n :  n + GWORDSIZE_BYTES - r); 
  *return_type = (GCommDatatype)n; 
  return TRUE;
#endif

} // end of method DataTypeFromStruct


//**********************************************************************************
//**********************************************************************************
// METHOD     : DataTypeCommit
// DESC       : Commits datatype
// ARGS       :
//              type: the GCommDatatype that will be passed in subsequent 
//                    GComm calls, is an MPI_Datatype when MPI is used, else
//                    it is the structure size (an GINT ), in bytes.
// RETURNS    : TRUE on success; else FALSE.
//**********************************************************************************
GBOOL GComm::DataTypeCommit(GCommDatatype *type)
{
#if defined(GEOFLOW_USE_MPI)
  GINT         i, ierrMPI, n  ;
  int          iMPILen;
  char         sMPIErr[GMAX_ERROR_STRING];
  ierrMPI=MPI_Type_commit((MPI_Datatype*)type);
  if ( ierrMPI != MPI_SUCCESS ) {
    MPI_Error_string(ierrMPI, sMPIErr, &iMPILen);
    std::cout << "GComm::DataTypeCommit: MPI_Type_commit failed: MPI_Err: " << sMPIErr << std::endl;
    return FALSE;
  }
#endif
  return TRUE;
} // end of method DataTypeCommit


//**********************************************************************************
//**********************************************************************************
// METHOD     : DataTypeResized
// DESC       : resizes true extent of datatype
// ARGS       :
//              ntype: new datatype
//              otype: old datatype
//              lb   : lower bound
//              ub   : upper bound
// RETURNS    : TRUE on success; else FALSE.
//**********************************************************************************
GBOOL GComm::DataTypeResized(GCommDatatype *ntype, GCommDatatype otype, GINT lb, GINT ub)
{
#if defined(GEOFLOW_USE_MPI) && !defined(G_MPI1)
  GINT         ierrMPI  ;
  int          iMPILen;
  char         sMPIErr[GMAX_ERROR_STRING];
  ierrMPI = MPI_Type_create_resized(otype, lb, ub, ntype);
  if ( ierrMPI != MPI_SUCCESS ) {
    MPI_Error_string(ierrMPI, sMPIErr, &iMPILen);
    std::cout << "GComm::DataTypeCommit: MPI_Type_commit failed: MPI_Err: " << sMPIErr << std::endl;
    return FALSE;
  }
#else
  *ntype = otype;
#endif
  return TRUE;
} // end of method DataTypeResized

//**********************************************************************************
//**********************************************************************************
// METHOD     : DataTypeFree
// DESC       : Frees derived data types created from DataTypeFromStruct method.
// ARGS       :
//              type: the datatype to be freed. Can be MPI or other
// RETURNS    : TRUE on success; else FALSE.
//**********************************************************************************
void GComm::DataTypeFree(GCommDatatype *type)
{
#if defined(GEOFLOW_USE_MPI)

  // Free MPI datatype:
  MPI_Type_free(type);
#endif

} // end of method DataTypeFree


//**********************************************************************************
//**********************************************************************************
// METHOD     : Address
// DESC       : Gets implementation--independent address of memory location 
// ARGS       : location: pointer location of interest
//              address : int address or Address Integer, AGINT 
// RETURNS    :
//**********************************************************************************
void GComm::Address(void *location, AGINT  *address)
{   
#if defined(GEOFLOW_USE_MPI)
//std::cout << " GComm::Address: location=" << location << " address=" << *address << std::endl;
#  if defined(G_MPI1)
  MPI_Address(location, address);
#  else
  MPI_Get_address(location, address);
#  endif
//std::cout << " GComm::Address: location=" << location << " address=" << *address << std::endl;
#else
  *address = *((AGINT*)(location));
#endif
} // end of method Address
  

//**********************************************************************************
//**********************************************************************************
// METHOD     : Synch(1)
// DESC       : Synchronozises all tasks
// ARGS       :
// RETURNS    :
//**********************************************************************************
void GComm::Synch(GC_COMM comm)
{
#if defined(GEOFLOW_USE_MPI)
  MPI_Barrier(comm);
#endif
}// end of Sync(1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : GCommData2Index
// DESC       : Get GD_DATATYPE index from GCommDatatype type
// ARGS       :
// RETURNS    :
//**********************************************************************************
GD_DATATYPE GComm::GCommData2Index(GCommDatatype ctype)
{
  GINT j=0;
#if defined(GEOFLOW_USE_MPI)
  while ( j<=GTYPE_NUM && ctype != GC_DATATYPE[j] ) j++;
  if ( j >= GTYPE_NUM ) assert(FALSE);
  return (GD_DATATYPE)j;
#else
  while ( j<=GTYPE_NUM && ctype != j ) j++;
  if ( j >= GTYPE_NUM ) assert(FALSE);
  return (GD_DATATYPE)j;
#endif
} // end, method GCommData2Index


//**********************************************************************************
//**********************************************************************************
// METHOD     : cout
// DESC       : Print all SPMD input strings to stdout
// ARGS       : string to output (may use ostringstream to create)
// RETURNS    : none.
//**********************************************************************************
void GComm::cout(GC_COMM comm, const char *pstr)
{
  GString xstr(pstr);

#if defined(GEOFLOW_USE_MPI)
  GString buff;
  GINT    myrank = GComm::WorldRank(comm);
  GINT    nprocs = GComm::WorldSize(comm);
  GSIZET  glen, len = xstr.length();
  GCommDatatype dtype = T2GCDatatype<GCHAR>();
  MPI_Status mstatus;

  GComm::Allreduce(&len, &glen, 1, T2GCDatatype<GSIZET>(), GC_OP_MAX, comm);

  buff.resize(glen);

//GBOOL GComm::BRecv(void *rbuff, GINT count, GCommDatatype rtype, GINT src, GC_COMM comm)
  if ( myrank == 0 ) {
    std::cout << myrank  << ": " << xstr << std::endl << std::flush;
    for ( GINT k=1; k<nprocs; k++ ) {
      MPI_Recv((void*)buff.data(), glen, dtype, k, 66666, comm, &mstatus); 
      std::cout << k << ": " << buff << std::endl << std::flush;
    }
  }
  else {
    MPI_Send((void*)pstr, len, dtype, 0, 66666, comm);
  }

#else
  std::cout << xstr << std::endl;
#endif
} // end, method cout


