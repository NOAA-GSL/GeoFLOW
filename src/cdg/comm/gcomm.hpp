//==================================================================================
// Module       : gcomm.hpp
// Date         : 3/1/18 (DLR)
// Description  : Namespace to 'wrap' MPI communication calls
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(GCOMM_HPP)
#define GCOMM_HPP

#include "gtypes.h"
#include "gcommdata_t.h"

#include <cstdlib>
#include <iostream>
#include <string>

#if defined(GEOFLOW_USE_MPI)
#  include "mpi.h"
#endif


#if defined(GEOFLOW_USE_MPI)
extern   GINT ngrstatus_;
extern   GINT ngsstatus_;
extern   MPI_Op G2MPI_OPTYPE[]  ;
extern   MPI_Request *mpi_recv_req_;
extern   MPI_Request *mpi_send_req_;
extern   MPI_Status  *grstatus_    ;
extern   MPI_Status  *gsstatus_    ;
#else
extern GINT ngrstatus_;
extern GINT ngsstatus_;
extern void *grstatus_;
extern void *gsstatus_;
#endif


namespace GComm
{

                       void     InitComm   (int *argc, char **argv[]);
                       void     TermComm   ();
                       GINT     WorldRank       (GC_COMM icomm=GC_COMM_WORLD);
                       GINT     WorldSize       (GC_COMM icomm=GC_COMM_WORLD);
                       GBOOL    ASendRecv  (void *RecvBuff, GINT nRecvBuff, GINT  *irecv, GINT RecvLen   , GCommDatatype rtype, GINT *source, GBOOL bUseSource, 
                                            void *SendBuff, GINT nSendBuff, GINT  *isend, GINT maxSendLen, GCommDatatype stype, GINT *dest, GC_COMM icomm=GC_COMM_WORLD   );
                       GINT Allreduce  (void *, void *, const GINT  count, GCommDatatype type, GC_OP op, GC_COMM icomm=GC_COMM_WORLD);
                       GINT Allgather(void *operand, GINT  sendcount, GCommDatatype stype, void *result , GINT  recvcount, GCommDatatype rtype, GC_COMM comm);
                       GINT Gatherv  (void *sbuff, GINT sendcount, GCommDatatype stype, void *rbuff, GINT *recvcount, GINT *disp, GCommDatatype gtype, GINT root, GC_COMM icomm=GC_COMM_WORLD);

                       GBOOL    BSend      (void *sbuff, GINT  buffcount, GCommDatatype stype, GINT dest, GC_COMM icomm=GC_COMM_WORLD  );
                       GBOOL    ISend      (void *sbuff, GINT  buffcount, GCommDatatype stype, GINT dest, void *hreq, GC_COMM icomm=GC_COMM_WORLD);
                       GBOOL    BRecv      (void *rbuff, GINT  buffcount, GCommDatatype stype, GINT dest, GC_COMM icomm=GC_COMM_WORLD);
                       GBOOL    IRecv      (void *rbuff, GINT  buffcount, GCommDatatype stype, GINT dest, void *hreq, GC_COMM icomm=GC_COMM_WORLD);
                       GBOOL    DataTypeFromStruct(AGINT  blk_ptr[], GCommDatatype blk_types[], GINT  n_type[],
                                const GINT  num_typ, GCommDatatype *return_type);
                       GBOOL    DataTypeCommit(GCommDatatype *type);
                       GBOOL    DataTypeResized(GCommDatatype *ntype, GCommDatatype otype, GINT lb, GINT ub);

                       void     DataTypeFree(GCommDatatype *type);
                       void     Address(void *location, AGINT  *address);
                       void     Synch(GC_COMM icomm=GC_COMM_WORLD);
                       GBOOL    BWaitAll(GCMHandle &handle);
                       GBOOL    BWaitAll(void *handle, GINT nhandle);
                       void     cout(GC_COMM comm, const char *);
  
                       GD_DATATYPE GCommData2Index(GCommDatatype ctype);

  // private data:

}

#endif

