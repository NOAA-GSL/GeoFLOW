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
// METHOD : Constructor
// DESC   : 
// ARGS   : GD_COMM object, defaults to GD_COMM_WORLD
// RETURNS: GGFX
//************************************************************************************
template<typename T>
GGFX<T>::GGFX(GC_COMM icomm)
:
comm_                 (icomm),
nprocs_               (1),
maxNodeVal_           (0),
rank_                 (0),
nglob_index_          (0),
bInit_                (FALSE),
bBinSorted_           (FALSE)
{
  rank_   = GComm::WorldRank(comm_);
  nprocs_ = GComm::WorldSize(comm_);

} // end of constructor (2) method


//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename T>
GGFX<T>::~GGFX()
{
}


//************************************************************************************
//************************************************************************************
// METHOD : Copy constructor
// DESC   : 
// ARGS   : GGFX
// RETURNS: none
//************************************************************************************
// Copy constructor method
template<typename T>
GGFX<T>::GGFX(const GGFX<T> &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GGFX
// RETURNS: GGFX
//************************************************************************************
template<typename T>
GGFX<T>  &GGFX<T>::operator=(const GGFX<T> &a)
{

  return *this;
 
} // end of = operator



//************************************************************************************
//************************************************************************************
// METHOD : init
// DESC   : Performs initialization of global gather-scatter operation, for
//          nodes sorted only by MPI task id. This init is meant to be
//          used with the doOp methods. 
// ARGS   : glob_index: global index list
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::init(GNIDBuffer &glob_index)
{
  GBOOL      bret;
  GString   serr = "GGFX<T>::init: ";

  nglob_index_ = glob_index.size();

  // Get node id dynamic range:
  GNODEID lmax = glob_index.max();
  GComm::Allreduce(&lmax, &maxNodeVal_, 1, T2GCDatatype<GNODEID>(), GC_OP_MAX, comm_);
 
  // Do initial sorting of all sortable data. 
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "MaxIndexDynamicRange=" << maxNodeVal_);
#endif
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "initSort...");
#endif
  bret = initSort(glob_index);
  if ( !bret ) {
    GPP(comm_,serr << "initSort failed");
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "initSort done.");
#endif

  bInit_ = TRUE;

  // With all comm data computed, find
  // inv multiplicity vector:
  initMult();

  return bInit_;

} // end of method init


//************************************************************************************
//************************************************************************************
// METHOD: initSort
// DESC  : Performs sorting for initialization of global gather-scatter operation.
// ARGS   : glob_index: global index list input to Init
//            Note:
//            iOpL2RTasks     : list of task ids to send rows of iOpIndices to;
//                             receive from
//            iOpL2RIndices   : matrix returned , 1 row for each task in iOpTasks, 
//                             containing columns of local indices into the glob_index array:
//                  row 0, task 0: index0 index1 index2....
//                  row 1, task 1: index0 index1 index2....
//                  row 2, task 2: index0 index1 index2....
//                             This data is used to send data.
//            nOpL2RIndices : number of OpSRIndices for each task
//            iOpL2LIndices: matrix of indices into local glob_index; each column
//                           represents a shared (global) index, and each row represents 
//                           the index to the same shared index 
//            nOpL2LIndices: for each column in iOpL2LIndices, the number
//                           of local indices representing the same shared node.
//                ...
// RETURNS: 
//          TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::initSort(GNIDBuffer  &glob_index)
{
  GString serr = "GGFX<T>::initSort: ";
  GBOOL   bret;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "glob_index=" << glob_index);
#endif

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Entering sort sequence...");
#endif
  GINT        max_numlocfilledbins; // max over tasks of # filled bins it has data for
  GINT        numlocfilledbins;     // no bins filled with local node data=max # tasks to recv from/send to
  GINT        max_numlocbinmem;     // max over tasks of # nodes in a bin
  GIMatrix    gBinMat;              // For each task, its filled bins and the # nodes in each
  GNIDBuffer  locWork;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "Calling binSort...");
#endif

  bret = binSort(glob_index, gBinMat, numlocfilledbins,
                 max_numlocfilledbins, max_numlocbinmem, gBinBdy_, locWork);
  if ( !bret ) {  
    GPP(comm_,serr << "binSort failed ");
    exit(1);
  } 

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "gBinMat =" << gBinMat);
#endif

  // Create/resize work buffer:
  GIBuffer   iRecvWorkTaskID;// list of tasks (not incl this one) to receive work from
  GIBuffer   iSendWorkTaskID;// list of tasks to send work to
  GNIDMatrix irWork;         // Wrk recv buff: (iRecvWorkTaskID.size + 1) x max_numlocbinmem (+1 for this rank)

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Calling createWorkBuffs...");
#endif

  GNIDMatrix isWork;         // Wrk send buff: iSendWorkTaskID.size x max_numlocbinmem
  bret = createWorkBuffs(gBinMat, max_numlocbinmem, iRecvWorkTaskID, irWork, iSendWorkTaskID, isWork);
  if ( !bret ) {  
    GPP(comm_, serr << "createWorkBuffs failed ");
    exit(1);
  } 

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Filling irWork...");
#endif
  // Add local work to the work recvd from other ranks:

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Calling doSendRecvWork...");
#endif

  // Fill bins with global nodes & send out, and receive work:
  bret = doSendRecvWork(glob_index, gBinBdy_, iRecvWorkTaskID, irWork, iSendWorkTaskID, isWork);
  if ( !bret ) {  
    GPP(comm_,serr << "doSendRecvWork failed ");
    exit(1);
  } 

  // At this point, we have a 'recv' work array, irWork that 
  // contains this bin's or rank's unsorted work, to sort. The sorted data
  // should be for each task, it's nodes that have a multiplicty>1.
  // This data should be sent back to each task, and each task
  // should then create a local index array that finds each of these
  // mult > 1 nodes in the original glob_index array.
  GNIDMatrix mySharedData; // recv buffer for sorted work data
  bret = doCommonNodeSort(glob_index, irWork, iRecvWorkTaskID, iSendWorkTaskID, mySharedData);
  // NOTE: iSendWorkTaskID on exit should contain the list of tasks that sorted work
  //       data is received from. These should be the same task ids that were used
  //       in sending this local task's data to the work tasks it identified
  if ( !bret ) {  
    GPP(comm_,serr << "doCommonNodeSort failed ");
    exit(1);
  } 

  // Finally, find list of local indices and for send/receive (SR) operations, and
  // for purely local operations:
  bret = extractOpData(glob_index, mySharedData);
  if ( !bret ) {  
    GPP(comm_,serr << "extractOpData failed ");
    exit(1);
  } 

  
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "initSort done." );
#endif
  
  return TRUE;

} // end of method initSort 


//************************************************************************************
//************************************************************************************
// METHOD : binSort 
// DESC   : Sorts nodelist into bins based on dynamic range, maxNodeVal_. 
//          Bin boundaries are also computed and stored in member data.
// ARGS   : nodelist         : input node list to sort
//          gBinMat          : result of this method: Contains matrix of 1 col for
//                             each rank in comm_, each col of the form:
//            1st nonempty bin id (rank), #nodes binned; 2nd nonempty bin, # nodes binned ... 
//                              Matrix is resized in this method.
//          numlocfilledbins : number of bins local nodelist added members to
//          max_numfilledbins: max number of bins filled by any task
//          max_numlocbinmem : max number of local bin members among all tasks
//          gBinBdy          : Bin boundary matrix, returned
//          locWork          : list of local node list that this task will work on
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::binSort(GNIDBuffer &nodelist, GIMatrix &gBinMat, 
                    GINT &numlocfilledbins, GINT &max_numfilledbins, GINT &max_numlocbinmem, 
                    GNIDMatrix &gBinBdy, GNIDBuffer &locWork)
{
  GString  serr = "GGFX<T>::binSort: ";
  GSIZET   i, j, n;

  locWork.resize(nodelist.size());
  locWork.set(-1);  // not all will be this task's work-nodes

  gBinBdy.resize(nprocs_,2);

  // Node list can be in any order:
  GNODEID nDel = (maxNodeVal_+1) / nprocs_;
  GNODEID nRem = (maxNodeVal_+1) % nprocs_;
  GNODEID nRange0 = 0;
  GNODEID nRange1;
  for ( i=0, n=0; i<nprocs_; i++ ) { // cycle over bins (tasks)
    nRange1  = nRange0 + nDel + ( (i<nRem || (i==nprocs_-1)) ? 1:0) - 1; // set range owned by task
    gBinBdy(i,0) = nRange0; gBinBdy(i,1) = nRange1;               // set global node bin ranges  
    nRange0 = nRange1 + 1;
  }

#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "gBinBdy=" << gBinBdy);
#endif

  // Find number nodes in each bin from local task:
  GIBuffer nInBin (nprocs_); nInBin   .set(0);
  GINT     ibin;
  GINT     bNodeInBin;
  GSIZET   index;
  GNODEID  inode;
  for ( j=0; j<nodelist.size(); j++ ) { // cycle over node ids--inefficient
    inode = nodelist[j];
    ibin = 0;
    while ( ibin < nprocs_ && 
           !(inode >= gBinBdy(ibin,0) && inode <= gBinBdy(ibin,1)) ) ibin++;
    bNodeInBin = ibin < nprocs_;
    if ( !bNodeInBin ) {
       GPP(comm_,serr << "Node id not found: inode=" << inode << " ibin=" << ibin);
       return FALSE;
    }
    nInBin[ibin] += bNodeInBin;
    if ( ibin == rank_ && bNodeInBin>0 ) {
      locWork[n++] = nodelist[j]; // local node ids managed by this rank 
    }
    nRange0 = nRange1 + 1;
  }

  // Find max number of nodes in each bin among all ranks:
  GINT lmax = nInBin.max();
  GComm::Allreduce(&lmax, &max_numlocbinmem, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "max_numlocbinmem=" << max_numlocbinmem);
#endif

  // Find max number of bins filled by any rank:
  numlocfilledbins = 0;
  for ( j=0; j<nprocs_; j++ ) numlocfilledbins += nInBin[j] > 0 ? 1 : 0;
  GComm::Allreduce(&numlocfilledbins, &max_numfilledbins, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);

  // Set bin data for each rank:
  // This contains the basic 'graph' data: which ranks must send data to and receive 
  // data from other ranks in order to sort nodes. Each col of the matrix 
  // contains the data for that column position's task, of the form:
  //   1st nonempty bin id (rank), #nodes binned in 1; 2nd nonempty bin, # nodes binned in 2... 
  // tells us that that task must send data to num. nodes to rank represented
  // by 1st nonempty bin, and no. nodes binned in 2nd bin to rank represented by
  // 2nd bin, etc.
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "max_numfilledbins=" << max_numfilledbins);
#endif
  gBinMat.resize(2*max_numfilledbins,nprocs_);
  gBinMat = -1;

  // Fill in the local data for this rank:
  GIBuffer loc_record(2*max_numfilledbins); 
  loc_record.set(-1);
  for ( j=0, n=0; j<nprocs_; j++ ) {
    if ( nInBin[j] > 0 ) {
       loc_record  [n] = j;
       loc_record[n+1] = nInBin[j];
       n+=2;
    }
  }
  GCommDatatype dtype = T2GCDatatype<GINT>();
  GComm::Allgather(loc_record.data(), 2*max_numfilledbins, dtype, 
                   gBinMat.data().data(), 2*max_numfilledbins, dtype, comm_);



#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "loc_record=" << loc_record);
  GPP(comm_,serr << "gBinMat =" << gBinMat);
#endif

  bBinSorted_ = TRUE;

  return TRUE;
} // end of method binSort


//************************************************************************************
//************************************************************************************
// METHOD : createWorkBuffs
// DESC   : Create or resize work buffers
//          
// ARGS   : gBinMat         : work matrix from binSort
//          max_numlocbinmem: max num local bin members from binSort
//          iRecvWorkTaskID: Task ids this task received work from; final entry should be
//                            this rank_
//          irWork          : Work buffer containing received work from all tasks;
//                            final column should contain data from this task.
//          iSendWorkTaskID : Task ids to send work to
//          isWork          : Work send buffer            
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::createWorkBuffs(GIMatrix &gBinMat, GINT max_numlocbinmem, GIBuffer &iRecvWorkTaskID, 
                            GNIDMatrix &irWork, GIBuffer &iSendWorkTaskID, GNIDMatrix &isWork)
{
  GString  serr = "GGFX<T>::createWorkBuffs: ";
  GSIZET   i, j;

  if ( !bBinSorted_ ) return FALSE;


  // Recall, Each row of the matrix contains the data for the task given by the
  // row position for the nodes that task has:
  //   1st nonempty bin id (rank), #nodes binned in 1; 2nd nonempty bin, # nodes binned in 2... 
  // or
  //   B1 NB1; B2 NB2; B3, NB3; ...

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Compute numRecvWorkTasks...");
#endif

  // Find no. tasks we will receive data from or send to:
  GBOOL    bret;
  GINT     numRecv = 0;// How many tasks recvd from
  GINT     gnumRecv;
  GINT     numSend = 0;// How many tasks sent to
  GINT     ibintask;
  GSIZET   iwhere;
  GIBuffer ittmp(nprocs_);  // Hold unique task ids

  // Find list of tasks to receive from (include rank_):
  for ( j=0; j<nprocs_; j++ ) { // j is task that might send work to this one
    for ( i=0; i<gBinMat.size(1) && gBinMat(i,j)!=-1; i+=2 ) {
      ibintask = gBinMat(i,j);
      if ( ibintask != rank_ ) continue;
      bret = ittmp.containsn(j,numRecv,iwhere); 
      if ( !bret ) { // if bin/task not in tmp buffer, add it
        ittmp[numRecv] = j;
        numRecv++;
      }
    }
  }
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "numRecv=" << numRecv);
#endif
  iRecvWorkTaskID.resize(numRecv);  
  iRecvWorkTaskID.set(ittmp.data(),numRecv); 

  // Find number tasks to send to (include rank_):
  for ( i=0; i<gBinMat.size(1) && gBinMat(i,rank_)!=-1; i+=2 ) {
    ibintask = gBinMat(i,rank_);
    bret     = ittmp.containsn(ibintask,numSend,iwhere); 
    if ( !bret ) { // if task not in buffer, add it
      ittmp[numSend] = ibintask;
      numSend++;
    }
  }
  iSendWorkTaskID.resize(numSend);  
 
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "iRecvWorkTaskID=" << iRecvWorkTaskID);
  GPP(comm_,serr << "numRecv=" << numRecv);
  GPP(comm_,serr << "numSend=" << numSend);
#endif
  GComm::Allreduce(&numRecv, &gnumRecv, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);

  // Create node 'receive' buffers for work received from other tasks, 
  // not including this task:
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "max_numlocbinmem=" << max_numlocbinmem
        << " gnumRecv=" << gnumRecv << " numSend=" << numSend);
#endif
  irWork.resize(max_numlocbinmem,gnumRecv); 
  irWork.set(-1); 

  // Create node 'send' buffers for work sent to other tasks:
  isWork.resize(max_numlocbinmem,numSend); 
  irWork.set(-1);

  return TRUE;

} // end, method createWorkBuffs


//************************************************************************************
//************************************************************************************
// METHOD : binWorkFill
// DESC   : Fills work bins, bins, with this task's binned node data to be sent
//          for processing
// ARGS   : nodelist: input node list to sort
//          gBinBdy : bin boundary ranges computed in binSort
//          sWork   : contains for each bin (task), the members of nodelist
//                    corresponding to this bin. Allocated here.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::binWorkFill(GNIDBuffer &nodelist, GNIDMatrix &gBinBdy, GNIDMatrix  &sWork, GIBuffer &iSendWorkTaskID)
{
  GString serr = "GGFX<T>::binWorkFill: ";
  GSIZET  i, iwhere, j, n;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "Entering binWorkFill...");
#endif
  GSZBuffer ntmp(nprocs_);
  ntmp.set(0);

  // Fill work send matrix & corresp send-tasks:
  sWork.set(-1);
  iSendWorkTaskID.set(-1);
  for ( i=0, n=0; i<nodelist.size(); i++ ) {
    for ( j=0; j<nprocs_; j++ ) { // cycle over bins (equiv to procs)
      if ( nodelist[i] >= gBinBdy(j,0) && nodelist[i] <= gBinBdy(j,1) ) {
        if ( !iSendWorkTaskID.containsn(j,n,iwhere) ) {
          iSendWorkTaskID[n] = j;
          sWork(ntmp[n],n) = nodelist[i];
          ntmp[n]++;
          n++;
        }
        else {
          sWork(ntmp[iwhere],iwhere) = nodelist[i];
          ntmp[iwhere]++;
        }
      }
    }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "done.");
#endif

  return TRUE;
 
} // end of method binWorkFill


//************************************************************************************
//************************************************************************************
// METHOD : doSendRecvWork
// DESC   : 
//          glob_index     : array of this task's node ids
//          gBinBdy        : bin boundarys, from binSort
//          iRecvWorkTaskID: array of task ids to receive work from
//          irWork         : Matrix of received work: num iRecvWorkTaskID.size+1 X Max no indices
//          iSendWorkTaskID: array of task ids to send work to
//          isWork         : Matrix of work to send: num iSendWorkTaskID.size+1 X Max no indices
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::doSendRecvWork( GNIDBuffer &glob_index      , GNIDMatrix &gBinBdy, 
                            GIBuffer   &iRecvWorkTaskID , GNIDMatrix &irWork, 
                            GIBuffer   &iSendWorkTaskID , GNIDMatrix &isWork)
{
  GString serr = "GGFX<T>::doSendRecvWork: ";
  GINT    i, j, m;
  GBOOL   bret=TRUE;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Calling binWorkFill...");
#endif

  // Fill work buffers:
  if ( !binWorkFill(glob_index, gBinBdy, isWork, iSendWorkTaskID) ) {
    GPP(comm_,serr << "binWorkFill failed");
    exit(1);
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "iSendWorkTaskID=" << iSendWorkTaskID );
  GPP(comm_,serr << "isWork=" << isWork);
  GPP(comm_,serr << "Doing ASenRecv call...");
#endif

  // Set send buffer indirection, and prevent
  // send and receive to the same task; place
  // this task's receive buff (if any) as final
  // real recv buffer.
  GINT     nsend;
  GINT     nrecv;
  GINT     maxsend = iSendWorkTaskID.size();
  GINT     maxrecv = iRecvWorkTaskID.size();
  GSIZET   iwhere;
  GIBuffer ibsend(isWork.size(2));         // send buffer indirection list
  GIBuffer itrecv(iRecvWorkTaskID.size()); // receive-from tasks
  GIBuffer itsend(iSendWorkTaskID.size()); // send-to tasks

  itrecv = iRecvWorkTaskID;
  itsend = iSendWorkTaskID;

  nrecv = itrecv.size();
  if ( iRecvWorkTaskID.contains(rank_,iwhere) ) {
    itrecv[iRecvWorkTaskID.size()-1] = rank_;
    for ( j=0, m=0; j<iRecvWorkTaskID.size(); j++ ) {
      if ( iRecvWorkTaskID[j] != rank_ ) {
        itrecv[m++] = iRecvWorkTaskID[j];
      }
    }
    nrecv = itrecv.size()-1;
  }
  iRecvWorkTaskID.set(itrecv.data(), itrecv.size());

  for ( j=0, nsend=0; j<maxsend; j++ ) {
    if ( iSendWorkTaskID[j] != rank_ ) {
      itsend[nsend] = iSendWorkTaskID[j]; // task id to send to
      ibsend[nsend] = j;                  // buffer id to send (column)
      nsend++;
    }
    else {
      // fill last (maxrecv) receive buffer with work it would send
      // to itself:
      for ( i=0; i<isWork.size(1); i++ ) irWork(i,maxrecv-1) = isWork(i,j);
    }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "nsend=" << nsend << " ibsend=" << ibsend);
  GPP(comm_,serr << "itsend=" << itsend);
  GPP(comm_,serr << "nrecv=" << nrecv << " itrecv=" << iRecvWorkTaskID );
#endif
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "done. Beginning irWork=" << irWork);
#endif

  // Send work data to respective procs; wait for receive:
  GCommDatatype dtype = T2GCDatatype<GNODEID>();
  bret = GComm::ASendRecv(irWork.data().data(),nrecv,NULLPTR      ,(GINT)irWork.size(1),dtype,itrecv.data(), TRUE, 
                          isWork.data().data(),nsend,ibsend.data(),(GINT)isWork.size(1),dtype,itsend.data(), comm_);
  if ( !bret ) {
    GPP(comm_,serr << "GComm::ASendRecv failed");
    exit(1);
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "done. Final irWork=" << irWork);
#endif


  return bret;

} // end, method doSendRecvWork


//************************************************************************************
//************************************************************************************
// METHOD : doCommonNodeSort
// DESC   : Sorts the nodes from the work matrkx, irWork into quantities that
//          are usable to actually send and receive data from tasks that share
//          physical nodes. The results of this work are sent to all the tasks
//          that also share these nodes. The final results consist of 
//          a list of buffers (matrix columns ) that contain local nodes that are shared
//          by other tasks, and those that are shared locally, as well as a list
//          of task ids that we must send to and receive from data at those node ids
//
// ARGS   : 
//          glob_index      : list of global indices in call to Init
//          irWork          : 'work' array, from initSort method. This contains iRecvWorkTaskID.size)+1
//                            columns indicating tasks whose ids are given in iRecvWorkTaskID, sent this
//                            task their bin data and each column is a list of node points that
//                            were sent to the bin represented by this task. We assume that the 
//                            results of this common node sort will be sent back to these same
//                            tasks, even if there are no shared nodes in the data sent back.
//                            in initSort.
//          iRecvWorkTaskID: list of tasks from which work was received. The final entry is this 
//                            task's rank. This is the same set of task ids to send sorted work
//                            back to.
//          iSendWorkTaskID : list of tasks to which local bin work was sent. This is the same
//                            set of tasks from which we will receive sorted work, but must contain 
//                            this task as well in final element
//          mySharedData   :  contains result of this call: a matrix of column records that give
//                            global node id and the tasks that share it, of the form
//                               nodeid0 NTask0 TaskID00 TaskID01 ...; nodeid1 NTask1 TaskID10 TaskID11...
//                            which tell the local task which nodes it has that are shared with 
//                            other tasks.
//                            
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::doCommonNodeSort(GNIDBuffer &glob_index, GNIDMatrix &irWork, 
                             GIBuffer &iRecvWorkTaskID, GIBuffer &iSendWorkTaskID,
                             GNIDMatrix &mySharedData)
{
  GString    serr = "GGFX<T>::doCommonNodeSort: ";
  GSIZET     i, j, k, nr, ns;
  GBOOL      bret=TRUE;
  GINT       itask, nl;

  if ( !bBinSorted_ ) return FALSE;

  GSIZET     ipos, iwhich, rd, nd, nrow, nt, nvals=0;
  GSIZET     mult, *vvals=0, *ivals=0;
  GNODEID    nid;
  GSZBuffer  irow; // holds distinct val row indices
  GSZBuffer  icol ;// holds distinct val col indices
  GSZBuffer  irowr;// holds distinct val row indices
  GNIDBuffer mval; // holds distinct vals in rows of matrix
  GNIDBuffer nworktmp;

  // For each task in iRecvWorkTaskID from which we received work data 
  // get list of its shared nodes, and which tasks they sit on. These are
  // constained in the sendShNodeWrk matrix: each row corresponds to a task
  // id in iRecvWorkTaskID, and contains the following info:
  //       nodeid0 NTask0 TaskID00 TaskID01 ...; nodeid1 NTask1 TaskID10 TaskID11...
  // where nodeid0 is the first shared node id; NTask0 is the number of tasks this node resides in,
  // TasnID0j are the list of task ids that the node resides in, then the same for
  // node id 1 etc.

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "irWork=" << irWork);
#endif

#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "irwork.size1=" << irWork.size(1)
        << " irwork.size2=" << irWork.size(2));
#endif
  // First, get sizes for this shared node matrix:
  ipos = 0;
  nworktmp.resize(irWork.size(1)*irWork.size(2));
  nworktmp.set(-1); nt = 0;
  for ( j=0, nrow=0; j<irWork.size(2); j++ ) { // cycle thru all tasks that sent work to this task
    nd = irWork.distinctcol_floor(j, mval, irowr, -1); // get distinct node ids for this col/task
    for ( i=0; i<nd; i++ ) { 
      nid = irWork(irowr[i],j);
      mult = irWork.multiplicity(nid, irow, icol); //
      if ( mult > 1 ) { // keep only nodes with mult > 1 
        // Keep only the unique column indices:
        rd     = icol.distinctrng(0,mult,1,vvals,ivals,nvals); 
        if ( !nworktmp.contains(nid,iwhich) ) { // don't allow duplication of nodes
          ipos  += rd + 2; // required size of each entry
          nworktmp[nt] = nid;
          nt++;
        }
      }
    } 
    nrow = MAX(nrow,ipos); // Max no. rows in shared node matrix
  }

  GSIZET lsz[2] = {nrow, iSendWorkTaskID.size()};
  GSIZET gsz[2];
  GComm::Allreduce(lsz, &gsz, 2, T2GCDatatype<GSIZET>() , GC_OP_MAX, comm_);
  
  GNIDMatrix sendShNodeWrk(gsz[0],1);

  // mySharedData column length must be at least 
  // the number of tasks work data was sent to:
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "gsz1=" << gsz[0] << " gsz2=" << gsz[1]);
#endif
  mySharedData.resize(gsz[0],gsz[1]);


  // Last, fill send buffer with sorted work data:
  ipos = 0;
  sendShNodeWrk.set(-1);
  nworktmp.set(-1); nt = 0;
  for ( j=0; j<irWork.size(2); j++ ) { // cycle thru all tasks that sent work to this task
    nd = irWork.distinctcol_floor(j, mval, irowr, -1); // get distinct node ids for this column/task
    for ( i=0; i<nd; i++ ) { 
      nid = irWork(irowr[i],j);
      mult = irWork.multiplicity(nid, irow, icol); 
      if ( mult > 1 ) { // only keep the nodes with mult > 1
        // Keep only the unique column indices (task ids):
        rd   = icol.distinctrng(0,mult,1,vvals,ivals,nvals); 
        // Fill column of shared node mat:
        //    nodeid0 NTask0 TaskID00 TaskID01 ...; nodeid1 NTask1 TaskID10 TaskID11...:
        if ( !nworktmp.contains(nid,iwhich) ) {
          sendShNodeWrk(ipos,0) = nid;
          sendShNodeWrk(ipos+1,0) = rd;
          for ( k=0; k<rd; k++ ) sendShNodeWrk(ipos+k+2,0) = iRecvWorkTaskID[vvals[k]];
          ipos += rd + 2;
          nworktmp[nt] = nid;
          nt++;
        }
      }
    } 
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "sendShNodeWrk=" << sendShNodeWrk);
#endif

  mySharedData.set(-1);
  GCommDatatype dtype = T2GCDatatype<GNODEID>();

  // NOTE: Here, when sending back sorted data, we send back to 
  //       the tasks we received work data from, and we receive 
  //       the tasks that we sent work data to:
  // Set indirection buffers for send & receive (don't include this rank):
  GIBuffer nttmp(MAX(iSendWorkTaskID.size(),iRecvWorkTaskID.size()));

  // Find tasks to send back to (not including rank_):
  for ( j=0, ns=0; j<iRecvWorkTaskID.size(); j++ ) {
     itask = iRecvWorkTaskID[j];
     if ( itask != rank_ && !nttmp.containsn(itask,ns,iwhich) ) { 
       nttmp[ns] = iRecvWorkTaskID[j]; // task to send to
       ns++;
     } 
  }
  GIBuffer isbuff(ns);   // indirection buffer for sends
  GIBuffer itsend(nttmp.data(),ns);

  isbuff.set(0); // point to first col only in sendShNodeWrk

  // Find tasks to recv from (not including rank_):
  for ( j=0, nr=0; j<iSendWorkTaskID.size(); j++ ) {
     itask = iSendWorkTaskID[j];
     if ( itask != rank_ && !nttmp.containsn(itask,nr,iwhich) ) {
       nttmp[nr] = iSendWorkTaskID[j]; // task to send to
       nr++;
     } 
  }
  GIBuffer irbuff(nr);   // indirection buffer for recvs
  GIBuffer itrecv(nttmp.data(),nr);


  // Find indirection array for recv buffers:
  for ( j=0, nr=0; j<(nl=iSendWorkTaskID.size()); j++ ) {
     if ( iSendWorkTaskID[j] != rank_ ) {
       irbuff[nr] = nr;                 // buffer to receive into
       nr++;
     } 
     else {
      // Fill return matrix with local data:
      for ( i=0; i<sendShNodeWrk.size(1); i++ ) 
        mySharedData(i,nl-1) = sendShNodeWrk(i,0);
     }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "irbuff=" << irbuff);
  GPP(comm_, serr << "isbuff=" << isbuff);

  GPP(comm_, serr << "itrecv=" << itrecv);
  GPP(comm_, serr << "itsend=" << itsend);
  GPP(comm_, serr << "....................doing ASendRecv");
#endif

  // Send this work task's data back to tasks that we recvd data from,
  // and gather all works tasks' data in MySharedData:
  bret = GComm::ASendRecv(
                          mySharedData .data().data(),irbuff.size(),NULLPTR      ,mySharedData .size(1),
                          dtype, itrecv.data(), TRUE, 
                          sendShNodeWrk.data().data(),isbuff.size(),isbuff.data(),sendShNodeWrk.size(1),
                          dtype ,itsend.data(), comm_);

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << " mySharedData=" << mySharedData);
#endif

  if ( !bret ) {
    GPP(comm_,serr << "GComm::ASendRecv failed on resultant work data");
    exit(1);
  }

  if ( vvals != NULLPTR ) delete [] vvals;
  if ( ivals != NULLPTR ) delete [] ivals;
 

  return TRUE;

} // end of method doCommonNodeSort


//************************************************************************************
//************************************************************************************
// METHOD : extractOpData
// DESC   : Extracts data required to carrty out global operations from the
//          matrix data returned from doCommonNodeSort
// ARGS   : glob_index   : global index list input to Init
//          mySharedDat a: matrix data returned from doCommonNodeSort with each column record giving
//                        the global node id, and the tasks that share it. Each record/column
//                        is specific to a corresp task from which work data was received by
//                        this task. Which these are is not important, specifically.
//          worktask0, col 0:   nodeid0 NTask0 TaskID00 TaskID01 ...; nodeid1 NTask1 TaskID10 TaskID11...
//          worktask1, col 1:   nodeid0 NTask0 TaskID00 TaskID01 ...; nodeid1 NTask1 TaskID10 TaskID11...
//                ...
//
//            Note: 
//            iOpL2RTasks_  : list of task ids to send rows of iOpIndices to. Do not
//                            include this rank.
//            iOpL2RIndices_: 1 column for each task in iOpTasks, 
//                            containing rows of local indices into the glob_index array.
//                            These indices point to a single global node id, and
//                            this data is used to send the local data to other tasks that
//                            own it. Form is:
//                  col 0, task 0: index0 index1 index2....
//                  col 1, task 1: index0 index1 index2....
//                  col 2, task 2: index0 index1 index2....
//            nOpL2RIndices_: length of OpL2RIndices record to be sent to each task
//            iOpR2LIndices_: same dimensions as iOpL2RIndices, but with a 'depth'
//                            array at each matrix element giving the local representations
//                            of that global id.
//            nOpR2LMult_   : multiplicity of each value in iOpR2LIndices_ matrix.
//            iOpL2LIndices_: matrix of indices into local glob_index; each column
//                            represents a shared (global) index, and each row represents 
//                            an index to the same shared index. Used for only local operations
//                            (on this task only).
//            nOpL2LIndices_: for each column in iOpL2LIndices, the number
//                            of local indices representing the same shared node.
// RETURNS: 
//          TRUE on success; else FALSE
//************************************************************************************
template<typename T>
GBOOL GGFX<T>::extractOpData(GNIDBuffer &glob_index, GNIDMatrix &mySharedData)
{
  GString serr = "GGFX<T>::extractOpData: ";
  GBOOL   bind, bret;
  GSIZET  i, j, k, m;

  // Find unique list of task IDs contained within mySharedData, that we will have
  // to send to/recv from. Find also the local indices into glob_index that represent
  // 'global' shared nodes.

  // First, get sizes for indirection arrays:
  GINT       itask, nt ;
  GSIZET     index, index1, mult, nl, nnidmax;
  GNODEID    nnid;
  GIBuffer   itasktmp(nprocs_); itasktmp.set(-1);
  GSZBuffer  nindtmp (nprocs_); nindtmp.set(0);
  GNIDBuffer igltmp(glob_index.size()); igltmp.set(0);
  GSZBuffer  ngltmp(glob_index.size()); ngltmp.set(0);
  for ( j=0, nt=0, nl=0; j<mySharedData.size(2); j++ ) { // cyle over all work tasks
    if ( mySharedData(0,j)==-1 ) continue; // column contain no usable data
    i = 0; 
    while ( i<mySharedData.size(1) && mySharedData(i,j)>-1 ) { // parse jth record
      nnid = mySharedData(i,j);       // shared node id
      for ( k=0; k<mySharedData(i+1,j); k++ ) { // self-ref data gives # tasks sharing this id
        itask = mySharedData(i+k+2,j);// get taskd id for ith shared node
        bind = glob_index.contains(nnid, index1);
        if ( itask != rank_ ) {       // if remote task != this rank
          if ( bind ) {
            if ( !itasktmp.containsn(itask, nt, index) )  {
              itasktmp[nt] = itask;   // add remote task id to list of unique ids
              nindtmp [nt]++;         // update # global node ids for remote task, itask
              nt++;                   // update number of remote tasks to send to/recv from
            }
            else {
              nindtmp[index]++;       // if task in list, update number shared nodes (not multiply represetned)
            }
          } 
        } 
        else {                        // task == this rank; not remote
          if ( bind && !igltmp.containsn(nnid, nl) 
             && (mult = glob_index.multiplicity(nnid)) > 0 ) {  // # times this node id represented locally
            igltmp[nl]  = nnid;      // add global shared node to list of global/shared nodes
            ngltmp[nl] += mult;      // update # local node ids for this shared global id
            nl++;                    // update number of global id in local matrix
          }
        }
      } // end, for k loop
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    } // end, while j loop
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "nigltmp=" << nl << " igltmp=" << igltmp);
  GPP(comm_, serr << "nnidmax=" << nnidmax <<  " nt=" << nt);
#endif

  nnidmax = nindtmp.maxn(nt);       // max # global indices for a remote task
  iOpL2RTasks_  .resize(nt);
  nOpL2RIndices_.resize(nt);
  iOpL2RIndices_.resize(nnidmax,nt); // of size # unique shared nodes X  # remote tasks 

  iOpL2RIndices_.set(999999);

  // Resize send, recv, & local buffers for operations:
  GSIZET lmax[2], gmax[2];
  lmax[0] = nnidmax; lmax[1] = nt; 
  GComm::Allreduce(lmax, gmax, 2, T2GCDatatype<GSIZET>() , GC_OP_MAX, comm_);
          
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "gmax1=" << gmax[0] << " gmax2=" << gmax[1]);
    GPP(comm_,serr << "ngltmp.maxn(nl)=" << ngltmp.maxn(nl)  << " nl=" << nl);
#endif
  sendBuff_.resize(gmax[0],gmax[1]);
  recvBuff_.resize(gmax[0],gmax[1]);    

  iOpL2RTasks_  .set(itasktmp.data(),nt);    // task ids to send/recv to/from
  nOpL2RIndices_.set(nindtmp.data (),nt);    // # indices to send/rcv for each task

  iOpL2LIndices_.resize(ngltmp.maxn(nl),nl); // local indices pointing to shared global nodes
  iOpL2LIndices_.set(999999);
  nOpL2LIndices_.resize(nl); 
  nOpL2LIndices_.set(ngltmp.data(),nl);

  // Find indirection arrays for local indices of shared nodes that will
  // be sent & recvd from remote tasks, and for these that are entirely local:

  // ... indirection arrays (L2R = 'local to remote') for remote sends & receives:
  // 1 column for each task to send to, and rows refer to the indices:
  nindtmp.set(0);
  for ( j=0; j<mySharedData.size(2); j++ ) {
    if ( mySharedData(0,j)==-1 ) continue; // column contains no usable data
    i = 0;
    while ( i<mySharedData.size(1) && mySharedData(i,j)>-1 ) { // parse jth record
      nnid = mySharedData(i,j); // shared node id
      for ( k=0; k<mySharedData(i+1,j); k++ ) { // self-ref data gives # tasks sharing this id
        itask = mySharedData(i+k+2,j);
        if ( itask == rank_ ) continue;
        if ( glob_index.contains(nnid,index1) ) {
          if ( iOpL2RTasks_.contains(itask, index) ) {
            iOpL2RIndices_(nindtmp[index],index) = index1; // assign local index to global id
            nindtmp[index]++;
          } 
        }
      }
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    }
  }

  // Re-order iOpL2RIndices s.t. global ids indices refer to
  // are re-ordered smallest to largest. By basing ordering on 
  // global ids, we ensure that the sent and received data
  // align:
  GNIDBuffer iglsort(nOpL2RIndices_.max()); 
  GSZBuffer  isort, ilocal(nOpL2RIndices_.max());
  for ( j=0, m=0; j<iOpL2RIndices_.size(2); j++ ) { // loop over tasks
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) { // get list of global ids
      iglsort[i] = glob_index[iOpL2RIndices_(i,j)];
      ilocal [i] = iOpL2RIndices_(i,j);     // copy this task's local ids
    }
    iglsort.range(0,nOpL2RIndices_[j]-1);   // restrict range
    iglsort.sortincreasing(isort);         // sort global list
    iglsort.range_reset();                 // reset range
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) { // re-order local indices
      iOpL2RIndices_(i,j) = ilocal[isort[i]];
    }
  }
  

  // ...indirection arrays (R2L = 'remote to local'): This is a matrix of same
  // size as iOpL2RIndices, but with a 'depth' that indicates the multiplicity
  // of each global id from each task. It is used to combine remote data with
  // local data:
  GSIZET  maxmult=0;
  for ( j=0; j<iOpL2RIndices_.size(2); j++ ) { // loop over all unique indices to find max sizes
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) {
      nnid = glob_index[iOpL2RIndices_(i,j)];
      maxmult = MAX(maxmult,glob_index.multiplicity(nnid)); 
    }
  }
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "iOpL2RIndices_.size1=" << iOpL2RIndices_.size(1)
        << " iOpL2RIndices_.size2=" << iOpL2RIndices_.size(2));
#endif
  iOpR2LIndices_.resize(iOpL2RIndices_.size(1),iOpL2RIndices_.size(2)); 
  nOpR2LMult_.resize(iOpL2RIndices_.size(1),iOpL2RIndices_.size(2)); 
  iOpR2LIndices_.set(999999); 
  nOpR2LMult_.set(0);
  for ( j=0; j<iOpR2LIndices_.size(2); j++ ) { 
    for ( i=0; i<iOpR2LIndices_.size(1); i++ ) iOpR2LIndices_(i,j).resize(maxmult); 
  }

  GSIZET *iwhere=0;
  GSIZET  nw=0;
  for ( j=0, m=0; j<iOpR2LIndices_.size(2); j++ ) { // loop over global ids
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) { // for each unique global index, find local indices
      nnid = glob_index[iOpL2RIndices_(i,j)]; // get 'global' id
      mult = glob_index.multiplicity(nnid, iwhere, nw);
      for ( k=0; k<mult; k++ ) iOpR2LIndices_(i,j)[k] = iwhere[k]; 
      nOpR2LMult_(i,j) = mult;
      m++;
    }
  }

  // ... indirection arrays (L2L = 'local to local ') for local operations:
  // 1 column for each 'global' node, and rows refer to the local indices that
  // represent multiplicity of that node id:
//ngltmp.set(0);
  for ( j=0, bret=TRUE; j<mySharedData.size(2) && bret; j++ ) {
    if ( mySharedData(0,j)==-1 ) break; // this & following rows contain no data
    i = 0;
    while ( bret && i<mySharedData.size(1) && mySharedData(i,j)>-1 ) { // parse jth record
      nnid = mySharedData(i,j); // shared node id
      for ( k=0; k<mySharedData(i+1,j); k++ ) { // self-ref data gives # tasks sharing this id
        itask = mySharedData(i+k+2,j);
        if ( itask != rank_ ) continue;
        bind = igltmp.containsn(nnid, nl, index); // where global id sits in local tmp array
        if ( bind ) {
          mult  = glob_index.multiplicity(nnid, iwhere, nw); // # times this node id represented
          for ( m=0; m<mult; m++ ) iOpL2LIndices_(m,index) = iwhere[m];
        }
      }
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    } 
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "iOpL2RTasks  ="  << iOpL2RTasks_ );
  GPP(comm_,serr << "nOpL2RIndices="  << nOpL2RIndices_);
  GPP(comm_,serr << "iOpL2RIndices="  << iOpL2RIndices_);
  GPP(comm_,serr << "nOpR2LIndices="  << nOpR2LMult_);
  GPP(comm_,serr << "iOpR2LIndices="  << iOpR2LIndices_);
  GPP(comm_,serr << "nOpL2LIndices="  << nOpL2LIndices_);
  GPP(comm_,serr << "iOpL2LIndices="  << iOpL2LIndices_);
#endif


  if ( iwhere != NULLPTR ) delete [] iwhere;
 
  return bret; 

} // end, method extractOpData


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

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "iOpL2LIndices=" << iOpL2LIndices_);
  GPP(comm_, serr << "Doing first localGS...");
#endif

  GTimerStart("ggfx_doop");

  
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

  GTimerStart("ggfx_doop_exch");

  // Perform a exchange of field data:
  bret = dataExchange(u);

  GTimerStop("ggfx_doop_exch");

  if ( !bret ) {
    std::cout << serr << "dataExchange.failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "dataExchange done.");
  GPP(comm_,serr << " recvBuff=" << recvBuff_);
#endif
  bret = localGS(u, iOpL2RIndices_, nOpL2RIndices_, op, &recvBuff_);
  if ( !bret ) {
    std::cout << serr << "localGS (2) failed " << std::endl;
    exit(1);
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "done.");
#endif

  GTimerStop("ggfx_doop");

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

  GTimerStart("ggfx_doop");

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


  GTimerStart("ggfx_doop_exch");

  // Perform a exchange of field data:
  bret = dataExchange(u, iind);

  GTimerStop("ggfx_doop_exch");

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

  GTimerStop("ggfx_doop");

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
  GLLONG  i, j, k;

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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iOpR2LIndices_(i,j)[k]] += (*qop)(i,j);
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iOpR2LIndices_(i,j)[k]] *= (*qop)(i,j);
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iOpR2LIndices_(i,j)[k]] = MAX(qu[iOpR2LIndices_(i,j)[k]],(*qop)(i,j));
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iOpR2LIndices_(i,j)[k]] = MIN(qu[iOpR2LIndices_(i,j)[k]],(*qop)(i,j));
            }
          }
        }
      }
      break;
    case GGFX_OP_SMOOTH:
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iOpR2LIndices_(i,j)[k]] += (*qop)(i,j);
            }
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
  GLLONG  i, j, k;


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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iind[iOpR2LIndices_(i,j)[k]]] += (*qop)(i,j);
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iind[iOpR2LIndices_(i,j)[k]]] *= (*qop)(i,j);
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iind[iOpR2LIndices_(i,j)[k]]] *= MAX(qu[iind[iOpR2LIndices_(i,j)[k]]],(*qop)(i,j));
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iind[iOpR2LIndices_(i,j)[k]]] *= MIN(qu[iind[iOpR2LIndices_(i,j)[k]]],(*qop)(i,j));
            }
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
            for ( k=0; k<nOpR2LMult_(i,j); k++ ) {
               qu[iind[iOpR2LIndices_(i,j)[k]]] += (*qop)(i,j);
            }
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
  GString serr = "GGFX<T>::initMult: ";

  assert(bInit_ && "Operator not initialized");

  GTVector<T> mult(nglob_index_);

  mult = 1.0;

  // Do DSS sum to find multiplicity:
  doOp(mult, GGFX_OP_SUM);


  // Compute 1/mult:
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << " mult.size=" << mult.size());
#endif
  imult_.resize(mult.size());
  for ( GSIZET j=0; j<mult.size(); j++ ) {
    imult_[j] = 1.0/mult[j];
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,"GGFX<T>::initMult: mult=" << mult);
  GPP(comm_,"GGFX<T>::initMult: imult=" << imult_);
#endif

} // end of method initMult


