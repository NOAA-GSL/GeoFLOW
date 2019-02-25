//==================================================================================
// Module       : ggfx.cpp
// Date         : 5/9/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a geometry--free global exchange (GeoFLOW Geometry-Free eXchange)
//                operator
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <math.h>
#include <type_traits>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"
 

//************************************************************************************
//************************************************************************************
// METHOD : Constructor
// DESC   : 
// ARGS   : GD_COMM object, defaults to GD_COMM_WORLD
// RETURNS: GGFX
//************************************************************************************
GGFX::GGFX(GC_COMM icomm)
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
GGFX::~GGFX()
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
GGFX::GGFX(const GGFX &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GGFX
// RETURNS: GGFX
//************************************************************************************
GGFX  &GGFX::operator=(const GGFX &a)
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
GBOOL GGFX::init(GNIDBuffer &glob_index)
{
  GBOOL      bret;
  GString   serr = "GGFX::init: ";

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
  bret = initSort(glob_index, iOpL2RTasks_, iOpL2RIndices_, nOpL2RIndices_, 
                  iOpR2LIndices_, nOpR2LIndices_, iOpL2LIndices_, nOpL2LIndices_);
  if ( !bret ) {
    GPP(comm_,serr << "initSort failed");
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
    GPP(comm_,serr << "initSort done.");
#endif

  bInit_ = TRUE;
  return bInit_;

} // end of method init


//************************************************************************************
//************************************************************************************
// METHOD: initSort
// DESC  : Performs sorting for initialization of global gather-scatter operation.
// ARGS   : glob_index: global index list input to Init
//          iOpL2RTasks     : returned; list of task ids to send rows of iOpIndices to;
//                           receive from
//          iOpL2RIndices   : matrix returned , 1 row for each task in iOpTasks, 
//                           containing columns of local indices into the glob_index array:
//                row 0, task 0: index0 index1 index2....
//                row 1, task 1: index0 index1 index2....
//                row 2, task 2: index0 index1 index2....
//                           This data is used to send data.
//          nOpL2RIndices : returned, number of OpSRIndices for each task
//          iOpL2LIndices: returned, matrix of indices into local glob_index; each column
//                         represents a shared (global) index, and each row represents 
//                         the index to the same shared index 
//          nOpL2LIndices: returned, for each column in iOpL2LIndices, the number
//                         of local indices representing the same shared node.
//                ...
// RETURNS: 
//          TRUE on success; else FALSE
//************************************************************************************
GBOOL GGFX::initSort(GNIDBuffer  &glob_index, GIBuffer &iOpL2RTasks, 
                     GSZMatrix &iOpL2RIndices, GSZBuffer &nOpL2RIndices, 
                     GSZMatrix &iOpR2LIndices, GSZBuffer &nOpR2LIndices, 
                     GSZMatrix &iOpL2LIndices, GSZBuffer &nOpL2LIndices)
{
  GString serr = "GGFX::initSort: ";
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
  bret = extractOpData(glob_index, mySharedData, iOpL2RTasks, iOpL2RIndices, nOpL2RIndices, 
                       iOpR2LIndices, nOpR2LIndices, iOpL2LIndices, nOpL2LIndices);
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
GBOOL GGFX::binSort(GNIDBuffer &nodelist, GIMatrix &gBinMat, 
                    GINT &numlocfilledbins, GINT &max_numfilledbins, GINT &max_numlocbinmem, 
                    GNIDMatrix &gBinBdy, GNIDBuffer &locWork)
{
  GString  serr = "GGFX::binSort: ";
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
  if ( rank_ == 0 ) 
    std::cout << serr << "gBinBdy=" << gBinBdy << std::endl;
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
GBOOL GGFX::createWorkBuffs(GIMatrix &gBinMat, GINT max_numlocbinmem, GIBuffer &iRecvWorkTaskID, 
                            GNIDMatrix &irWork, GIBuffer &iSendWorkTaskID, GNIDMatrix &isWork)
{
  GString  serr = "GGFX::createWorkBuffs: ";
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
GBOOL GGFX::binWorkFill(GNIDBuffer &nodelist, GNIDMatrix &gBinBdy, GNIDMatrix  &sWork, GIBuffer &iSendWorkTaskID)
{
  GString serr = "GGFX::binWorkFill: ";
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
GBOOL GGFX::doSendRecvWork( GNIDBuffer &glob_index      , GNIDMatrix &gBinBdy, 
                            GIBuffer   &iRecvWorkTaskID , GNIDMatrix &irWork, 
                            GIBuffer   &iSendWorkTaskID , GNIDMatrix &isWork)
{
  GString serr = "GGFX::doSendRecvWork: ";
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
  // recv buffer.
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
      itsend[nsend] = iSendWorkTaskID[j];
      ibsend[nsend] = j;
      nsend++;
    }
    else {
      for ( i=0; i<isWork.size(1); i++ ) irWork(i,maxrecv-1) = isWork(i,j);
    }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "nsend=" << nsend << " ibsend=" << ibsend);
  GPP(comm_,serr << "itsend=" << itsend);
  GPP(comm_,serr << "nrecv=" << nrecv << " itrecv=" << iRecvWorkTaskID );
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
  GPP(comm_, serr << "done. irWork=" << irWork);
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
GBOOL GGFX::doCommonNodeSort(GNIDBuffer &glob_index, GNIDMatrix &irWork, 
                             GIBuffer &iRecvWorkTaskID, GIBuffer &iSendWorkTaskID,
                             GNIDMatrix &mySharedData)
{
  GString    serr = "GGFX::doCommonNodeSort: ";
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
//          iOpL2RTasks   : returned; list of task ids to send rows of iOpIndices to. Do not
//                          include this rank.
//          iOpL2RIndices : returned; 1 column for each task in iOpTasks, 
//                          containing rows of local indices into the glob_index array.
//                          These indices point to a single global node id, and
//                          this data is used to send the local data to other tasks that
//                          own it. Form is:
//                col 0, task 0: index0 index1 index2....
//                col 1, task 1: index0 index1 index2....
//                col 2, task 2: index0 index1 index2....
//          nOpL2RIndices : returned; length of OpIndices record to be sent to each task
//          iOpL2LIndices : returned; matrix of indices into local glob_index; each column
//                          represents a shared (global) index, and each row represents 
//                          an index to the same shared index. Used for only local operations
//                          (on this task only).
//          nOpL2LIndices : returned; for each column in iOpL2LIndices, the number
//                          of local indices representing the same shared node.
// RETURNS: 
//          TRUE on success; else FALSE
//************************************************************************************
GBOOL GGFX::extractOpData(GNIDBuffer &glob_index, GNIDMatrix &mySharedData, GIBuffer &iOpL2RTasks, 
                          GSZMatrix &iOpL2RIndices, GSZBuffer &nOpL2RIndices,
                          GSZMatrix &iOpR2LIndices, GSZBuffer &nOpR2LIndices,
                          GSZMatrix &iOpL2LIndices, GSZBuffer &nOpL2LIndices)
{
  GString serr = "GGFX::extractOpData: ";
  GBOOL   bind, bret;
  GSIZET  i, j, k, m;

  // Find unique list of task IDs contained within mySharedData, that we will have
  // to send to/recv from. Find also the local indices into glob_index that represent
  // 'global shared nodes.

  // First, get sizes for indirection arrays:
  GINT      itask, nt ;
  GSIZET    index, index1, mult, nl, nnidmax;
  GNODEID   nnid;
  GIBuffer  itasktmp(nprocs_); itasktmp.set(-1);
  GSZBuffer nindtmp (nprocs_); nindtmp.set(0);
  GSZBuffer igltmp(glob_index.size()); igltmp.set(0);
  GSZBuffer ngltmp(glob_index.size()); ngltmp.set(0);
  for ( j=0, nt=0, nl=0; j<mySharedData.size(2); j++ ) {
    if ( mySharedData(0,j)==-1 ) continue; // column contain no usable data
    i = 0; 
    while ( i<mySharedData.size(1) && mySharedData(i,j)>-1 ) { // parse jth record
      nnid = mySharedData(i,j);        // shared node id
      for ( k=0; k<mySharedData(i+1,j); k++ ) { // self-ref data gives # tasks sharing this id
        itask = mySharedData(i+k+2,j); // get taskd id for ith shared node
        bind = glob_index.contains(nnid, index1);
        if ( itask != rank_ ) {        // if remote task != this rank
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
        else {                      // task == this rank; not remote
          if ( bind && !igltmp.containsn(nnid, nl) 
             && (mult = glob_index.multiplicity(nnid)) > 1 ) {  // # times this node id represented locally
            igltmp[nl] = nnid;      // add global shared node to list of global/shared nodes
            ngltmp [nl]+=mult;      // update # local node ids for this shared global id
            nl++;                   // update number of global id in local matrix
          }
        }
      } // end, for k loop
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    } // end, while j loop
  }


  nnidmax = nindtmp.maxn(nt);       // max # global indices for a remote task
  iOpL2RTasks  .resize(nt);
  nOpL2RIndices.resize(nt);
  iOpL2RIndices.resize(nnidmax,nt); // of size # unique remote tasks X max # unique shared nodes

  iOpL2RIndices.set(999999);

  // Resize send, recv, & local buffers for operations:
  GSIZET lmax[2], gmax[2];
  lmax[0] = nnidmax; lmax[1] = nt; 
  GComm::Allreduce(lmax, gmax, 2, T2GCDatatype<GSIZET>() , GC_OP_MAX, comm_);
          
  sendBuff_.resize(gmax[0],gmax[1]);
  recvBuff_.resize(gmax[0],gmax[1]);    

  iOpL2RTasks  .set(itasktmp.data(),nt);    // task ids to send/recv to/from
  nOpL2RIndices.set(nindtmp.data (),nt);    // # indices to send/rcv for each task

  iOpL2LIndices.resize(ngltmp.maxn(nl),nl); // local indices pointing to shared global nodes
  iOpL2LIndices.set(999999);
  nOpL2LIndices.resize(nl); 
  nOpL2LIndices.set(ngltmp.data(),nl);

  // Find indirection arrays for local indices of shared nodes that will
  // be sent & recvd from remote tasks, and for these that are entirely local:

  // ... indirection arrays (G2R = 'global to remote') for remote sends & receives:
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
          if ( iOpL2RTasks.contains(itask, index) ) {
            iOpL2RIndices(nindtmp[index],index) = index1; // assign local index to global id
            nindtmp[index]++;
          } 
        }
      }
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    }
  }


  // ...indirection arrays (R2L = 'remote to local'): apply remote data to local
  // indices: 1 column for each global index, and rows indicate indices of other 
  // representations of that global index; sort of an inversion of iOpL2RIndices:
  GSIZET  maxmult=0;
  for ( j=0; j<iOpL2RIndices.size(2); j++ ) { // loop over all unique indices to find max sizes
    for ( i=0; i<nOpL2RIndices[j]; i++ ) {
      nnid = glob_index[iOpL2RIndices(i,j)];
      maxmult = MAX(maxmult,glob_index.multiplicity(nnid)); 
    }
  }
  iOpR2LIndices.resize(maxmult,nOpL2RIndices.sum()); iOpR2LIndices.set(999999); 
  nOpR2LIndices.resize(iOpR2LIndices.size(2)); 

  GSIZET *iwhere=0;
  GSIZET  nw=0;
  for ( j=0, m=0; j<iOpL2RIndices.size(2); j++ ) { // loop over unique global indices
    for ( i=0; i<nOpL2RIndices[j]; i++ ) { // for each unique global index, find local indices
      nnid = glob_index[iOpL2RIndices(i,j)];
      mult = glob_index.multiplicity(nnid, iwhere, nw);
      for ( k=0; k<mult; k++ ) iOpR2LIndices(k,m) = iwhere[k]; 
      nOpR2LIndices[m] = mult;
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
          for ( m=0; m<mult; m++ ) iOpL2LIndices(m,index) = iwhere[m];
        }
      }
      i += mySharedData(i+1,j) + 2; // skip to next global node id
    } 
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "iOpL2RTasks  ="  << iOpL2RTasks );
  GPP(comm_,serr << "nOpL2RIndices="  << nOpL2RIndices);
  GPP(comm_,serr << "iOpL2RIndices="  << iOpL2RIndices);
  GPP(comm_,serr << "nOpR2LIndices="  << nOpR2LIndices);
  GPP(comm_,serr << "iOpR2LIndices="  << iOpR2LIndices);
  GPP(comm_,serr << "nOpL2LIndices="  << nOpL2LIndices);
  GPP(comm_,serr << "iOpL2LIndices="  << iOpL2LIndices);
#endif


  if ( iwhere != NULLPTR ) delete [] iwhere;
 
  return bret; 

} // end, method extractOpData


//************************************************************************************
//************************************************************************************
// METHOD     : initMult
// DESCRIPTION: initialize for smoothing
// ARGUMENTS  : 
// RETURNS    : 
//************************************************************************************
void GGFX::initMult()
{
  GFVector mult(nglob_index_);

  mult = 1.0;

  // Do DSS sum to find multiplicity:
  doOp<GFLOAT>(mult, GGFX_OP_SUM);

  // Compute 1/mult:
  imult_.resize(mult.size());
  for ( GSIZET j=0; j<imult_.size(); j++ ) {
    imult_[j] = 1.0/mult[j];
  }

} // end of method initMult



//************************************************************************************
//************************************************************************************
// METHOD : getTimes
// DESC   : Retrieves timer results
// ARGS   : 
// RETURNS: 
//************************************************************************************
GDOUBLE GGFX::getTimes(GGFX_Timer_t type)
{
    return timer_data_[type];

} // end of method getTimes


//************************************************************************************
//************************************************************************************
// METHOD     : resetTimes
// DESCRIPTION: Reset timing variables
// ARGUMENTS  : 
// RETURNS    : 
//************************************************************************************
void GGFX::resetTimes()
{
  timer_data_ = 0.0;

} // end of method resetTimes

