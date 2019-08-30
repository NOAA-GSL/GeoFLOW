//==================================================================================
// Module       : gbitblock
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a bit-block data type for manipulating bits in a block
//                of arbitrary size.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <math.h>
#include "gbitblock.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with total no. bits desired
// ARGS   : GINT NumOfBits: bit block size
// RETURNS: none
//**********************************************************************************
GBitBlock::GBitBlock(GINT  NumOfBits)
:
iBlock_                (NULL),
bitFieldLength_        (NumOfBits),
nBits_                 (BITSPERBYTE*sizeof(GUSHORT )),
numIntSubBlocks_       (0)
{
  GUSHORT  iOffset = NumOfBits%(nBits_) > 0 ? 1 : 0;

  numIntSubBlocks_ = (NumOfBits)/(nBits_) + iOffset;
  iBlock_ = new GUSHORT  [numIntSubBlocks_];
  memset(iBlock_,0,numIntSubBlocks_*sizeof(GUSHORT ));

} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GBitBlock::~GBitBlock()
{
  if ( iBlock_ != NULL ) 
  {
    delete [] iBlock_;
    iBlock_ = NULL;
  }
}


//**********************************************************************************
//**********************************************************************************
// Copy constructor method
GBitBlock::GBitBlock(const GBitBlock &a)
{
  GINT  i;

  // copy member data:
  bitFieldLength_ = a.bitFieldLength_;
  numIntSubBlocks_= a.numIntSubBlocks_;
  iBlock_ = new GUSHORT  [numIntSubBlocks_];
  for ( i=0; i<numIntSubBlocks_; i++ )
     *(iBlock_+i) =  *(a.iBlock_+i); 
   
} // end of copy constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : operator=
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : GBitBlock &
//**********************************************************************************
GBitBlock &GBitBlock::operator=(const GBitBlock &a)
{
  GINT  i;

  bitFieldLength_ = a.bitFieldLength_;
  numIntSubBlocks_= a.numIntSubBlocks_;
  if ( iBlock_ != NULL) delete [] iBlock_;
  iBlock_ = new GUSHORT  [numIntSubBlocks_];
  for ( i=0; i<numIntSubBlocks_; i++ )
     *(iBlock_+i) =  *(a.iBlock_+i);

  return *this;  

} // end of = operator


//**********************************************************************************
//**********************************************************************************
// METHOD     : () operator  
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : GUSHORT  bit value (0 or 1)
//**********************************************************************************

GUSHORT  GBitBlock::operator()(const GINT  j) 
{
  GINT  i=j, iOffset, ksubblock;

  // (position 0 is right-most end)

  // find subblock:
  if ( j<0 || j>bitFieldLength_-1 )
  {
    std::cout << "GBitBlock::() : invalid bit index" << std::endl;
    exit(1);
  }

  iOffset = (i+1)%(nBits_) > 0 ? 1 : 0;
  ksubblock = (i+1) / (nBits_) + iOffset - 1;
  
  return getbits_basic(*(iBlock_+ksubblock), i-ksubblock*nBits_, 1);

} // end of operator () 

//**********************************************************************************
//**********************************************************************************
// METHOD     : [] operator  
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : GUSHORT  bit value (0 or 1)
//**********************************************************************************
  
GUSHORT  GBitBlock::operator[](const GINT  j)
{ 
  GINT  i=j, iOffset, ksubblock;
  
  // (position 0 is right-most end)

  // find subblock:
  if ( j<0 || j>bitFieldLength_-1 )
  {
    std::cout << "GBitBlock::[] : invalid bit index" << std::endl;
    exit(1);
  }
  
  iOffset = (i+1)%(nBits_) > 0 ? 1 : 0;
  ksubblock = (i+1) / (nBits_) + iOffset - 1;
  
  return getbits_basic(*(iBlock_+ksubblock), i-ksubblock*nBits_, 1);

} // end of operator []


//**********************************************************************************
//**********************************************************************************
// METHOD     : operator<<
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : ostream &
//**********************************************************************************
std::ostream &operator<<(std::ostream &str, GBitBlock &a)
{

  str << "\n";

  GINT  i;
  str << " Bits: " << std::endl;
  for ( i=a.getBlockSize_InBits()-1; i>=0; i-- )
    str << a(i) << "  ";
  str << std::endl << std::endl;

  str << " GUSHORT  Blocks: " << std::endl;
  for ( i=0; i<a.getBlockSize_InBlocks(); i++ )
    str << *(a.getBlock()+i) << " ";
  str << std::endl << std::endl;

  return str;
} // end of stream operator


//**********************************************************************************
//**********************************************************************************
// METHOD     : setBits
// DESCRIPTION:
// ARGUMENTS  : j  : which bit
//              val: 0 or 1
// RETURNS    : none.
//**********************************************************************************
void GBitBlock::setBits(GINT  j, const GUSHORT val )
{
  GINT  i=j, iOffset, ksubblock;
  GUSHORT  yval=0, iret;
  
  if ( val != 0 ) yval = 1;

  // (position 0 is right--most end)

  // find subblock:
  if ( j<0 || j>bitFieldLength_-1 )
  {
    std::cout << "GBitBlock::setBits: invalid bit index" << std::endl;
    exit(1);
  }
//if ( j < 0 ) i = 0;
//if ( j > bitFieldLength_-1 ) i = bitFieldLength_-1;

  iOffset = (i+1)%(nBits_) > 0 ? 1 : 0;
  ksubblock = (i+1) / (nBits_) + iOffset - 1;

  iret  = setbits_basic(*(iBlock_+ksubblock), i-ksubblock*nBits_, 1, yval);
  *(iBlock_+ksubblock) = iret;

} // end of method setBits



//**********************************************************************************
//**********************************************************************************
// METHOD     : getbits_basic
// DESCRIPTION:
// DESCRIPTION: Basic getbits routine, privates; gets at sub-block level
// RETURNS    : bits for subblock of type GUSHORT 
//**********************************************************************************
GUSHORT  GBitBlock::getbits_basic(GUSHORT  x, int p, int n)
{
  GUSHORT  iret;

  iret = ( (x >> (p-n+1)) & ~(~0 << n) ); 

  return iret;

} // end of method getbits_basic


//**********************************************************************************
//**********************************************************************************
// METHOD     : setbits_basic
// DESCRIPTION: Basic setBits routine, privates; sets at sub-block level
// ARGUMENTS  :
// RETURNS    : subblock GUSHORT  value
//**********************************************************************************
GUSHORT  GBitBlock::setbits_basic(GUSHORT  x, int p, int n, GUSHORT  y)
{
 
  GUSHORT  ytmp = 0;
  GUSHORT  iret ;

  ytmp = getbits_basic(y,n-1,n);

  // (position 0 is right--most end)


  // zero n bits starting at p by building
  // mask and &'ing with x; then turn on
  // y's bits, left shifted to p:
  // zero'ing mask= ( (~0 << p+1) | ~(~0 << (p-n+1)) )

  iret =  ( ( (~0 << (p+1)) | ~(~0 << (p-n+1)) ) & x ) | ( ytmp << p );

  return iret;

} // end of method setbits_basic

//**********************************************************************************
//**********************************************************************************
// METHOD     : setBlock (1)
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GBitBlock::setBlock(GUSHORT   *newblock, GINT  num_subblocks)
{
  if ( newblock == NULL || num_subblocks < 1 ) return FALSE;
 
  GINT  n, i;

  n = num_subblocks;
   
  delete [] iBlock_;
  iBlock_ = NULL;

  iBlock_ = new GUSHORT  [n];
  if ( iBlock_ == NULL ) return FALSE;

  for ( i=0; i<n; i++ )
    iBlock_[i] = newblock[i];

  numIntSubBlocks_ = n;
  bitFieldLength_  = numIntSubBlocks_ * (nBits_);

  return TRUE;
} // end of method setBlock (1)

  
//**********************************************************************************
//**********************************************************************************
// METHOD     : setBlock (2)
// DESCRIPTION:
// ARGUMENTS  :
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GBitBlock::setBlock(GTVector<GUSHORT >  *newblock, GINT  num_subblocks)
{
  if ( newblock == NULL || num_subblocks < 1 ) return FALSE;
  
  GINT  n, i;

  n = MIN(num_subblocks, newblock->size());
    
  delete [] iBlock_;
  iBlock_ = NULL;

  iBlock_ = new GUSHORT  [n];
  if ( iBlock_ == NULL ) return FALSE;

  for ( i=0; i<n; i++ )
    iBlock_[i] = (*newblock)[i];
  
  numIntSubBlocks_ = n;
  bitFieldLength_  = numIntSubBlocks_ * (nBits_);

  return TRUE;
} // end of method setBlock (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : setBlock (3)
// DESCRIPTION: 
// Take bitspersubblock right-most bits of each of num_subblocks sublocks 
// of newblock (each sublock is size_subblocks GUSHORT s long), and contatenate
// them into a single block.
// ARGUMENTS  : 
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GBitBlock::setBlock(GUSHORT   *newblock, GINT  num_subblocks, GINT  size_subblocks, GINT  bitspersubblock)
{
  if ( newblock == NULL ) return FALSE;
  if ( bitspersubblock > size_subblocks*nBits_ ) return FALSE;
 
  GUSHORT  iret, yval=1;
  GINT  i,j, k, ibit, iOffset, ksubblock;

  delete [] iBlock_;
  iBlock_ = NULL;
  
  bitFieldLength_ = num_subblocks * bitspersubblock;
  

  iOffset = bitFieldLength_%(nBits_) > 0 ? 1 : 0;
  numIntSubBlocks_ = (bitFieldLength_)/(nBits_) + iOffset;
  iBlock_ = new GUSHORT  [numIntSubBlocks_];
  if ( iBlock_ == NULL ) return FALSE;
 
  i = 0;
  for ( k=0; k<num_subblocks*size_subblocks; k+=size_subblocks )
  {
    for ( j=0; j<bitspersubblock; j++ )
    {
      ksubblock = i / nBits_;
      yval      = getbits_basic(*(newblock+k+j/nBits_), j%nBits_, 1);
      ibit      = i - ksubblock*nBits_;
      iret      = setbits_basic(*(iBlock_+ksubblock), ibit, 1, yval);
      *(iBlock_+ksubblock) = iret;
      i++;
    }
  }

  return TRUE;
} // end of method setBlock (3)


//**********************************************************************************
//**********************************************************************************
// METHOD     : setBlock (4) 
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : TRUE on success; else FALSE
//**********************************************************************************
GBOOL GBitBlock::setBlock(GUSHORT  *hi_block, GUSHORT  *lo_block, GINT  nbsize)
{
  if ( hi_block == NULL || lo_block == NULL ) return FALSE;

  GINT  n, i;

  n = 2*nbsize;

  delete [] iBlock_;
  iBlock_ = NULL;

  iBlock_ = new GUSHORT  [n];
  if ( iBlock_ == NULL ) return FALSE;

  for ( i=0; i<n/2; i++ ) {
    iBlock_    [i] = lo_block[i];
    iBlock_[i+n/2] = hi_block[i];
  }

  numIntSubBlocks_ = n;
  bitFieldLength_  = numIntSubBlocks_ * (nBits_);

  return TRUE;
} // end of method setBlock (4)

//**********************************************************************************
//**********************************************************************************
// METHOD     : hiWord
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : none.
//**********************************************************************************
void *GBitBlock::hiWord()
{
  if ( numIntSubBlocks_%2 != 0 ) return NULL;
  return (iBlock_+numIntSubBlocks_/2);
}


//**********************************************************************************
//**********************************************************************************
// METHOD     : loWord
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : none.
//**********************************************************************************
void *GBitBlock::loWord()
{
  if ( numIntSubBlocks_%2 != 0 ) return NULL;
  return (iBlock_);
} // end, loWord method

//**********************************************************************************
//**********************************************************************************
// METHOD     : getBlock 
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : value of GUSHORT  sublock
//**********************************************************************************
GUSHORT  *GBitBlock::getBlock()
{
  return iBlock_;
} // end of method getBlock


//**********************************************************************************
//**********************************************************************************
// METHOD     : getBlockSize_InBits_
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : GINT  blocksize (num bits)
//**********************************************************************************
GINT  GBitBlock::getBlockSize_InBits()
{
 
  return bitFieldLength_;

} // end of method getBlockSize_InBits


//**********************************************************************************
//**********************************************************************************
// METHOD     : getBlockSize_InBlocks
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : GINT  number of GUSHORT  subblocks
//**********************************************************************************
GINT  GBitBlock::getBlockSize_InBlocks()
{
 
  return numIntSubBlocks_;

} // end of method getBlockSize_InBlocks


//**********************************************************************************
//**********************************************************************************
// METHOD     : reset
// DESCRIPTION: 
// ARGUMENTS  : 
// RETURNS    : none.
//**********************************************************************************
void GBitBlock::reset()
{
  memset(iBlock_,0,numIntSubBlocks_*sizeof(GUSHORT ));
} // end of method reset


//**********************************************************************************
//**********************************************************************************
// METHOD     : transferBits
// DESCRIPTION: Transfer bit block to the nelems in_elements of size elem_size. 
// ARGUMENTS  : 
// RETURNS    : none.
//**********************************************************************************
void GBitBlock::transferBits(void *in_elems, GINT  nelems, size_t elem_size)
{
#if 1
  memset(in_elems, 0, nelems*elem_size);
  memcpy(in_elems, iBlock_, MIN(numIntSubBlocks_*sizeof(GUSHORT),nelems*elem_size));
#else
  for ( GINT j=0; j<MIN(numIntSubBlocks_*sizeof(GUSHORT),nelems*elem_size); j++ ) {
    memcpy( (GUSHORT*)in_elems+j*sizeof(GUSHORT), iBlock_+j , sizeof(GUSHORT));
  }
#endif

} // end of method transferBits


