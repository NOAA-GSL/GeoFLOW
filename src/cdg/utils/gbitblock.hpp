//==================================================================================
// Module       : gbitblock
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a bit-block data type for manipulating bits in a block
//                of arbitrary size.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(GBITBLOCK_HPP)
#define GBITBLOCK_HPP
#include "gtypes.h"
#include <string.h>
#include "gtvector.hpp"


class GBitBlock 
{
public:
                  GBitBlock(GINT  Size);
                 ~GBitBlock();
                  GBitBlock(const GBitBlock &a);
GBitBlock         &operator=(const GBitBlock &);
GUSHORT           operator()(const GINT  bit_posn);
GUSHORT           operator[](const GINT  bit_posn);
friend std::ostream&     operator<<(std::ostream&, GBitBlock&);

void              setBits(GINT  bit_index, const GUSHORT  yTo);
//void            setBits(const GUSHORT  yTo);
GBOOL             setBlock(GUSHORT  *new_block, GINT  num_subblocks);
GBOOL             setBlock(GUSHORT  *new_block, GINT  num_subblocks, GINT  size_subblock, GINT  desired_bitspersubblock);
GBOOL             setBlock(GTVector<GUSHORT > *new_block, GINT  num_subblocks);
GBOOL             setBlock(GUSHORT  *hi_block, GUSHORT  *lo_block, GINT  n);
GUSHORT          *getBlock();
void             *loWord();
void             *hiWord();
GINT              getBlockSize_InBits();
GINT              getBlockSize_InBlocks();
void              reset();
void              transferBits(void *in_elems, GINT  nelems, size_t elem_size);


private:
// Private methods:
GUSHORT          getbits_basic(GUSHORT , int posn, int numbits);
GUSHORT          setbits_basic(GUSHORT , int posn, int numbits, GUSHORT  y);

// Private data:
GUSHORT          *iBlock_;
GINT              bitFieldLength_;
GINT              nBits_;   
GINT              numIntSubBlocks_;

};
#endif

