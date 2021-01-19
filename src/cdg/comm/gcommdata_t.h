//==================================================================================
// Module       : gcommdata_t.hpp
// Date         : 3/1/18 (DLR)
// Description  : Provides global data types for GeoFLOW comm code. These types
//                correspond in kind with the GeoFLOW scientific processing portion,
//                as contained in gtypes.h
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(GCOMMDATA_T_H)
#define GCOMMDATA_T_H

#include "gtypes.h"
#include <type_traits>
#include <cassert>


//--------------------------------------------------------------------------------
// Comm Datatypes
//--------------------------------------------------------------------------------
#if defined(GEOFLOW_USE_MPI)
#include <stdint.h>
#include <limits.h>
#include "mpi.h"

#define GC_COMM         MPI_Comm
#define GC_COMM_WORLD   MPI_COMM_WORLD

#define GC_OP_TYPE      MPI_Op
#define GCommDatatype   MPI_Datatype
#define AGINT           MPI_Aint

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "Invalid SIZE_MAX"
#endif

//  Just take the MPI defs...
//   NOTE: Make sure there is a 1-to-1 correspondence
//         between these types, and GD_**** types in gtypes.h:
const GCommDatatype GC_DATATYPE[] = 
                        {MPI_BYTE     , MPI_C_BOOL            , MPI_UNSIGNED_CHAR, MPI_BYTE              ,
                         MPI_SHORT    , MPI_UNSIGNED_SHORT    , MPI_INT          , MPI_UNSIGNED          ,
                         MPI_LONG     , MPI_UNSIGNED_LONG     , MPI_LONG_LONG    , MPI_INT               ,
                         MPI_LONG     , MPI_UNSIGNED_LONG_LONG, MPI_LONG_LONG    , my_MPI_SIZE_T         ,
                         MPI_LONG_LONG, MPI_FLOAT             , MPI_DOUBLE       , MPI_LONG_DOUBLE       };

// Make sure there is a 1-1 corresp between GC_Optype and GC_OP
// indices in gtypes.h:
const GC_OP_TYPE GC_Optype[] =
                        {MPI_MAX   ,MPI_MIN   ,MPI_SUM    ,MPI_PROD  ,
                         MPI_LAND  ,MPI_BAND  ,MPI_LOR    ,MPI_BOR   ,
                         MPI_LXOR  ,MPI_BXOR};
#else
//  Else, define the comm datatypes as GD_DATATYPE enum types:
#define GC_COMM         GINT
#define GC_COMM_WORLD   1

#define GC_OP_TYPE      GINT
#define GCommDatatype   GINT
#define AGINT           GINT

const GCommDatatype GC_DATATYPE[] = 
                        {GD_GBYTE  ,GD_GBOOL  ,GD_GCHAR   ,GD_GUCHAR ,
                         GD_GSHORT ,GD_GUSHORT,GD_GINT    ,GD_GUINT  ,
                         GD_GLONG  ,GD_GULONG ,GD_GLLONG  ,GD_GWORD  ,
                         GD_GDWORD ,GD_GKEY   ,GD_GNODEID ,GD_GSIZET ,
                         GD_GFPOS  ,GD_GFLOAT ,GD_GDOUBLE ,GD_GQUAD  };

// Make sure there is a 1-1 corresp between GC_Optype and GC_OP
// indices in gtypes.h:
const GC_OP_TYPE GC_Optype[] =
                        {GC_OP_MAX ,GC_OP_MIN  ,GC_OP_SUM  ,GC_OP_PROD,
                         GC_OP_LAND,GC_OP__BAND,GC_OP_LOR  ,GC_OP_BOR ,
                         GC_OP_LXOR,GC_OP_BXOR};
                         
#endif


#if !defined(T2GCDATATYPE_DEFC)
#define T2GCDATATYPE_DEFC
//************************************************************************************
//************************************************************************************
// METHOD     : T2GCDatatype
// DESCRIPTION: Convert template type parameter to GCommDatatype at compile time
// ARGUMENTS  : 
// RETURNS    : 
//************************************************************************************
template<typename T> 
GCommDatatype T2GCDatatype() {
    GCommDatatype ctype; 
    if      ( std::is_same<T,  GBYTE>::value ) {
      ctype = GC_DATATYPE[GD_GBYTE] ; }
    else if ( std::is_same<T, GBOOL>::value ) {
      ctype = GC_DATATYPE[GD_GBOOL]; }
    else if ( std::is_same<T, GCHAR>::value ) {
      ctype = GC_DATATYPE[GD_GCHAR]; }
    else if ( std::is_same<T, GUCHAR>::value ) {
      ctype = GC_DATATYPE[GD_GUCHAR]; }
    else if ( std::is_same<T, GSHORT>::value ) {
      ctype = GC_DATATYPE[GD_GSHORT]; }
    else if ( std::is_same<T, GUSHORT>::value ) {
      ctype = GC_DATATYPE[GD_GUSHORT]; }
    else if ( std::is_same<T,   GINT>::value ) {
      ctype = GC_DATATYPE[GD_GINT] ; }
    else if ( std::is_same<T,  GUINT>::value ) {
      ctype = GC_DATATYPE[GD_GUINT]; }
    else if ( std::is_same<T,  GLONG>::value ) {
      ctype = GC_DATATYPE[GD_GLONG]; }
    else if ( std::is_same<T, GULONG>::value ) {
      ctype = GC_DATATYPE[GD_GULONG]; }
    else if ( std::is_same<T, GLLONG>::value ) {
      ctype = GC_DATATYPE[GD_GLLONG]; }
    else if ( std::is_same<T,  GWORD>::value ) {
      ctype = GC_DATATYPE[GD_GWORD]; }
    else if ( std::is_same<T, GDWORD>::value ) {
      ctype = GC_DATATYPE[GD_GDWORD]; }
    else if ( std::is_same<T,   GKEY>::value ) {
      ctype = GC_DATATYPE[GD_GKEY]; }
    else if ( std::is_same<T,GNODEID>::value ) {
      ctype = GC_DATATYPE[GD_GNODEID]; }
    else if ( std::is_same<T, GSIZET>::value ) {
      ctype = GC_DATATYPE[GD_GSIZET]; }
    else if ( std::is_same<T, GFPOS>::value ) {
      ctype = GC_DATATYPE[GD_GFPOS]; }
    else if ( std::is_same<T, GFLOAT>::value ) {
      ctype = GC_DATATYPE[GD_GFLOAT]; }
    else if ( std::is_same<T,GDOUBLE>::value ) {
      ctype = GC_DATATYPE[GD_GDOUBLE]; }
    else if ( std::is_same<T,  GQUAD>::value ) {
      ctype = GC_DATATYPE[GD_GQUAD];   }
    else {
      assert(FALSE);
    }
    return ctype;
  };  // end, T2GCDatatype

#endif

#if !defined(T2GDATATYPE_DEF)
#define T2GDATATYPE_DEF
//************************************************************************************
//************************************************************************************
// METHOD     : T2GDatatype
// DESCRIPTION: Convert template type parameter to GD_DATATYPE at compile time
// ARGUMENTS  : 
// RETURNS    : 
//************************************************************************************
template<typename T> 
GD_DATATYPE T2GDatatype() {
    GD_DATATYPE gtype; 
    if      ( std::is_same<T,  GBYTE>::value ) {
      gtype = GD_GBYTE ; }
    else if ( std::is_same<T, GBOOL>::value ) {
      gtype = GD_GBOOL; }
    else if ( std::is_same<T, GCHAR>::value ) {
      gtype = GD_GCHAR; }
    else if ( std::is_same<T, GUCHAR>::value ) {
      gtype = GD_GUCHAR; }
    else if ( std::is_same<T, GSHORT>::value ) {
      gtype = GD_GSHORT; }
    else if ( std::is_same<T, GUSHORT>::value ) {
      gtype = GD_GUSHORT; }
    else if ( std::is_same<T,   GINT>::value ) {
      gtype = GD_GINT ; }
    else if ( std::is_same<T,  GUINT>::value ) {
      gtype = GD_GUINT; }
    else if ( std::is_same<T,  GLONG>::value ) {
      gtype = GD_GLONG; }
    else if ( std::is_same<T, GULONG>::value ) {
      gtype = GD_GULONG; }
    else if ( std::is_same<T, GLLONG>::value ) {
      gtype = GD_GLLONG; }
    else if ( std::is_same<T,  GWORD>::value ) {
      gtype = GD_GWORD; }
    else if ( std::is_same<T, GDWORD>::value ) {
      gtype = GD_GDWORD; }
    else if ( std::is_same<T,   GKEY>::value ) {
      gtype = GD_GKEY; }
    else if ( std::is_same<T,GNODEID>::value ) {
      gtype = GD_GNODEID; }
    else if ( std::is_same<T, GSIZET>::value ) {
      gtype = GD_GSIZET; }
    else if ( std::is_same<T, GFPOS>::value ) {
      gtype = GD_GFPOS; }
    else if ( std::is_same<T, GFLOAT>::value ) {
      gtype = GD_GFLOAT; }
    else if ( std::is_same<T,GDOUBLE>::value ) {
      gtype = GD_GDOUBLE; }
    else if ( std::is_same<T,  GQUAD>::value ) {
      gtype = GD_GQUAD;   }
    else {
      assert(FALSE);
    }
    return gtype;
  };  // end, T2GDatatype

#endif

// GCMHandle: for calls to underlying exchange layer (in gcomm)
#if !defined(GCMHANDLE_DEF)
#define GCMHANDLE_DEF
class GCMHandle {
public :
GINT   bowned_; // Does this instance own its data?
GINT   nposts_; // Number posts = size of mhandle_ array
void  *mhandle_;
GCMHandle() { nposts_ = 0; mhandle_ = NULLPTR; bowned_ = 1; }
~GCMHandle() { 
#if defined(GEOFLOW_USE_MPI)
if ( bowned_ && mhandle_!=NULLPTR ) delete [] (MPI_Request*)mhandle_; 
#endif
}
#if 0
void operator=(GCMHandle &ch) { 
  GINT i;
  nposts_=ch.nposts_; 
  if ( mhandle_!=NULLPTR ) delete [] mhandle_; mhandle_ = NULLPTR;
  if ( ch.mhandle_ != NULLPTR ) { 
    for ( i=0; i<nposts_; i++ ) mhandle_[i] = ch.mhandle_[i];
  }
}  
#endif
friend std::ostream&  operator<<(std::ostream&os, GCMHandle &h) {
  os   << "nposts    =" << h.nposts_           << std::endl
       << "mhandle   =" << (int*) (h.mhandle_) << std::endl;
  return os;
}
};
#endif



#endif
