//==================================================================================
// Module       : gtypes.hpp
// Date         : 1/1/18 (DLR)
// Description  : Basic types, defs
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GTYPES_HPP)
#define _GTYPES_HPP

#include <cstddef>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "configure.hpp"



// Following is a list of preprocessor variables that may be set:
// _G_AUTO_CREATE_DEV : Auto-copy/create classes on device
// _G_AUTO_UPDATE_DEV : Auto-update data on device after computation
// _G_BOUNDS_CHK      : Do bounds checking
// _GLAPACK           : Set if using external Lapack API
// _G_VEC_CACHE_SIZE  : Sets vector (BLAS II) op cache blocking factor
// _G_MAT_CACHE_SIZE  : Sets vector op (BLAS III) cache blocking factor
// GEOFLOW_USE_OPENACC     : Set when using OpenACC
// GEOFLOW_USE_GPTL        : Set when using GPTL 

#if !defined(_G_TYPES_DEF)
#define _G_TYPES_DEF

//#define _G_IS2D
//#define _G_IS3D

// Set standard float type:
#define GFTYPE double

// Set standard 'compute' type (e.g., for basis function); may
// be the same as GFTYPE, but should not be of lower precision:
#define GCTYPE long double

// Define string type:
typedef std::string GString;


// Globally-defined defs:
// _G_AUTO_CREATE_DEV: Auto-create class data and class on device if ACC defined
// _G_BOUNDS_CHK     : Do data bounds check

// Basic data types:
#define GBYTE      unsigned char
#define GBOOL      bool
#define GUCHAR     unsigned char
#define GCHAR      char
#define GSHORT     short
#define GUSHORT    unsigned short
#define GINT       int
#define GUINT      unsigned int
#define GLONG      long
#define GULONG     unsigned  long
#define GLONGLONG  long long
#define GLLONG     GLONGLONG
#define GWORD      int
#define GDWORD     long
#define GKEY       size_t
#define GNODEID    long long
#define GSIZET     size_t
#define GFPOS      size_t
#define GFLOAT     float
#define GDOUBLE    double
#define GQUAD      long double


#if !defined(GBOOL)
#define GBOOL bool
#endif
#if !defined(TRUE)
#define TRUE  true
#endif
#if !defined(FALSE)
#define FALSE false
#endif

// Miscellaneous defs:
#if defined(__PGIC__) || defined(__GNUC__)
//#define NULLPTR    NULL
  #define NULLPTR    nullptr
#else
  #define NULLPTR    nullptr
#endif
#define GNULL_HANDLE    -1
#define PI               3.14159265358979323846264338328

// Miscellaneous macros:
#if !defined(_GDISTANCE_DEF)
#define _GDISTANCE_DEF
  #if defined(_G_IS1D)
    #define GDISTANCE(p1,p2) sqrt( ((p2)[0]-(p1)[0])*((p2)[0]-(p1)[0]) ) 
  #elif defined(_G_IS2D)
    #define GDISTANCE(p1,p2) sqrt( ((p2)[0]-(p1)[0])*((p2)[0]-(p1)[0])   \
                                  +((p2)[1]-(p1)[1])*((p2)[1]-(p1)[1])  )
  #elif defined(_G_IS3D)
    #define GDISTANCE(p1,p2) sqrt( ((p2)[0]-(p1)[0])*((p2)[0]-(p1)[0])   \
                                  +((p2)[1]-(p1)[1])*((p2)[1]-(p1)[1])   \
                                  +((p2)[2]-(p1)[2])*((p2)[2]-(p1)[2]) )
  #endif
#endif


#ifndef ABS
#define ABS(a)   ( a<0?-a:a )
#endif

#ifndef MIN
#define MIN(a,b)   ( a<b?a:b )
#endif

#ifndef MAX
#define MAX(a,b)   ( a>b?a:b )
#endif

#ifndef GSIGN
#define GSIGN(a,b) ( b>=0 ? ( a>=0 ? a:-a ) : ( a>=0 ? -a : a ) )
#endif


// Define datatypes (there is a corresponding 'communication' 
// set in gcommdata_t.h):
#if !defined(_G_DATATYPE_DEFTYPE_DEF)
#define  _G_DATATYPE_DEFTYPE_DEF
enum GD_DATATYPE        {GD_GBYTE=0,GD_GBOOL  ,GD_GCHAR   ,GD_GUCHAR ,
                         GD_GSHORT ,GD_GUSHORT,GD_GINT    ,GD_GUINT  ,
                         GD_GLONG  ,GD_GULONG ,GD_GLLONG  ,GD_GWORD  ,
                         GD_GDWORD ,GD_GKEY   ,GD_GNODEID ,GD_GSIZET ,
                         GD_GFPOS  ,GD_GFLOAT ,GD_GDOUBLE ,GD_GQUAD  };
#define GTYPE_NUM GD_GQUAD+1

const GINT GD_DATATYPE_SZ[] = 
                        {sizeof  (GBYTE),sizeof  (GBOOL),sizeof  (GCHAR),sizeof (GUCHAR),
                         sizeof (GSHORT),sizeof(GUSHORT),sizeof   (GINT),sizeof  (GUINT),
                         sizeof  (GLONG),sizeof (GULONG),sizeof (GLLONG),sizeof  (GWORD),
                         sizeof (GDWORD),sizeof   (GKEY),sizeof(GNODEID),sizeof (GSIZET),
                         sizeof  (GFPOS),sizeof (GFLOAT),sizeof(GDOUBLE),sizeof  (GQUAD)};
#endif

// Define reduction operations (there is a 'communication' set
// in gcommdata_t.h):
#if !defined(_GC_OP_DEF)
#define _GC_OP_DEF
enum GC_OP              {GC_OP_MAX=0,GC_OP_MIN ,GC_OP_SUM ,GC_OP_PROD,
                         GC_OP_LAND ,GC_OP_BAND,GC_OP_LOR ,GC_OP_BOR ,
                         GC_OP_LXOR ,GC_OP_BXOR};
#endif

#if !defined(_G_BDYTYPE_DEF)
#define _G_BDYTYPE_DEF
enum GBdyType                    { GBDY_DIRICHLET=0 , GBDY_INFLOW , GBDY_NOSLIP , GBDY_0FLUX , GBDY_OUTFLOW , GBDY_SPONGE , GBDY_PERIODIC , GBDY_NONE };
const char * const sGBdyType [] ={"GBDY_DIRICHLET"  ,"GBDY_INFLOW","GBDY_NOSLIP","GBDY_0FLUX","GBDY_OUTFLOW","GBDY_SPONGE","GBDY_PERIODIC","GBDY_NONE"};
#define GBDY_MAX GBDY_NONE
// DIRICHLET  : Acts on state directly; may be variable in space, but constant in time
// INFLOW     : Acts on state directly; may be variable in space, and time 
// NOSLIP     : Acts only on velocity, and is 0 on bdy, constant in time
// 0FLUX      : Acts only to velocity, normal compoment to surf is 0; time variable
// PERIODIC   : Acts on entire state; no user configuration required
// OUTFLOW    : Acts on state directly; variable in space, time; TBD
// SPONGE     : Acts on state directly; variable in space, time
#endif

#if !defined(_G_ELEMTYPE_DEF)
#define _G_ELEMTYPE_DEF
enum GElemType           {GE_REGULAR=0, GE_DEFORMED, GE_2DEMBEDDED}; // regular, deformed, embedded 2d
#define GE_MAX 3
#endif

#if !defined(_G_STEPPERTYPE_DEF)
#define _G_STEPPERTYPE_DEF
enum GStepperType        { GSTEPPER_EXRK=0 , GSTEPPER_BDFAB , GSTEPPER_BDFEXT };
const char * const sGStepperType[] =  
                         {"GSTEPPER_EXRK"  ,"GSTEPPER_BDFAB","GSTEPPER_BDFEXT"};
#define GSTEPPER_MAX 3
#endif

#if !defined(_G_VECTORTYPE_DEF)
#define _G_VECTORTYPE_DEF
enum GVectorType        {GVECTYPE_PHYS=0 , GVECTYPE_CONTRAVAR, GVECTYPE_COVAR};
#endif

#if !defined(_G_STATE_COMP_TYPE_DEF)
#define _G_STATE_COMP_TYPE_DEF
enum GStateCompType                    {GSC_KINETIC=0, GSC_MAGNETIC , GSC_DENSITYT , GSC_DENSITY , GSC_MASSFRAC , GSC_ENERGY , GSC_TEMPERATURE , GSC_PRESCRIBED , GSC_NONE };
const char * const sGStateCompType [] ={"GSC_KINETIC","GSC_MAGNETIC","GSC_DENSITYT","GSC_DENSITY","GSC_MASSFRAC","GSC_ENERGY","GSC_TEMPERATURE","GSC_PRESCRIBED","GSC_NONE"};
#define GSC_MAX GSC_NONE
#endif


// Variables used to set dimension & operational stack size:
#undef GDIM
#define G_MEMLOCNULL             -1     // Memblk index NULL value (<0)
#if defined(_G_IS1D)
  const GUSHORT GDIM = 1;               // GeoFLOW dimensionality
#elif defined(_G_IS2D)
  const GUSHORT GDIM = 2;               // GeoFLOW dimensionality
#elif defined(_G_IS3D)
  const GUSHORT GDIM = 3;               // GeoFLOW dimensionality
#else
  # error "_G_IS1D, _G_IS2D or _G_IS3D must be defined!"
#endif

// Math defs
#if !defined(sec)
  #define sec(a_rad)   (1.0/cos(a_rad))
#endif
#if !defined(csc)
  #define csc(a_rad)   (1.0/sin(a_rad))
#endif
#if !defined(cot)
  #define cot(a_rad)   (cos(a_rad)/sin(a_rad))
#endif

#if !defined FUZZYEQ
  #define FUZZYEQ(x,y,eps) ( (y <= (x+eps)) && (y >= (x-eps)) ) 
#endif

// Misc. defs
#define GMAX_ERROR_STRING 1024
#define BITSPERBYTE        8
#if !defined(GWORDSIZE_BITS)
  #define GWORDSIZE_BITS (sizeof(GWORD)*BITSPERBYTE)
#endif
#if !defined(GWORDSIZE_BYTES)
  #define GWORDSIZE_BYTES (GWORDSIZE_BITS / BITSPERBYTE)
#endif


// Macros for 'packing' unsigned ints:
#define _WLO 4
#if !defined(LOMASK)
  #define LOMASK      ( ~( (~(GUINT)0) << _WLO ) )
#endif
#if !defined(HIMASK)
  #define HIMASK      ( ~LOMASK ) 
#endif
#if !defined(GET_LOWORD)
  #define GET_LOWORD(a) ( a & LOMASK )
#endif
#if !defined(GET_HIWORD)
  #define GET_HIWORD(a) ( a >> _WLO )
#endif
#if !defined(SET_LOWORD)
  #define SET_LOWORD(a,b)   (a &= HIMASK);( a |= (b & LOMASK) )
#endif
#if !defined(SET_HIWORD)
  #define SET_HIWORD(a,b) (a &= LOMASK);( a |= ( b << _WLO ) )
#endif


#if !defined GError
  #define GError() printf("Error: %s; line: %d\n",__FILE__,__LINE__);
#endif

#define GPP(comm,a)                                   \
    do {                                              \
        std::stringstream oss_global_;                \
        oss_global_ << a << std::endl;                \
        GComm::cout(comm,oss_global_.str().c_str());  \
    } while(0)


#endif // _G_TYPES_DEF

#endif // !defined(_GTYPES_HPP)

