//==================================================================================
// Module       : gtvector_decl.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types.
//                NOTE: Partial support for generalized 'index' that
//                      yields, start, stop, stride, pad indices.
//                The implementation file, gtvector.ipp, is not included
//                in this file, just the declarations.
//                NOTE: Math operations that accept a GTVector object
//                      argument may be treated as constant, by having
//                      size of input argument = 1, and taking as the
//                      operand only the first element.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(_GTVECTOR_DECL_HPP)
#define _GTVECTOR_DECL_HPP

#include <cstdlib>
#include <limits>
#include <iostream>
#include <vector>

#include "gtypes.h"
#include "gindex.hpp"
#include "cff_blas.h"
#include "gcomm.hpp"

#include "tbox/assert.hpp"

#if !defined(_G_VEC_CACHE_SIZE)
  # define _G_VEC_CACHE_SIZE 16
#endif


template <class T> class GTVector
{
  public:

    GTVector<T>();
    GTVector<T>(GSIZET n);
    GTVector<T>(GIndex &gin);
    GTVector<T>(GTVector<T> &obj);
    GTVector<T>(const GTVector<T> &obj);
    GTVector<T>(T *, GSIZET n, GSIZET istride=1);
    GTVector<T>(T *, GSIZET n, GSIZET istride, GBOOL bmanaged=TRUE);
   ~GTVector<T>();
    
    T *data();
    const T *data() const;
    GSIZET size() const;        // Return used buffer size
    GSIZET capacity() const;    // Return total available buffer size
    void   reserve(GSIZET n);   // Set capacity before having to initialize
    void   resize(GSIZET n);    // Resize data buffer
    void   resize(GIndex &);    // Resize data buffer with gindex
    void   resizem(GSIZET n);   // Resize data buffer only if n > current
    void   clear();             // Set capacity to 0
    void   push_back(T const &);// Push new element to end of data buffer
    T     &back();              // Get reference to last element
    T     &back() const;        // Get reference to last element

    void range(GLONG ibeg, GLONG end);      // Set range of vector within capacity
    void range_reset();                     // Reset range of vector 
    GIndex &getIndex() ;                    // Return generalized index member

    
    GBOOL              operator==(const GTVector<T> &b);
    
    GTVector<T>       &operator=(const GTVector<T> &b);
    
    GTVector<T>       &operator=(const std::vector<T> &b);
    
    void               operator=(T b);
    
    void               pointProd(const GTVector<T> &fact, GTVector<T> &ret);
    
    void               pointProd(const T a, const GTVector<T> &fact, GTVector<T> &ret);
    
    void               apointProd(const T a, const GTVector<T> &fact);
    
    void               pointProd(const GTVector<T> &);
    
    void               constProd(const T a, GTVector<T> &ret);
    
    void transpose(GSIZET n);

  
    
    void               set(T b);
    
    void               set(T *b, GSIZET n);

    void               floor(T b);
    // Device//accelerator data methods:
    void updatehost();
    void updatedev();

    inline T& operator[](const GSIZET i) {
      ASSERT_MSG(!( i+gindex_.beg() > gindex_.end() ), "i = " << i);
      return data_[i+gindex_.beg()];
    };

    inline const T& operator[](const GSIZET i) const {
      ASSERT_MSG(!( i+gindex_.beg() > gindex_.end() ), "i = " << i);
      return data_[i+gindex_.beg()];
    };

    
    GTVector       operator+(const T b);
    
    GTVector       operator-(const T b);
    
    GTVector       operator*(const T b);

     
    GTVector       operator+(const GTVector &b);
    
    GTVector       operator-(const GTVector &b);
     
    GTVector       operator*(const GTVector &b);
     
    T              dot(const GTVector &b);
     
    T              gdot(const GTVector &b, GC_COMM comm);
     
    T              gdot(const GTVector &b, const GTVector &c, GC_COMM comm);


    
    void               operator+=(const T b);
    
    void               operator-=(const T b);
    
    void               operator*=(const T b);
    
    void               operator/=(const T b);

     
    void               operator+=(const GTVector &b);
    
    void               operator-=(const GTVector &b);
    
    void               operator*=(const GTVector &b);
    
    void               operator/=(const GTVector &b);


    
        GBOOL isfinite();

    
        GBOOL isfinite(GSIZET &iwhere);
    
        T max();
    
        T amax();
    
        T amaxdiff(T tiny);
    
        T maxn(GSIZET n);
    
        GSIZET imax();
    
        T min();
    
        T amin();
    
        T amindiff(T tiny);
    
        T minn(GSIZET n);
    
        GSIZET imin();
    
        T sum();
    
        T sum(GSIZET ibeg, GSIZET iend);
    
        T infnorm();
    
        T Eucnorm();
    
        void rpow(GDOUBLE p);
    
        void abs();
    
inline GLONG findfirst(T val);    // Find first occurrence of val, else -1
    
inline GLONG findlast(T val);     // Find last occurrence of val, else -1
    
inline GBOOL onlycontains(T val); // Buffer contains only val?
    
inline GBOOL contains(T val);     // Buffer contains val?
    
inline GSIZET contains(T val, GSIZET *&iwhere, GSIZET &nw);
    
inline GBOOL containsn(T val, GSIZET n); // check n elements
    
inline GBOOL  contains(T val, GSIZET  &index);
    
inline GBOOL  containsn(T val, GSIZET n, GSIZET  &index);
    
inline GBOOL  contains_floor(T val, GSIZET &index, T floor, GSIZET istart=0);
    
inline GBOOL  contains_ceil(T val, GSIZET  &index, T ceil  , GSIZET istart=0);
    
inline GSIZET distinctrng(GSIZET istart, GSIZET n, GSIZET is,  T *&vals, GSIZET *&index, GSIZET  &n_distinct, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET distinctrng(GSIZET istart, GSIZET n, GSIZET is,  GSIZET *&index, GSIZET  &n_distinct, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET distinctrng_floor(GSIZET istart, GSIZET n, GSIZET is, T *&vals, GSIZET *&index, GSIZET  &n_distinct, T floor, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET distinctrng_floor(GSIZET istart, GSIZET n, GSIZET is, GSIZET *&index, GSIZET  &n_distinct, T floor, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET distinct(GSIZET *&index, GSIZET  &n_distinct, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET distinct_floor(GSIZET *&index, GSIZET  &n_distinct, T floor, T * const &tunique, GSIZET * const &itmp);
    
inline GSIZET multiplicity(T val);
    
inline GSIZET multiplicity_s(T val, GSIZET istart, GSIZET &ifound);
    
inline GSIZET multiplicity(T val, GSIZET *&index, GSIZET &n);
    
inline GSIZET multiplicity_s(T val, GSIZET istart, GSIZET *&index, GSIZET &n, GSIZET &ifound);
    
inline GSIZET multiplicity_floor(T val, T floor);
    
inline GSIZET multiplicity_ceil(T val, T ceil);
    
inline GSIZET multiplicity_floor(T val, GSIZET *&index, GSIZET &n, T floor);

    
       void   sortincreasing();
    
       void   sortincreasing(GTVector<GSIZET> &);
    
       void   sortdecreasing();
    
       void   sortdecreasing(GTVector<GSIZET> &);

    
       void concat(T *array, GSIZET n);

  private:

    GIndex gindex_; // gen. index object, changeable
    GIndex gindex_keep_; // gen. index object, stored
    T     *data_;
    GSIZET n_;
    GINT   icsz_;  // GTVector cache-blocking factor

    GBOOL  bdatalocal_; // tells us that data_ is owned by caller


   
inline GLLONG partitions2l(T *a, GLLONG start, GLLONG end);
   
inline GLLONG partitions2l(T *a, GSIZET *isort, GLLONG start, GLLONG end);
   
inline GLLONG partitionl2s(T *a, GLLONG start, GLLONG end);
   
inline GLLONG partitionl2s(T *a, GSIZET *isort, GLLONG start, GLLONG end);
   
inline void   quicksortl2s(T *a, GLLONG start, GLLONG end);
   
inline void   quicksortl2s(T *a, GSIZET *isort, GLLONG start, GLLONG end);
   
inline void   quicksorts2l(T *a, GLLONG start, GLLONG end);
   
inline void   quicksorts2l(T *a, GSIZET *isort, GLLONG start, GLLONG end);

   
   GTVector add_impl_(const GTVector &b, std::true_type);
   
   GTVector add_impl_(const GTVector &b, std::false_type);
  
   
   GTVector sub_impl_(const GTVector &b, std::true_type);
   
   GTVector sub_impl_(const GTVector &b, std::false_type);

  
   GTVector mul_impl_(const GTVector &b, std::true_type);
  
   GTVector mul_impl_(const GTVector &b, std::false_type);

  
   GTVector dot_impl_(const GTVector &b, std::true_type);
  
   GTVector dot_impl_(const GTVector &b, std::false_type);
};


//
// Writing to ostream doesn't need access to internal data so
// don't bother making it a friend class
template<typename T>
std::ostream &operator<<(std::ostream &os, GTVector<T> &obj) {
  if ( obj.data() == NULLPTR || obj.size() == 0 ) return os;
  os << obj[0];
  for ( GLONG j=1; j<obj.size(); j++ ) {
    os << " " << obj[j];
  }
  os << std::endl;
  return os;
};

typedef GTVector<GFTYPE>   GVector;
typedef GTVector<GINT>     GIVector;
typedef GTVector<GINT>     GIBuffer;
typedef GTVector<GNODEID>  GNIDBuffer;
typedef GTVector<GSIZET>   GSZBuffer;
typedef GTVector<GFLOAT>   GFVector;
typedef GTVector<GDOUBLE>  GDVector;
typedef GTVector<GDOUBLE>  GGVector;
typedef GTVector<GQUAD>    GQVector;

#include "gtvector.ipp"

#endif

