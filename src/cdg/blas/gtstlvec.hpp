//==================================================================================
// Module       : gtstlvec.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types. Derives from std::vector.
//                NOTE: Consider adding generallized 'index' to mimic
//                      STL iterator (which can't be used bec of ACC).
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : std::vector
//==================================================================================

#include "gtypes.h"
#include <cstdlib>
#include <iostream>
#include <vector>


template <class T> class GTSTLVec;
template<class T> std::ostream &operator<<(std::ostream &, GTSTLVec<T> &);


#if !defined(_GTSTLVEC_HPP)
#define _GTSTLVEC_HPP


template <class T> class GTSTLVec: public std::vector<T>
{
  public:
    GTSTLVec<T>();
    GTSTLVec<T>(GSIZET n);
    GTSTLVec<T>(GTSTLVec<T> &obj);
    GTSTLVec<T>(const GTSTLVec<T> &obj);
    GTSTLVec<T>(T *, GSIZET n);
   ~GTSTLVec<T>();
    
    
    void resize(GSIZET n);
    void SetCacheBlockingFactor(GINT icsz);

    // Device//accelerator data methods:
    void updatehost();
    void updatedev();

/*
inline   T &operator[](const GSIZET i) {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTSTLVec<T>::operator[]: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return data_[i];
    };

inline    T operator[](const GSIZET i) const {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTSTLVec<T>::operator[] const: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      T ret = data_[i];
      return ret;
    };
*/

inline    T &operator()(const GSIZET i) {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTSTLVec<T>::operator(): ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return (*this)[i];
    };

inline    T operator()(const GSIZET i) const {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTSTLVec<T>::operator() const: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      T ret = (*this)[i];
      return ret;
    };

    #pragma acc routine vector
    GTSTLVec<T>       &operator=(const GTSTLVec<T> &b);
    #pragma acc routine vector
    GTSTLVec<T>       &operator=(T b);
    #pragma acc routine vector
    void               pointProd(GTSTLVec<T> &ret, const GTSTLVec<T> &);
    #pragma acc routine vector
    void               constProd(GTSTLVec<T> &ret, const T);
    #pragma acc routine vector
    void               add      (GTSTLVec<T> &ret, const GTSTLVec<T> &,  const T a, const T b);
    #pragma acc routine vector
    void               sub      (GTSTLVec<T> &ret, const GTSTLVec<T> &,  const T a, const T b);


    #pragma acc routine vector
    void               operator+=(T b);
    #pragma acc routine vector 
    void               operator+=(GTSTLVec<T> &b);
    #pragma acc routine vector
    void               operator-=(GTSTLVec<T> &b);
    #pragma acc routine vector
    void               operator-=(T b);
    #pragma acc routine vector
    void               operator*=(T b);
    #pragma acc routine vector
    void               operator*=(GTSTLVec<T> &b);
    const GTSTLVec<T>  range(GSIZET ibeg, GSIZET nn); // get range from beg index

    // Define ostream op here:
    friend std::ostream &operator<<(std::ostream &os, GTSTLVec<T> &obj) {
      T *d=obj.data_;
      for ( GSIZET j=0; j<obj.n_; j++ ) {
        os << d[j] << " ";
      }
      return os;
  }; // end, operator<<


  private:
    T     *data_;
    GSIZET n_;
    GINT   icsz_; // GTSTLVec cache-blocking factor 

};


#endif
