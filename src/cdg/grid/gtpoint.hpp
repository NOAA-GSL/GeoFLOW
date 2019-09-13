//==================================================================================
// Module       : gtpoint
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining template 'point' object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtypes.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <typeinfo>
#include <vector>
#include "gtvector.hpp"


template <typename T> class GTPoint;
template<typename T> std::ostream &operator<<(std::ostream &, GTPoint<T> &);

#if !defined(GTPOINT_HPP)
#define GTPOINT_HPP

template<typename T> class GTPoint
{
private:
  GINT           gdim_;
  T              eps_;
  GTVector<T*>   px_;


public:
                  GTPoint(); 
                  GTPoint(GINT dim); 
                  GTPoint(GINT idim, T eps); 
                  GTPoint(const GTPoint<T> &e); 
                 ~GTPoint();

  T x1;
  T x2;
  T x3;
  T x4;

         void      resize(GINT dim);
         void      setBracket(T eps);
         T         getBracket();
  inline GTVector<T*> &data() 
  { return px_;}

  inline GINT     size() 
  { return gdim_;}
  inline GINT     dim() 
  { return gdim_;}


  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GBOOL b;     
    for ( GINT i=0,b=TRUE;i<gdim_;i++) { 
      b = b && FUZZYEQ(*px_[i], pp[i], eps_); }  // do check for 'equality'
    return b; }

  template<typename U=T> typename std::enable_if<!std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GBOOL b; GTPoint<T> pq = pp;
    for ( GINT i=0,b=TRUE;i<gdim_;i++) b = b && ( pq[i] == *px_[i] );
    return b; }

  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, T>::type
  norm()  // Euclidean norm
  { T xn=0; for (GINT i=0; i<gdim_; i++ ) xn += ((*px_[i])*(*px_[i])); return sqrt(xn); }

  inline GBOOL    operator!=(const GTPoint<T> &p) 
  { return !this->operator==(p); }

  inline void  operator=(const GTPoint<T> &p)
  { eps_ = p.eps_; x1 = p.x1; x2 = p.x2; x3 = p.x3; x4 = p.x4;}

  inline void  operator=(const std::vector<T> &p)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] = p[j];
    if ( gdim_>0) x1 = p[0]; if ( gdim_>1) x2 = p[1]; if (gdim_>2) x3 = p[2]; if ( gdim_>3) x4 = p[3];}

  inline void  operator=(const GTVector<T> &p)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] = p[j];
    if ( gdim_>0) x1 = p[0]; if ( gdim_>1) x2 = p[1]; if (gdim_>2) x3 = p[2]; if ( gdim_>3) x4 = p[3];}

  inline GTPoint<T> &operator=(T f)
  { x1 = f;  x2 = f;  x3 = f;  x4 = f; return *this;}

  inline void assign(GTVector<GTVector<T>> &v, GSIZET i)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] = v[j][i];
    if ( gdim_>0) x1 = v[0][i]; if ( gdim_>1) x2 = v[1][i]; if (gdim_>2) x3 = v[2][i]; if ( gdim_>3) x4 = v[3][i];}

  inline GTPoint<T> operator-(const T f)
  {
      GTPoint ret(*this);
      ret -= f;
      return ret;
  }

  inline GTPoint<T> operator-(const GTPoint<T> &a)
  {
      GTPoint ret(*this);
      ret -= a;
      return ret;
  }

  inline void  operator-=(const T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] -= a; }

  inline void  operator-=(const GTPoint<T> &a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] -= a[j]; }

  inline GTPoint<T> operator+(const T f)
  {
      GTPoint ret(*this);
      ret += f;
      return ret;
  }

  inline GTPoint<T> operator+(const GTPoint<T> &a)
  {
      GTPoint ret(*this);
      ret += a;
      return ret;
  }

  inline void  operator+=(const T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] += a; }

  inline void  operator+=(const GTPoint<T> &a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] += a[j]; }

  inline GTPoint<T> operator*(const T a)
  {
      GTPoint ret(*this);
      ret *= a;
      return ret;
  }

  inline void  operator*=(const T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] *= a; }

  inline void  operator/=(const T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] /= a; }

  inline T &operator()(const GINT i)
  { 
#if defined(_G_BOUNDS_CHK)
    if ( i<0 || i>=gdim_ ) {
      while(1);  
      std::cout << "GTPoint::(): access error; bad index: " << i << std::endl; exit(1); }
#endif
    return *px_[i]; }

  inline T& operator[](const GINT i)
  { 
#if defined(_G_BOUNDS_CHK)
    if ( i<0 || i>=gdim_ ) { 
      while(1);
      std::cout << "GTPoint::[]: access error; bad index: " << i << std::endl; exit(1); }
#endif
    return *px_[i]; }

  inline const T& operator[](const GINT i) const
  { 
#if defined(_G_BOUNDS_CHK)
    if ( i<0 || i>=gdim_ ) { 
      while(1);
      std::cout << "GTPoint::[]: access error; bad index: " << i << std::endl; exit(1); }
#endif
    return *px_[i]; }

  inline T mag() { T v=0.0; for ( GINT j=0; j<gdim_; j++ ) v += pow(*px_[j],2.0); return sqrt(v); }

  friend std::ostream &operator<<(std::ostream &str, GTPoint<T> &obj) {
    str << "("  << *obj.px_[0];
    for ( GINT i=1; i<obj.gdim_; i++ ) {
    str << ", " << *obj.px_[i];} str << ")";
    return str;
  } // end, operator<<


};

typedef GTPoint<GFTYPE>  GFPoint;
typedef GTPoint<GINT>    GIPoint;

#include "gtpoint.ipp"

#endif
