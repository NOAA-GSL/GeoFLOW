//==================================================================================
// Module       : gtpoint
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining template 'point' object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(GTPOINT_HPP)
#define GTPOINT_HPP

#include "gtypes.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <typeinfo>
#include <vector>
#include "gtvector.hpp"


template <typename U> class GTPoint;
template<typename U> std::ostream &operator<<(std::ostream &, GTPoint<U> &);


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
  inline GINT     dim() const 
  { return gdim_;}


  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GBOOL b=TRUE, f;     
    for ( auto i=0;i<gdim_;i++) { 
      f = FUZZYEQ(*px_[i], pp[i], eps_);
      b = b && f; 
    }  // do check for 'equality'
    return b; }

  template<typename U=T> typename std::enable_if<!std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GBOOL b=TRUE; GTPoint<T> pq = pp;
    for ( auto i=0;i<gdim_;i++) b = b && ( pq[i] == *px_[i] );
    return b; }

  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, T>::type
  norm()  // Euclidean norm
  { T xn=0; for ( auto i=0; i<gdim_; i++ ) xn += ((*px_[i])*(*px_[i])); return sqrt(xn); }

  inline void unit()  // Make into unit vector
  { T xn=this->norm(); for ( auto i=0; i<gdim_; i++ ) (*px_[i]) /= xn; }

  inline GBOOL    operator!=(const GTPoint<T> &p) 
  { return !this->operator==(p); }

  inline void  operator=(const GTPoint<T> &p)
  { assert(p.gdim_ == gdim_ );
    eps_ = p.eps_; for ( auto j=0; j<gdim_; j++ ) *px_[j] = p[j];
    if ( gdim_>0) x1 = p[0]; if ( gdim_>1) x2 = p[1]; if (gdim_>2) x3 = p[2]; if ( gdim_>3) x4 = p[3];}

  inline void  operator=(const std::vector<T> &p)
  { assert(p.size() >= gdim_ );
    for ( auto j=0; j<gdim_; j++ ) *px_[j] = p[j];
    if ( gdim_>0) x1 = p[0]; if ( gdim_>1) x2 = p[1]; if (gdim_>2) x3 = p[2]; if ( gdim_>3) x4 = p[3];}

  inline void  operator=(const GTVector<T> &p)
  { assert(p.size() == gdim_ );
    for ( auto j=0; j<gdim_; j++ ) *px_[j] = p[j];
    if ( gdim_>0) x1 = p[0]; if ( gdim_>1) x2 = p[1]; if (gdim_>2) x3 = p[2]; if ( gdim_>3) x4 = p[3];}

  inline GTPoint<T> &operator=(T f)
  { x1 = f;  x2 = f;  x3 = f;  x4 = f; return *this;}

  inline void assign(const GTVector<GTVector<T>> &v, GSIZET i)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] = 0.0;
    for ( auto j=0; j<v.size(); j++ ) *px_[j] = v[j][i];
  }

  inline void assign(const GTVector<GTVector<T>*> &v, GINT nv, GSIZET i)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] = 0.0;
    for ( auto j=0; j<nv; j++ ) *px_[j] = (*v[j])][i];
  }

  inline GTPoint<T> operator-(T f) {
      GTPoint ret(*this);
      ret -= f;
      return ret;
  }

  inline GTPoint<T> operator-(GTPoint<T> &a){
      GTPoint ret(*this);
      ret -= a;
      return ret;
  }

  inline void  operator-=(T a)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] -= a; }

  inline void  operator-=(const GTPoint<T> &a)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] -= a[j]; }

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
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] += a; }

  inline void  operator+=(const GTPoint<T> &a)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] += a[j]; }

  inline GTPoint<T> operator*(const T a)
  {
      GTPoint ret(*this);
      ret *= a;
      return ret;
  }

  inline T dot(const GTPoint<T> &a)
  {
      T ret = (*px_[0]) * a[0];
      for ( auto j=1; j<gdim_; j++ ) ret += (*px_[j]) * a[j];
      return ret;
  }

  inline T cross2(const GTPoint<T> &p)
  {
      // Compute ret_z = this X p for 2d vectors
      // (there's only a z-component)
       return (x1*p.x2 - x2*p.x1);
  }

  inline void cross(const GTPoint<T> &p, GTPoint<T> &ret)
  {
      assert(gdim_ == 3 && p.dim() == 3 && ret.dim() == 3 );

      // Compute ret = this X p
      ret.x1 = x2 * p.x3 - x3 * p.x2;
      ret.x2 = x3 * p.x1 - x1 * p.x3;
      ret.x3 = x1 * p.x2 - x2 * p.x1;
  }

  inline void  operator*=(const T a)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] *= a; }

  inline void  operator/=(const T a)
  { for ( auto j=0; j<gdim_; j++ ) *px_[j] /= a; }

  inline T &operator()(const GINT i)
  { 
    ASSERT_MSG(!( i<0 || i>=gdim_ ), "i = " << i);
    return *px_[i]; }

  inline T& operator[](const GINT i)
  { 
    ASSERT_MSG(!( i<0 || i>=gdim_ ), "i = " << i);
    return *px_[i]; }

  inline T operator[](const GINT i) const
  { 
    ASSERT_MSG(!( i<0 || i>=gdim_ ), "i = " << i);
    return *px_[i]; }

  inline T mag() { T v=0.0; for ( auto j=0; j<gdim_; j++ ) v += pow(*px_[j],2.0); return sqrt(v); }

  friend std::ostream &operator<<(std::ostream &str, GTPoint<T> &obj) {
    str << "("  << *obj.px_[0];
    for ( auto i=1; i<obj.gdim_; i++ ) {
    str << ", " << *obj.px_[i];} str << ")";
    return str;
  } // end, operator<<


};

typedef GTPoint<GFTYPE>  GFPoint;
typedef GTPoint<GINT>    GIPoint;

#include "gtpoint.ipp"

#endif
