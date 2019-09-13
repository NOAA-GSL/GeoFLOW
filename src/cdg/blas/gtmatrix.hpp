//==================================================================================
// Module       : gtmatrix_decl.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template matrix object, whose data is composed of
//                a regular array of contiguous data, ordered like a Fortran
//                matrix in row-major order (row data changing fastest over column data).
//                The implementation file, gtvector.ipp, is not included
//                in this file, just the declarations.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(_GTMATRIX_DECL_HPP)
#define _GTMATRIX_DECL_HPP
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <iostream>
#include "gtypes.h"
#include "gcommdata_t.h"
#include "gtvector.hpp"

template<class T> class GTMatrix;
template<class T> std::ostream &operator<<(std::ostream &, const GTMatrix<T> &);


#if !defined(GTMATRIX_ROTATE)
#define GTMATRIX_ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
        a(k,l)=h+s*(g-h*tau);
#endif


#if !defined(_G_MAT_CACHE_SIZE)
  # define _G_MAT_CACHE_SIZE 16
#endif



template<class T> class GTMatrix
{
public:
                          GTMatrix();
                          GTMatrix(const GSIZET , const GSIZET );
                          GTMatrix(const GSIZET);
                          GTMatrix(T *array, GSIZET  n1, GSIZET  n2);
                          GTMatrix(GTMatrix<T> &m, GSIZET *ind, GSIZET  nn);
                          GTMatrix(const GTMatrix<T> &);

                         ~GTMatrix();

         GTMatrix<T>       &operator=(const GTMatrix<T> &) ;
         void               operator=(T);
         void               operator=(const GTVector<T> &);


         void               operator+=(const T);
         void               operator-=(const T);
         void               operator*=(const T);            // multiplication by const Mat *= const

         void               operator+=(const GTMatrix &);
         void               operator-=(const GTMatrix &);

         GTMatrix        operator*(const T);             // right multiplication: Matrix * const

         GTVector<T>        operator*(const GTVector<T> &); // right multiplication: Matrix * Vector

         GTMatrix        operator+(const GTMatrix &); // addition
         GTMatrix        operator-(const GTMatrix &); // subtraction
         GTMatrix        operator*(const GTMatrix &); // right multiplication: Matrix * Matrix




         GBOOL              transpose(GTMatrix<T> &);
         GBOOL              transpose(T *&, GSIZET nx, GSIZET ny);
//       void               transpose();
         GBOOL              inverse(GTMatrix<T> &);
         GBOOL              inverse(GTMatrix<T> &, GSIZET n1, GSIZET n2);
         GTMatrix<T>       &inverse();
         GBOOL              isSymmetric();
         void               setSingularZero(GDOUBLE tol);
         void               setCacheBlockingFactor(GINT icsz);
         void               set(T);
         void               set(T *, GSIZET n1, GSIZET n2);
         GSIZET             distinctrow( GSIZET irow, T *&val, GSIZET *&icol, GSIZET &n);
         GSIZET             distinctrow( GSIZET irow, GTVector<T> &val, GSZBuffer &icol);
         GSIZET             distinctrow_floor( GSIZET irow, T *&val, GSIZET *&icol, GSIZET &n, T floor);
         GSIZET             distinctrow_floor( GSIZET irow, GTVector<T> &val, GSZBuffer &icol, T floor);
         GSIZET             distinctcol( GSIZET icol, T *&val, GSIZET *&irow, GSIZET &n);
         GSIZET             distinctcol( GSIZET icol, GTVector<T> &val, GSZBuffer &irow);
         GSIZET             distinctcol_floor( GSIZET icol, T *&val, GSIZET *&irow, GSIZET &n, T floor);
         GSIZET             distinctcol_floor( GSIZET icol, GTVector<T> &val, GSZBuffer &irow, T floor);
         GSIZET             multiplicity   (T val);
         GSIZET             multiplicity   (T val, GSIZET *&irow, GSIZET *&icol, GSIZET &n);
         GSIZET             multiplicity   (T val, GSZBuffer &irow, GSZBuffer &icol);
         GSIZET             distinctr(GSIZET row, GSIZET, GSIZET *&col, GSIZET &n);


         // Device/accelerator data methods
         void              updatehost(); 
         void              updatedev(); 

         inline T  &operator()(const GSIZET i, const GSIZET j) {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::&(): ";
           if ( i >= n1_ || j >= n2_ || j < 0  ) {
             std::cout << serr << "access error: i=" << i << "; j=" << j << std::endl;
             while(1);
             exit(1);
           }
           #endif
      //   return data_[j+i*n2_]; 
           return data_[i+j*n1_]; 
         }; 

         inline T &operator()(const GSIZET i) {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::&(): ";
           if ( i >= n1_*n2_ ) {
             std::cout << serr << "access error: i=" << i << std::endl;
             while(1);
             exit(1);
           }
           #endif
           return data_[i];
         }; 

         inline T operator()(const GSIZET i, const GSIZET j) const {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::() const: ";
           if ( i >= n1_ || j >= n2_ || j < 0  ) {
             std::cout << serr << "access error: i=" << i << "; j=" << j << std::endl;
             exit(1);
           }
           #endif
           T ret = data_[i+j*n1_];
           return ret;
         }; 

         inline T operator()(const GSIZET i) const {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::() const: ";
           if ( i >= n1_*n2_ ) {
             std::cout << serr << "access error: i=" << i << std::endl;
             exit(1);
           }
           #endif
           T ret = data_[i];
           return ret;
         }; 


         inline GTVector<T> &data() { return this->data_; }
         inline const GTVector<T> &data() const { return this->data_; }
         inline T *rowdata(GSIZET i) { return (this->data_.data()+i*n2_); }
         inline T *coldata(GSIZET j) { return (this->data_.data()+j*n1_); }

         GSIZET            size(GINT idir) const;  // _total_ size in idir 
         GSIZET            dim(GINT idir) const;   // matrix dims in idir
         GBOOL             resize(GSIZET  Nx, GSIZET  Ny);
         GBOOL             resizem(GSIZET  Nx, GSIZET  Ny);
         void              zero();
         void              createIdentity();
         GSIZET            multiplity(T val);
         void              multiplity(T val, GSIZET *&index);
  
         friend std::ostream &operator<<(std::ostream &os, GTMatrix<T> &obj ) {
           GLONG i,j;
           if ( obj.dim(1) == 0 || obj.dim(2) == 0 ) return os;
           os << std::endl << "[ ";
           for ( i=0; i<obj.dim(1); i++ ) {
             if ( i == 0 ) os << "["; else if ( i > 0 ) os << "  [" ;
             for ( j=0; j<obj.dim(2)-1; j++ ) {
               os
//                 << std::setiosflags(std::ios::scientific)
                   << std::setw(8)
                   << std::setprecision(5)
                   << std::setiosflags(std::ios::fixed)
                   << obj(i,j)
                   << ", ";
               }
               os
//                 << std::setiosflags(std::ios::scientific)
                   << std::setw(8)
                   << std::setprecision(5)
                   << std::setiosflags(std::ios::fixed)
                   << obj(i,j) ; 
             if ( i < obj.dim(1)-1 ) os << " ]; " << std::endl;
             else os  << " ]  ]" ;
           }
           return os;
           } // end of << operator 




GBOOL             ludcmp (T **&a, GSIZET  n, GLLONG  *&indx, T *d);
GBOOL             wludcmp(T **&a, GSIZET  n, GLLONG  *&indx, T *d);
GBOOL             lubksb (T **&A, GSIZET  n, GLLONG  *&indx, T b[]);
GBOOL             improve(T **&a, T **alud, GSIZET n, GLLONG *&indx, T b[], T x[]);
GBOOL             svdcmp(T **a, GSIZET m, GSIZET n, T w[], T **v);
GBOOL             svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1);
GBOOL             svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1, GSIZET n1, GSIZET n2);
void              svbksub(T **u, T w[], T **v, GSIZET n, GSIZET m, T b[], T x[]);
void              choldc(T **u, T w[], GSIZET n);
void              choldc(GTMatrix<T> &a, GTVector<T> &p);
void              jacobi(GTVector<T> &d, GTMatrix<T> &v, GSIZET &nrot, GTMatrix<T> &a, GTVector<T> &b, GTVector<T> &z);


// Private methods:
private:
T                 dpythag(T a, T b);
GSIZET            isamax(GSIZET n, T *sx, GSIZET  incx);
GTVector<T>      matvec_impl_(const GTVector<T> &obj, std::true_type);
GTVector<T>      matvec_impl_(const GTVector<T> &obj, std::false_type);
GTMatrix<T>      matmat_impl_(const GTMatrix<T> &obj, std::true_type);
GTMatrix<T>      matmat_impl_(const GTMatrix<T> &obj, std::false_type);

// Private data:
GSIZET            n1_;
GSIZET            n2_;
GINT              icsz_;  // GTMatrix cache-blocking factor

GTVector<T>       data_;
GDOUBLE           singzero_;

};


typedef GTMatrix<GFTYPE>  GMatrix;
typedef GTMatrix<GINT>    GIMatrix;
typedef GTMatrix<GNODEID> GNIDMatrix;
typedef GTMatrix<GLLONG>  GLLMatrix;
typedef GTMatrix<GSIZET>  GSZMatrix;
typedef GTMatrix<GFLOAT>  GFMatrix;
typedef GTMatrix<GDOUBLE> GDMatrix;
typedef GTMatrix<GQUAD>   GQMatrix;

#include "gtmatrix.ipp"

#endif
