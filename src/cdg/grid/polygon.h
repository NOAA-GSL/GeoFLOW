//==================================================================================
// Module       : polygons.h 
// Date         : 8/31/18 (DLR)
// Description  : Header defining polygons used in GeoFLOW.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtpoint.hpp"

#include "tbox/tracer.hpp"

#if !defined(_GTRIANGLE_DEF)
# define _GTRIANGLE_DEF
template<typename T> class GTriangle {
public:
GTVector<GTPoint<T>*>  v;
GTPoint<T>   v1;
GTPoint<T>   v2;
GTPoint<T>   v3;
 GTriangle<T>() { 
   GEOFLOW_TRACE(); v.resize(3); v1.resize(GDIM); v2.resize(GDIM); v3.resize(GDIM);
   v[0] = &v1; v[1] = &v2; v[2] = &v3; }
 GTriangle<T>(GINT dim) { 
   GEOFLOW_TRACE(); v.resize(3); v1.resize(dim); v2.resize(dim); v3.resize(dim);
   v[0] = &v1; v[1] = &v2; v[2] = &v3; }
 ~GTriangle() {GEOFLOW_TRACE();}
 void operator=(const GTriangle &t) {GEOFLOW_TRACE(); v=t.v; v1=t.v1; v2=t.v2; v3=t.v3; }
 inline GTPoint<T> &operator[](const GINT i)
  { return *v[i]; }
 void resize(GINT sz) {GEOFLOW_TRACE(); v1.resize(sz); v2.resize(sz); v3.resize(sz); }

 friend std::ostream &operator<<(std::ostream &os, GTriangle<T> &obj) {
   os << "{" << obj.v1;
   for ( GLONG j=1; j<3; j++ ) {
     os << "; " << *(obj.v[j]);
   }
   os << "}" << std::endl;
   return os; } // end, operator<<
};
#endif

#if !defined(_GQUAD_DEF)
# define _GQUAD_DEF
template<typename T> class GQuad{
public:
GTVector<GTPoint<T>*>  v;
GTPoint<T>   v1;
GTPoint<T>   v2;
GTPoint<T>   v3;
GTPoint<T>   v4;
 GQuad<T>() {GEOFLOW_TRACE(); v.resize(4); v1.resize(3); v2.resize(3); v3.resize(3);
   v4.resize(3); v[0]=&v1; v[1]=&v2; v[2]=&v3; v[3]=&v4; }
 GQuad<T>(GINT dim) {GEOFLOW_TRACE(); v.resize(4); v1.resize(dim); v2.resize(dim); v3.resize(dim);
   v4.resize(dim); v[0]=&v1; v[1]=&v2; v[2]=&v3; v[3]=&v4; }
 ~GQuad() {GEOFLOW_TRACE();}
 void operator=(const GQuad &t) {GEOFLOW_TRACE(); v=t.v; v1=t.v1; v2=t.v2; v3=t.v3; v4=t.v4; }
 inline GTPoint<T> &operator[](const GINT i)
  { return *v[i]; }

 void resize(GINT sz) {GEOFLOW_TRACE(); v1.resize(sz); v2.resize(sz); v3.resize(sz);
      v4.resize(sz); }

 friend std::ostream &operator<<(std::ostream &os, GQuad<T> &obj) {
   os << "{" << obj.v1;
   for ( GLONG j=1; j<4; j++ ) {
     os << "; " << *(obj.v[j]);
   }
   os << "}" << std::endl;
   return os; } // end, operator<<

};
#endif

 
#if !defined(_GHEX_DEF)
# define _GHEX_DEF
template<typename T> class GHex{
public:
GTVector<GTPoint<T>*>  v;
GTPoint<T>   v1;
GTPoint<T>   v2;
GTPoint<T>   v3;
GTPoint<T>   v4;
GTPoint<T>   v5;
GTPoint<T>   v6;
GTPoint<T>   v7;
GTPoint<T>   v8;
 GHex<T>() {GEOFLOW_TRACE(); v.resize(8); v1.resize(3); v2.resize(3); v3.resize(3);
 v4.resize(3); v5.resize(3); v6.resize(3); v7.resize(3); v8.resize(3);
   v[0]=&v1; v[1]=&v2; v[2]=&v3; v[3]=&v4; v[4]=&v5; v[5]=&v6; v[6]=&v7; v[7]=&v8; }
 GHex<T>(GINT dim) {GEOFLOW_TRACE(); v.resize(8); v1.resize(dim); v2.resize(dim); v3.resize(dim);
 v4.resize(dim); v5.resize(dim); v6.resize(dim); v7.resize(dim); v8.resize(dim);
   v[0]=&v1; v[1]=&v2; v[2]=&v3; v[3]=&v4; v[4]=&v5; v[5]=&v6; v[6]=&v7; v[7]=&v8;}
 ~GHex() {GEOFLOW_TRACE();}
 void operator=(const GHex &t) {GEOFLOW_TRACE(); v=t.v; v1=t.v1; v2=t.v2; v3=t.v3; v4=t.v4;
   v5=t.v5; v6=t.v6; v7=t.v7; v8=t.v8; }
 void resize(GINT sz) {GEOFLOW_TRACE(); v1.resize(sz); v2.resize(sz); v3.resize(sz);
      v4.resize(sz); v5.resize(sz); v6.resize(sz); v7.resize(sz); v8.resize(sz); }

 friend std::ostream &operator<<(std::ostream &os, GHex<T> &obj) {
   os << "{" << obj.v1;
   for ( GLONG j=1; j<8; j++ ) {
     os << "; " << *(obj.v[j]);
   }
   os << "}" << std::endl;
   return os; } // end, operator<<

};
#endif

 
