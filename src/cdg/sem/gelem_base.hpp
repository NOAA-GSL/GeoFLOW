//==================================================================================
// Module       : gelem_base
// Date         : 6/1/18 (DLR)
// Description  : Base class forming interfaces for all allowed 1D/2D/3D (lin/quad/hex) 
//                elements
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//               
//
//          V7 ________________ V6<->P1_
//            /|     M6       /|
//       M7  / |             / |
//          /  |    F5   M5 /  |
//         /M11|  M4   F2  /   |                     z 
//     V4 /____|__________/    | M10                       y
//        |    |         | V5  |                     |     
//        |    |     M2  | F1  |                     |   /
//        | V3 |_________|_____| V2                  |  /
//    M8  | F3 / F0      |    /                      | /
//        |   /       M9 |   /                       |/________  x
//        |  / M3   F4   |  / M1                     
//        | /            | /                       
//        |/_____________|/                      
// P0_<-> V0     M0     V1
//
// Faces are labeled s.t. F0-F3 correspond to orientation of edges on bottom plane;
// F4 and F5 are bottom and top faces, respectively. Nodes are labeled implicitly starting
// from V0, increasing fastest with x, then y, then z
//
// [Note that in 2d, we use just the bottom plane.]
//==================================================================================
#if !defined(GELEM_BASE_HPP)
#define GELEM_BASE_HPP 

#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gshapefcn_base.hpp"
#include "gtpoint.hpp"


typedef GTVector<GTVector<GTPoint<GFTYPE>>>  GVVFPoint;
typedef GTVector<GTVector<GINT>>             GVVInt;
typedef GTVector<GTVector<GUINT>>            GVVUInt;
typedef GTVector<GTVector<GBdyType>>         GVVBTyp;
typedef GTVector<GTVector<GFTYPE>>           GVVFType;
typedef GTVector<GFTYPE>                     GVFType;
typedef GTMatrix<GTVector<GFTYPE>>           GMVFType;


class GElem_base 
{
public:
                            enum ElemNodeType {VERTEX=0, EDGE, FACE, VOL};
                            GElem_base(GINT dim);
                            GElem_base(GINT dim, GElemType etype, GNBasis<GCTYPE,GFTYPE> *b1, GNBasis<GCTYPE,GFTYPE> *b2, GNBasis<GCTYPE,GFTYPE> *b3=NULLPTR);
                            GElem_base(GINT dim, GElemType etypa, GNBasis<GCTYPE,GFTYPE> *b[], GINT nb);
                            GElem_base(GINT dim, GElemType etype, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b);
                            GElem_base(const GElem_base &) = default;
                           ~GElem_base();

        void                operator=(const GElem_base&);

                   

virtual void                set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b); // must be called first
virtual void                set_elemtype(GElemType);          // set element type
virtual void                init(GVVFType &xNodes);           // must be called second
        void                set_nodes(GVVFType &xNodes);      // if modifying coords

virtual void                interp(GTVector<GFTYPE> &ufrom, GTVector<GTVector<GFTYPE>*> &xito, GTVector<GFTYPE> &uto);
virtual void                interp(GTVector<GFTYPE> &ufrom, GTVector<GTVector<GFTYPE>*> &xito, GTVector<GFTYPE> &uto, GTVector<GTMatrix<GFTYPE>> &matv, GTVector<GTMatrix<GFTYPE>> &matu, GTVector<GFTYPE> &tmp);
#if 0
virtual GBOOL               point_in_elem(GFTYPE x, GFTYPE y, GFTYPE z=0.0);
virtual GBOOL               point_in_elem(GFPoint p[], GINT n);
virtual GBOOL               any_point_in_elem(GTVector<GFPoint> &p);
virtual GFTYPE              integrate(GTVector<GFTYPE>  &v, GFTYPE *multiplicity=NULLPTR);
virtual GBOOL               differentiate(GTVector<GFTYPE>  &dv, GTVector<GFTYPE>  &v, GINT  idir);
virtual GBOOL               differentiateD(GTVector<GFTYPE> &dv, GTVector<GFTYPE> &v, GINT  idir, GNBasis<GCTYPE,GFTYPE>*> *dbasis);
virtual GBOOL               differentiateWithMass(GTVector<GFTYPE> &dv, GTVector<GFTYPE> &v, GTVector<GFTYPE> &tmp, GINT  idir);
virtual GBOOL               differentiateWeak(GTVector<GFTYPE> &dv, GTVector<GFTYPE> &v, GTVector<GFTYPE> &tmp, GINT  idir);
#endif

        void                set_elemid(GINT id){elemid_ = id;}
        void                set_rootid(GINT id){rootid_ = id;}
        void                resize(GTVector<GINT> &);

inline  GLONG              &igbeg(){ return igbeg_;} // vol beg index 
inline  GLONG              &igend(){ return igend_;} // vol end index
inline  GLONG              &ibbeg(){ return ibbeg_;} // bdy beg index
inline  GLONG              &ibend(){ return ibend_;} // bdy end index
inline  GLONG              &ifbeg(){ return ifbeg_;} // face beg index
inline  GLONG              &ifend(){ return ifend_;} // face end index
inline  GINT                nnodes(){ return Ntot_;} // no. vol nodes
inline  GINT                nfnodes(){ return Nftot_;} // no. face nodes 
inline  GTVector<GINT>     &dim(){ return N_;} 
inline  GINT                dim(GINT i){ return N_[i];} 
inline  GTVector<GINT>     &size(){ return N_;} 
inline  GINT                size(GINT i){ return N_[i];} 
inline  GINT                order(GINT i){ return N_[i]-1;} 
inline  GElemType           elemtype(){ return elemtype_;} 
inline  GINT                nvertices(){ return nVertices_;} 
inline  GINT                nedges(){ return nEdges_;} 
inline  GINT                nfaces(){ return nFaces_;}
inline  GINT                elemid(){ return elemid_;}
inline  GINT                rootid(){ return rootid_;}

inline  GVVInt             &vedge(){ return ivedge_;}
inline  GTVector<GINT>     &vedge(GINT i){ return ivedge_[i];}
inline  GVVInt             &fedge(){ return ifedge_;}
inline  GTVector<GINT>     &fedge(GINT i){ return ifedge_[i];}
inline  GVVInt             &eface(){ return ieface_;}
inline  GTVector<GINT>     &eface(GINT i){ return ieface_[i];}
inline  GVVInt             &vface(){ return ivface_;}
inline  GTVector<GINT>     &vface(GINT i){ return ivface_[i];}
inline  GTVector<GFPoint>  &edgeCentroid(){ return edgeCentroid_;}
inline  GFPoint            &edgeCentroid(GINT i){ return edgeCentroid_[i];}
inline  GTVector<GFPoint>  &faceCentroid(){ return faceCentroid_;}
inline  GFPoint            &faceCentroid(GINT i){ return faceCentroid_[i];}
inline  GFPoint            &elemCentroid(){ return elemCentroid_;}

inline  GTVector<GNBasis<GCTYPE,GFTYPE>*> 
                           &gbasis(){ return gbasis_;}
inline  GNBasis<GCTYPE,GFTYPE>
                           *gbasis(GINT i){ return  gbasis_[i];}
inline  GVVFType           &xNodes(){ return xNodes_;}
inline  GTVector<GFTYPE>   &xNodes(GINT i){ return xNodes_[i];}
inline  GTVector<GTVector<GFTYPE>*>
                           &xiNodes(){ return xiNodes_;}
inline  GTVector<GFTYPE>
                           &xiNodes(GINT i){ return *xiNodes_[i];}
inline  GTVector<GFPoint>  &xVertices(){ return xVertices_;}
inline  GFPoint            &xVertices(GINT i){ return xVertices_[i];}

inline  GFTYPE              volume(){ return volume_;} 
#if 0
inline  GTMatrix<GTVector<GFTYPE>>
                           &dXidX(){ return dXidX_; };
inline  GTVector<GFTYPE>   *dXidX(GINT i,GINT j){ return &dXidX_(i,j);}
inline  GTMatrix<GTVector<GFTYPE>*>
                           &gij(){ return gij_; };
inline  GTVector<GFTYPE>   *gij(GINT i,GINT j){ return gij_(i,j);}
inline  GTVector<GFTYPE>   &Jac(){ return Jac_; }
inline  GVVFType           &faceJac(){ return faceJac_; }
inline  GTVector<GFTYPE>   &faceJac(GINT i){ return faceJac_[i];}
#endif

        GBOOL               isvert(GSIZET ind, GUINT &ivert);
inline  GVVInt             &vert_indices(){ return vert_indices_;}
inline  GTVector<GINT>     &vert_indices(GINT i){ return vert_indices_[i];}
inline  GVVInt             &edge_indices(){ return edge_indices_;}
inline  GTVector<GINT>     &edge_indices(GINT i){ return edge_indices_[i];}
inline  GVVInt             &face_indices(){ return face_indices_;}
inline  GTVector<GINT>     &face_indices(GINT i){ return face_indices_[i];}
inline  GTVector<GFTYPE>   &face_mass(){ return face_mass_;}
//inline  GVVFType           &face_mass(){ return face_mass_;}
//inline  GTVector<GFTYPE>   &face_mass(GINT i){ return face_mass_[i];}
inline  GVVUInt            &face_desc(){ return face_desc_;}
inline  GTVector<GUINT>    &face_desc(GINT i){ return face_desc_[i];}
inline  GTVector<GINT>     &bdy_indices(){ return bdy_indices_;}
inline  GTVector<GBdyType> &bdy_types(){ return bdy_types_;}

virtual void                dogeom1d (GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &facejac); 
virtual void                dogeom2d (GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &facejac); 
virtual void                dogeom3d (GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &facejac); 

friend std::ostream&        operator<<(std::ostream&, GElem_base &);    // Output stream operator
 

protected:
virtual void                set_size();   // set data depended on basis and GDIM
virtual void                set_size1d(); 
virtual void                set_size2d(); 
virtual void                set_size3d(); 
virtual void                build_elem(); // create all geometry quantities dependent on xNodes_
virtual void                build_elem1d();
virtual void                build_elem2d();
virtual void                build_elem3d();

virtual void                get_indirect  (GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                                          GVVInt &edge_ind, GVVInt &face_ind, GVFType &face_mass, GVVUInt &face_desc);
virtual void                get_indirect1d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                                          GVVInt &edge_ind, GVVInt &face_ind, GVFType &face_mass, GVVUInt &face_desc);
virtual void                get_indirect2d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                                          GVVInt &edge_ind, GVVInt &face_ind, GVFType &face_mass, GVVUInt &face_desc);
virtual void                get_indirect3d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                                          GVVInt &edge_ind, GVVInt &face_ind, GVFType &face_mass, GVVUInt &face_desc);
        void                Jac(GMVFType &rij, GTVector<GFTYPE> &jac, GBOOL &pChk, GINT *ind, GINT nind);
//virtual void                Jac_embed(GMVFType &G, GTVector<GFTYPE> &jac, GBOOL &pChk, GINT *pind, GINT nind);
        void                inv(GMVFType &G, GMVFType &iG);
/*
        void                set_faceNormal2d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GINT> &iind, GTVector<GTVector<GFTYPE>> &bdyNormal);
        void                set_faceNormal3d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GINT> &iind, GTVector<GTVector<GFTYPE>> &bdyNormal);
        void                set_bdyNormal2d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GINT> &iind, GTVector<GTVector<GFTYPE>> &bdyNormal);
        void                set_bdyNormal3d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GINT> &iind, GTVector<GTVector<GFTYPE>> &bdyNormal);
*/
                                        

GINT                    dim_;           // elem dimension 1, 2, 3)
GINT                    Ntot_;          // total no. node points in element
GINT                    Nftot_;         // total no. face node points in element
GTVector<GINT>          N_;             // no. _nodes_ in each direction
GLONG                   igbeg_;         // holds starting global index 
GLONG                   igend_;         // holds ending global index
GLONG                   ibbeg_;         // holds starting global bdy node index 
GLONG                   ibend_;         // holds ending global bdy node index
GLONG                   ifbeg_;         // holds starting global face node index 
GLONG                   ifend_;         // holds ending global face node index
GBOOL                   bInitialized_;  // element initialized?
GBOOL                   bbasis_;        // basis set?
GElemType               elemtype_;      // elem type
GINT                    elemid_;        // elem id 
GINT                    rootid_;        // elem's root id
GINT                    nVertices_;     // # vertices 
GINT                    nEdges_;        // # edges 
GINT                    nFaces_;        // # faces

// Element id data:
GVVInt                  ivedge_;        // vertex indices comprising edges: vert id X edge
GVVInt                  ifedge_;        // face indices defining edges: face id X edge id
GVVInt                  ivface_;        // vert indices comprising faces: vert id X face id
GVVInt                  ieface_;        // edge indices comprising faces: index X face id


// Element grid data:
GTVector<GNBasis<GCTYPE,GFTYPE>*>      
                        gbasis_       ; // expansion basis
GVVFType                xNodes_       ; // phys loc of nodes
GTVector<GTVector<GFTYPE>*>
                        xiNodes_      ; // reference nodes [-1,1], in tensort prod form
GTVector<GFPoint>       xVertices_    ; // phys loc of vertices, in local coords
GTVector<GFTYPE>        Ni_           ; // tensor product shape function

// Geometry data:
GFTYPE                  volume_;        // volume (3d), area(2d)
GShapeFcn_base<GFTYPE>  *gshapefcn_;     // shape fcn object; depends on elemtype
GTVector<GFPoint>       edgeCentroid_;  // phys loc of edge centroids, global coords
GTVector<GFPoint>       faceCentroid_;  // phys loc of face centroids, global coords
GFPoint                 elemCentroid_;  // phys loc of elem centroid, global coords
#if 0
GTMatrix<GTVector<GFTYPE>*>
                        gij_;           // Sum_k dxi_i/dx_k dxi_j/dx_k;                
GTMatrix<GTVector<GFTYPE>>
                        dXidX_;         // dxi_i/dx_k 
GTVector<GFTYPE>        Jac_;           // volume Jacobian |Gij|
GVVFType                faceJac_;       // face Jaobians |Gij| on edges
GTVector<GVVFType>      bdyNormal_;    // normal to face at each node point (2d & 3d)
#endif

// Indirection indices:
GVVInt                  vert_indices_;  // all indices comprising vertices
GVVInt                  edge_indices_;  // all indices comprising edges          
GVVInt                  face_indices_;  // all indices comprising faces
GVVUInt                 face_desc_;     // 'description' of face nodes
GVFType                 face_mass_;     // weights at each face_indices node, indexed by vol index
GIBuffer                bdy_indices_;   // global bdy indices
GTVector<ElemNodeType>  bdy_nodetype_;  // node type of bdy nodes
GTVector<GBdyType>      bdy_types_;     // bdy types on any global bdy nodes
};

#endif
