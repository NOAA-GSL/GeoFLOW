//==================================================================================
// Module       : gindex.hpp
// Date         : 5/10/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a generalized, 'global' index object
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================


#if !defined(GINDEX_HPP)
#define GINDEX_HPP
#include <iostream>
#include "gtypes.h"


class GIndex
{
// Private data:
private:

        GLONG              ibeg_;
        GLONG              iend_;
        GLONG              istride_;
        GLONG              ibase_;
        GLONG              npad_;
        GLONG              nlocal_;
        GLONG              nglob_;


public:
                           GIndex(GLONG ng, GLONG nl, GLONG ib, GLONG ie, GLONG inc, GLONG npad) 
                           {nglob_=ng; nlocal_=MAX(nl,ie-ib+1); ibeg_=ib; iend_=ie; istride_=inc; npad_=npad; ibase_=0;}  
                           GIndex() {ibeg_=0;iend_=0;nglob_=0;nlocal_=0;ibase_=0; istride_=1; npad_=0;}
                           GIndex(const GIndex &i) {nglob_=i.nglob_;  nlocal_=i.nlocal_; ibeg_=i.ibeg_; iend_=i.iend_; 
                                                     istride_=i.istride_; npad_=i.npad_; ibase_=i.ibase_;}
  virtual                 ~GIndex(){};

         inline GIndex   &operator()(GLONG ib, GLONG ie, GLONG inc, GLONG ibs) {
                           ibeg_ = ib; iend_ = ie; 
                           istride_=inc; ibase_=ibs; return *this;}
         inline GIndex   &operator()(GLONG ng, GLONG nl, GLONG ib, GLONG ie, GLONG inc, GLONG npad) {
                           nglob_ = ng>=0?ng:0; ibeg_ = ib; iend_ = ie; 
                           nlocal_=MAX(nl,iend_-ibeg_+1); istride_=inc; npad_=npad; return *this;}
         inline GIndex    &operator=(const GIndex &i) {
                           nglob_=i.nglob_; nlocal_=i.nlocal_; ibeg_=i.ibeg_; iend_=i.iend_; 
                           istride_=i.istride_; npad_=i.npad_;  ibase_=i.ibase_; return *this;}
         inline GBOOL      operator==(const GIndex &i) { 
                             return (ibeg_  ==i.ibeg_  &&iend_ ==i.iend_  &&istride_==i.istride_&&
                                     nlocal_==i.nlocal_&&nglob_==i.nglob_ &&ibase_  ==i.ibase_);  }
         inline GIndex    &operator+(GLONG i) {ibeg_ += i; iend_ += i; return *this;}


         friend std::ostream&   operator<<(std::ostream&str, const GIndex&i){
             str << "(" << i.szglobal() << "::" << i.szlocal() << ":" << i.beg() << ":" << i.end() << ":" << 
             i.stride() << ":" << i.pad() << ":" << i.base() << ")"; return str; }
         friend std::ostream&   operator<<(std::ostream&str, const GIndex*i){
             str << "(" << i->szglobal() << "::" << i->szlocal() << ":" << i->beg() << ":" << i->end() << ":" << 
             i->stride() << ":" << i->pad() << ":" << i->base() << ")"; return str; }

    
                   
         inline GLONG           &beg() {return ibeg_;}
         inline const GLONG     &beg() const {return ibeg_;}
         inline GLONG           &end() {return iend_;}
         inline const GLONG     &end() const {return iend_;}
         inline GLONG           &szglobal() {return nglob_;}
         inline const GLONG     &szglobal() const {return nglob_;}
         inline GLONG           &szlocal () {return nlocal_;}
         inline const GLONG     &szlocal () const {return nlocal_ ;}
         inline GLONG           &stride() {return istride_;}
         inline const GLONG     &stride() const {return istride_;}
         inline GLONG           &pad() {return npad_;}
         inline const GLONG     &pad() const {return npad_;}
         inline GLONG           &base() {return ibase_;}
         inline const GLONG     &base() const {return ibase_;}

        // Device//accelerator data methods:
        void updatehost(){};
        void updatedev (){};

};

#endif
