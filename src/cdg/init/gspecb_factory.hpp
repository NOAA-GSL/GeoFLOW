//==================================================================================
// Module       : gspecb_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECB_FACTORY_HPP)
#define _GSPECB_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gspecb.hpp"
#include "gspecb_user.hpp" // include user-provided spec functions


using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

class GSpecBFactory
{
  public:


        static GBOOL dospec(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, const GINT id,  GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy);

  private:
}; // class GSpecBFactory


#include "gspecb_factory.ipp"

#endif 
