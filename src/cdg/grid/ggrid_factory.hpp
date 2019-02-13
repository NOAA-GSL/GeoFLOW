//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_FACTORY_HPP)
#define _GGRID_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"


class GGridFactory
{
  public:

	static GGrid *build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM comm);

  private:
        static GBdyType      str2bdytype(const GString &stype);

}; // namespace geoflow

//#include "ggrid_factory.ipp"

#endif 
