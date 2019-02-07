//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(GGRID_FACTORY_HPP)
#define GGRID_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"


class GGridFactory
{
 public:

	static GGrid *build(const tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM comm);


} // namespace geoflow

#include "integrator_factory.ipp"

#endif /* SRC_PDEINT_INTEGRATOR_FACTORY_HPP_ */
