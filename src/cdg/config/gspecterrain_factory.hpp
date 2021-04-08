//==================================================================================
// Module       : gspecterrain_factory.hpp
// Date         : 3/11/20 (DLR)
// Description  : GeoFLOW terrain specification factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECTERRAIN_FACTORY_HPP)
#define _GSPECTERRAIN_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "pdeint/equation_base.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gterrain_specbox_user.hpp"
#include "gterrain_specsph_user.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GSpecTerrainFactory
{
  public:
        using Types         = TypePack;
        using State         = typename Types::State;
        using Grid          = typename Types::Grid;
        using GridBox       = typename Types::GridBox;
        using GridIcos      = typename Types::GridIcos;
        using Ftype         = typename Types::Ftype;


	static GBOOL spec(const PropertyTree& ptree, Grid &grid, State &utmp, State &xb, GBOOL &bterr);

  private:
	static GBOOL spec_box   (const PropertyTree& ptree, Grid &grid, State &utmp, GString stype, State &xb);
	static GBOOL spec_sphere(const PropertyTree& ptree, Grid &grid, State &utmp, GString stype,  State &xb);

}; // end, class GSpecTerrainFactory

#include "gspecterrain_factory.ipp"

#endif 
