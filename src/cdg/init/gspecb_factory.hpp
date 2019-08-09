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

namespace geoflow {
namespace pdeint {


template<typename EquationType>
class GSpecBFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


        static GBOOL init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, const GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy);

  private:
}; // class GSpecBFactory




} // namespace pdeint
} // namespace geoflow


#include "gspecb_factory.ipp"

#endif 
