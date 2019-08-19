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

using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

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


        static GBOOL dospec(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, const GINT id,  GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy);

  private:
}; // class GSpecBFactory


#include "gspecb_factory.ipp"

#endif 
