/*
 * equation_factory.ipp
 *
 *  Created on: July 8, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

#include "gburgers.hpp"

namespace geoflow {
namespace pdeint {


template<typename ET>
typename EquationFactory<ET>::EqnBasePtr
EquationFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid, State& utmp){


	// Set the default eqution type
	const std::string default_equation = "none";

	// Get the equation name
	std::string equation_name = ptree.getValue("pde_name", default_equation);

        // Set traits from prop tree:
        typename EquationBase<ET>::Traits traits;

	// Set the default state components to force:
	std::vector<int> comps, default_comps;

        PropertyTree eqn_ptree = ptree.getPropertyTree(equation_name);
        PropertyTree stp_ptree = ptree.getPropertyTree("stepper_props");

     
	// Create the equation and cast to base type
	EqnBasePtr base_ptr;
	if( "pde_burgers" == equation_name ){
		using EqnImpl = NullEquation<Equation>;

                traits.doheat    = eqn_ptree.getValue<bool>    ("doheat",false);
                traits.bpureadv  = eqn_ptree.getValue<bool>    ("bpureadv",false);
                traits.bconserved= eqn_ptree.getValue<bool>    ("bconserved",false);
                traits.bforced   = eqn_ptree.getValue<bool>    ("use_forcing",false);
                traits.variabledt= stp_ptree.getValue<bool>    ("variable_dt",false);
                traits.courant   = stp_ptree.getValue<double>  ("courant",0.5);
                traits.itorder   = stp_ptree.getValue<int>     ("time_deriv_order",4);
                traits.inorder   = stp_ptree.getValue<int>     ("extrap_order",2);
                traits.steptype  = stp_ptree.getValue<int>     ("stepping_method","GSTEPPER_EXRK");
                for ( auto i=0; i<GDIM; i++ ) default_comps.push_back(i);
                comps            = eqn_ptree.getArray<int>     ("forcing_comp",default_comps);
                traits.iforced.resize(comps.size);
                traits.iforced   = comps; // traits.iforced may be a different d.structure
		// Allocate equation Implementation
		std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(traits, grid, utmp));

		// Set back to base type
		base_ptr = eqn_impl;
	}
	else {
		EH_ERROR("Requested equation not found: " << equation_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

