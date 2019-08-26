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

        int itmp;

	// Set the default eqution type
	const std::string default_equation = "none";

	// Get the equation name
	std::string equation_name = ptree.getValue("pde_name", default_equation);

        // Set traits from prop tree:
        typename GBurgers<ET>::Traits btraits;

	// Set the default state components to force:
	std::vector<int> comps, default_comps;

        PropertyTree eqn_ptree = ptree.getPropertyTree(equation_name);
        PropertyTree stp_ptree = ptree.getPropertyTree("stepper_props");
        PropertyTree dis_ptree = ptree.getPropertyTree("dissipation_traits");

     
	// Create the equation and cast to base type
	EqnBasePtr base_ptr;
	if( "pde_burgers" == equation_name ){
		using EqnImpl = GBurgers<ET>;

                btraits.doheat    = eqn_ptree.getValue<bool>  ("doheat",false);
                btraits.bpureadv  = eqn_ptree.getValue<bool>  ("bpureadv",false);
                btraits.bconserved= eqn_ptree.getValue<bool>  ("bconserved",false);
                btraits.bforced   = eqn_ptree.getValue<bool>  ("use_forcing",false);
                btraits.variabledt= stp_ptree.getValue<bool>  ("variable_dt",false);
                btraits.courant   = stp_ptree.getValue<double>("courant",0.5);
                btraits.itorder   = stp_ptree.getValue<int>   ("time_deriv_order",4);
                btraits.inorder   = stp_ptree.getValue<int>   ("extrap_order",2);
                btraits.ssteptype = stp_ptree.getValue<std::string>
                                                             ("stepping_method","GSTEPPER_EXRK");
                btraits.nu        = dis_ptree.getValue<double>("nu");
                for ( auto i=0; i<GDIM; i++ ) default_comps.push_back(i);
                comps            = eqn_ptree.getArray<int>   ("forcing_comp",default_comps);
                btraits.iforced.resize(comps.size());
                btraits.iforced   = comps; // traits.iforced may be a different d.structure
		// Allocate equation Implementation
		std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(grid, btraits, utmp));

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

