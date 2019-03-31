/*
 * stirrer_factory.hpp
 *
 *  Created on: MAr 27, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename StirrerFactory<ET>::StirBasePtr
StirrerFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid){


	// Set the default stirrer type
	const std::string default_stirrer = "none";

	// Get the type of stirrer
	const std::string stirrer_name = ptree.getValue("stirrer_name", default_stirrer);

        // Set traits from prop tree:
        typename StirrerBase<ET>::Traits traits;

        // Get whether 'stirring' correlation interval should be by 
        // cycle or time:
        std::string stype = ptree.getValue<std::string>("corr_itype","none");
        if      ( "cycle" == stype )  traits.itype = StirrerBase<ET>::STIR_CYCLE;
        else if ( "time"  == stype )  traits.itype = StirrerBase<ET>::STIR_TIME;
        else EH_ERROR("Invalid stirrer correlation type specified");

        traits.state_index   = ptree.getArray<int>        ("state_index");    // state ids to 'stir' [0, 1, 2...]
        traits.corr_cycle    = ptree.getValue<size_t>     ("corr_cycle");     // correlation cycle
        traits.corr_time     = ptree.getValue<double>     ("corr_time");      // correlation time 
     
	// Create the stirrer and cast to base type
	StirBasePtr base_ptr;
	if( "none" == stirrer_name ){
		using StirImpl = NullStirrer<Equation>;

		// Allocate stirrer Implementation
		std::shared_ptr<StirImpl> stir_impl(new StirImpl(traits, grid));

		// Set back to base type
		base_ptr = stir_impl;
	}
	else {
		EH_ERROR("Requested stirrer not found: " << stirrer_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

