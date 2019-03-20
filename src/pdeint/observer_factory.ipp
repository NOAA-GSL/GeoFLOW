/*
 * observer_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename ObserverFactory<ET>::ObsBasePtr
ObserverFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid){


	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = ptree.getValue("observer_name", default_observer);

        // Set traits from prop tree:
        typename ObserverBase<ET>::Traits traits;

        std::string stype = ptree.getValue<std::string>("itype");
        if      ( "cycle" == stype )  traits.itype = ObserverBase<ET>::OBS_CYCLE;
        else if ( "time " == stype )  traits.itype = ObserverBase<ET>::OBS_TIME;
        else EH_ERROR("Invalid observer type specified");

        traits.state_index   = ptree.getArray<int>        ("state_index");
        traits.state_names   = ptree.getArray<std::string>("state_names");
        traits.cycle_interval= ptree.getValue<size_t>     ("cycle_interval");
        traits.time_interval = ptree.getValue<double>     ("time_interval");
     
	// Create the observer and cast to base type
	ObsBasePtr base_ptr;
	if( "none" == observer_name ){
		using ObsImpl = NullObserver<Equation>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(traits, grid));

		// Set any parameters we need to set
		// NA

		// Set back to base type
		base_ptr = obs_impl;
	}
        else if( "gout_simple_observer" == observer_name ) {
		using ObsImpl = GOutSimpleObserver<Equation>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(traits, grid));

		// Set back to base type
		base_ptr = obs_impl;
        }
	else {
		EH_ERROR("Requested observer not found: " << observer_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

