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
        Observer_base::Traits traits;
        traits.itype = static_cast<ObsType>(ptree.getValue("itype"));
        traits.state_index = ptree.getArray<int>("state_index");
        traits.state_names= ptree.getArray<std::string>("state_names");
        traits.cycle_interval= ptree.getArray<size_t>("cycle_interval");
        traits.time_interval= ptree.getArray<double>("cycle_interval");
     
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
		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new GOutSimpleObserver(traits, grid));

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

