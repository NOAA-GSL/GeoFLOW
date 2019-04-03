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

        // Get whether 'observation' cadence should be by cycle or time:
        std::string stype = ptree.getValue<std::string>("cadence_itype","none");
        if      ( "cycle" == stype )  traits.itype = ObserverBase<ET>::OBS_CYCLE;
        else if ( "time"  == stype )  traits.itype = ObserverBase<ET>::OBS_TIME;
        else EH_ERROR("Invalid observer type specified");

        std::vector<std::string> default_names;
        std::vector<int> default_ids;
        traits.state_index   = ptree.getArray<int>        ("state_index",default_ids);    // state ids to 'observe' [0, 1, 2...]
        traits.state_names   = ptree.getArray<std::string>("state_names",default_names);  // state names 
        traits.cycle_interval= ptree.getValue<size_t>     ("cycle_interval", 10);         // cadence for cycle type
        traits.time_interval = ptree.getValue<double>     ("time_interval", 1.0);         // cadence for time type
        traits.dir           = ptree.getValue<std::string>("directory",".");  // directory
     
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
        else if( "posixio_observer" == observer_name ) {
		using ObsImpl = GPosixIOObserver<Equation>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(traits, grid));

		// Set back to base type
		base_ptr = obs_impl;
        }
        else if( "global_diag_basic" == observer_name ) {
		using ObsImpl = GGlobalDiag_basic<Equation>;

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

