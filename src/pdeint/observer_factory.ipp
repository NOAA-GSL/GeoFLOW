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
ObserverFactory<ET>::build(const tbox::PropertyTree& ptree, const std::string obsname, EqnBasePtr& equation, Grid& grid, IOBasePtr &pIO){

        // Note: the entry ptree is the main property tree!


	// Get the prop tree specified by obsname:
        PropertyTree obstree = ptree.getPropertyTree(obsname);

	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = obstree.getValue("observer_name", default_observer);

        // Set traits from prop tree:
        ObsTraitsType obstraits;
    
        // Set default derived quantity info:
        get_traits(ptree, obsname, obstraits);


	ObsBasePtr base_ptr;
	if( "none" == observer_name 
         || ""     == observer_name ){
		using ObsImpl = NullObserver<ET>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid));

		// Set back to base type
		base_ptr = obs_impl;
	}
    else if( "gio_observer" == observer_name ) {
		using ObsImpl = GIOObserver<ET>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid, pIO, obstraits));

                // Set some pIO traits from observer:
                pIO->get_traits().idir = obstraits.idir;
                pIO->get_traits().odir = obstraits.odir;
//              obs_impl->setIO(pIO);

		// Set back to base type
		base_ptr = obs_impl;
        }
    else if( "burgers_diag" == observer_name ) {
		using ObsImpl = GBurgersDiag<ET>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid, obstraits));

		// Set back to base type
		base_ptr = obs_impl;
        }
    else if( "gmconv_diag" == observer_name ) {
		using ObsImpl = GMConvDiag<ET>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid, obstraits));

		// Set back to base type
		base_ptr = obs_impl;
        }
	else {
		EH_ERROR("Requested observer not found: " << observer_name);
	}

	return base_ptr;
}


template<typename ET>
void ObserverFactory<ET>::get_traits(const tbox::PropertyTree& ptree, const std::string obsname, ObsTraitsType& traits){

        // Note: the entry ptree is the main property tree!


	// Get the prop tree specified by obsname:
        PropertyTree obstree = ptree.getPropertyTree(obsname);

	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = obstree.getValue("observer_name", default_observer);

        // Set default derived quantity info:
        PropertyTree dqtree;
        std::vector<int>         defi;
        std::vector<std::string> defq;
        std::vector<std::string> dqnames;
        

        // Get whether 'observation' cadence should be by cycle or time:
        std::string stype = obstree.getValue<std::string>("cadence_type","none");
        if      ( "cycle" == stype )  traits.itype = ObserverBase<ET>::OBS_CYCLE;
        else if ( "time"  == stype )  traits.itype = ObserverBase<ET>::OBS_TIME;
        else if ( "none" != stype ) {
          cout << "ObserverFactory<ET>::build: stype=" << stype << endl;
          EH_ERROR("Invalid observer type specified");
        }

        std::vector<std::string> defst_names = {"u1","u2","u3","u4","u5","u6","u7","u8","u9"};
        std::vector<std::string> defgr_names = {"xgrid","ygrid","zgrid"};;
        std::vector<int> def_ids;
        traits.treat_as_1d   = obstree.getValue<bool>       ("treat_as_1d",false);        // treat-as-1d flag
        traits.state_index   = obstree.getArray<int>        ("state_index",defi);         // state ids to 'observe' [0, 1, 2...]
        traits.state_names   = obstree.getArray<std::string>("state_names",defst_names);  // state names 
        traits.agg_state_name= obstree.getValue<std::string>("agg_state_name","state");  // state names 
        traits.grid_names    = obstree.getArray<std::string>("grid_names",defgr_names);   // grid comp names 
        traits.agg_grid_name = obstree.getValue<std::string>("agg_grid_name","grid");   // grid comp names 
        traits.cycle_interval= obstree.getValue<size_t>     ("cycle_interval", 10);       // cadence for cycle type
        traits.time_interval = obstree.getValue<double>     ("time_interval", 1.0);       // cadence for time type
        traits.freq_fact     = obstree.getValue<double>     ("interval_freq_fact", 1.0);  // freq factor relative to, say restart
        traits.idir          = obstree.getValue<std::string>("idir",".");          // input directory
        traits.odir          = obstree.getValue<std::string>("odir",".");         // outputdirectory
        traits.start_ocycle  = obstree.getValue<size_t>     ("start_ocycle",0);           // starting output cycle 
//      traits.start_cycle   = obstree.getValue<size_t>     ("start_cycle",0);            // start evol cycle
        traits.start_time    = obstree.getValue<double>     ("start_time",0.0);           // start evol time
     
	// Set traits that depend on observer type:
        if( "gio_observer" == observer_name ) {
		using ObsImpl = GIOObserver<ET>;

                // Fill derived quantities strutures, if any:
                dqnames       = obstree.getArray<std::string> ("derived_quantities",defq);  // list of derived quantities to output
                traits.derived_quantities.resize(dqnames.size());
                for ( auto j=0; j<dqnames.size(); j++ ) {
                  dqtree = ptree.getPropertyTree(dqnames[j]); // get prop tree for named quantity
                  traits.derived_quantities[j].icomponents = dqtree.getArray<int>("state_index");
                  traits.derived_quantities[j].snames      = dqtree.getArray<std::string>("names"  ,defq);
                  traits.derived_quantities[j].smath_op    = dqtree.getValue<std::string>("mathop" ,"");
                }

                // Fill state-derived quantities strutures, if any:
                dqnames       = obstree.getArray<std::string> ("state_derived_quantities",defq);  // list of derived quantities to output
                traits.state_derived_quantities.resize(dqnames.size());
                for ( auto j=0; j<dqnames.size(); j++ ) {
                  dqtree = ptree.getPropertyTree(dqnames[j]); // get prop tree for named quantity
//                traits.state_derived_quantities[j].icomponents = dqtree.getArray<int>("state_index");
                  traits.state_derived_quantities[j].snames      = dqtree.getArray<std::string>("names"  ,defq);
                  traits.state_derived_quantities[j].smath_op    = dqtree.getValue<std::string>("mathop" ,"");
                }
        }

}


} // namespace pdeint
} // namespace geoflow

