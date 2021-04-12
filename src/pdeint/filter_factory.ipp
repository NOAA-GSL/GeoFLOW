/*
 * filter_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename FilterFactory<ET>::FilterBasePtr
FilterFactory<ET>::build(const tbox::PropertyTree& eqn_ptree, const std::string filterblock, Grid& grid){

        // Note: the entry eqn_ptree is the equation property tree!


	// Get the prop tree specified by filterblock. 
        // If it doesn't exist in eqn_ptree, it's an error:
	std::string filter_type;
        PropertyTree ftree;
        if ( !eqn_ptree.isPropertyTree(filterblock) ) {
          cout << "FilterFactory::build: PropertyTree does not exist: " << filterblock << endl;
          exit(1);
        }
        else {
          ftree = eqn_ptree.getPropertyTree(filterblock);
	  // Get the type of filter
  	  filter_type = ftree.getValue<std::string>("filter_type", "none");
        }


	FilterBasePtr base_ptr;
	if( "none" == filter_type 
         || ""     == filter_type ) {
		using FilterImpl = NullFilter<Types>;

		// Allocate filter Implementation
                std::shared_ptr<FilterImpl> filter_impl(new FilterImpl());

		// Set back to base type
		base_ptr = filter_impl;

	}
        else if( "boyd_filter" == filter_type ) {
		using FilterImpl = GBoydFilter<Types>;
                typename GBoydFilter<ET>::Traits traits;

	        traits.ifilter  = ftree.getValue<int>("ifilter", 2);
	        traits.mufilter = ftree.getValue<double>("muifilter", 0.05);

		// Allocate filter Implementation
		std::shared_ptr<FilterImpl> filter_impl(new FilterImpl(traits, grid));
		// Set back to base type
		base_ptr = filter_impl;
        }
        else if( "projection_filter" == filter_type ) {
		using FilterImpl = GProjectionFilter<Types>;
                typename GProjectionFilter<ET>::Traits traits;

	        traits.istate = ftree.getArray   <int>("istate");
	        traits.pdelta = ftree.getArray   <int>("pdelta");
	        traits.alpha  = ftree.getArray<double>("alpha");

		// Allocate filter Implementation
		std::shared_ptr<FilterImpl> filter_impl(new FilterImpl(traits, grid));
		// Set back to base type
		base_ptr = filter_impl;

        }
	else {
		EH_ERROR("Requested filter type not found: " << filter_type);
	}

	return base_ptr;
}



} // namespace pdeint
} // namespace geoflow

