/*
 * mixer_factory.hpp
 *
 *  Created on: Mar 27, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename MixerFactory<ET>::MixBasePtr
MixerFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid){

	// Set the default mixer type
	const std::string default_mixer = "none";

	// Get the type of mixer
	const std::string mixer_name = ptree.getValue("mixer_name", default_mixer);

        // Set traits from prop tree:
        typename MixerBase<ET>::Traits traits;

        // Get whether 'mixing' correlation interval should be by 
        // cycle or time:
        std::string stype = ptree.getValue<std::string>("corr_itype","none");
        if      ( "cycle" == stype )  traits.itype = MixerBase<ET>::STIR_CYCLE;
        else if ( "time"  == stype )  traits.itype = MixerBase<ET>::STIR_TIME;
        else EH_ERROR("Invalid mixer correlation type specified");

        traits.state_index   = ptree.getArray<int>        ("state_index");    // state ids to 'mix' [0, 1, 2...]
        traits.corr_cycle    = ptree.getValue<size_t>     ("corr_cycle");     // correlation cycle
        traits.corr_time     = ptree.getValue<double>     ("corr_time");      // correlation time 
     
	// Create the mixer and cast to base type
	MixBasePtr base_ptr;
	if( "none" == mixer_name ){
		using MixImpl = NullMixer<Equation>;

		// Allocate mixer Implementation
		std::shared_ptr<MixImpl> mix_impl(new MixImpl(traits, grid));

		// Set back to base type
		base_ptr = mix_impl;
	}
	else {
		EH_ERROR("Requested mixer not found: " << mixer_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

