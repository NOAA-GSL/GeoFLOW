//==================================================================================
// Module       : ggfx.hpp
// Date         : 1/6/21 (BTF)
// Description  : Encapsulates the methods and data associated with
//                a geometry--free global exchange (GeoFLOW Geometry-Free eXchange)
//                operator
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(GGFX_HPP)
#define GGFX_HPP

#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gcomm.hpp"


#include "boost/mpi.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include "tbox/spatial.hpp"
#include "tbox/pio.hpp"

#include <array>
#include <limits>
#include <map>
#include <set>
#include <vector>


// GGFX reduction operation defs:
#if !defined(GGFX_OP_DEF)
#define GGFX_OP_DEF
enum GGFX_OP {GGFX_OP_SUM=0, GGFX_OP_PROD, GGFX_OP_MAX, GGFX_OP_MIN, GGFX_OP_SMOOTH};
#endif


template<typename ValueType>
class GGFX {

public:
	GGFX()                             = default;
	GGFX(const GGFX& other)            = default;
	GGFX(GGFX&& other)                 = default;
	~GGFX()                            = default;
	GGFX& operator=(const GGFX& other) = default;
	GGFX& operator=(GGFX&& other)      = default;

	template<typename Coordinates>
	GBOOL init(const ValueType tolerance, Coordinates& xyz);

	template<typename ValueArray>
	GBOOL doOp(ValueArray& u,  GGFX_OP op);

	GC_COMM  getComm() const;

	template<typename CountArray>
	void get_imult(CountArray& imult) const;

private:
	using rank_type  = int;
	using size_type  = std::size_t;
	using value_type = ValueType;

	std::map<rank_type, std::set<size_type>>                      send_map_; // [Rank][1:Nsend] = Local Index
	std::map<rank_type, std::map<size_type, std::set<size_type>>> recv_map_; // [Rank][1:Nrecv][1:Nshare] = Local Index

	std::map<rank_type, std::vector<value_type>>             send_buffer_;      // [Rank][1:Nsend] = Value sending
	std::map<rank_type, std::vector<value_type>>             recv_buffer_;      // [Rank][1:Nrecv] = Value received
	std::vector<std::vector<size_type>>                      reduction_buffer_; // [1:Nlocal][1:Nreduce] = Value to reduce
};



namespace detail_extractor {
template<typename PairType>
struct pair_extractor {
	using bound_type = typename PairType::first_type;
	using key_type   = typename PairType::second_type;

	const bound_type& operator()(PairType const& pair) const {
		return pair.first;
	}
};
}


template<typename T>
template<typename Coordinates>
GBOOL
GGFX<T>::init(const T tolerance, Coordinates& xyz){
	using namespace geoflow::tbox;

	// Types
	using index_bound_type     = tbox::spatial::bound::Box<value_type,GDIM>;
	using index_key_type       = std::size_t; // local array index
	using index_value_type     = std::pair<index_bound_type,index_key_type>;   // (bound, key)
	using index_extractor_type = detail_extractor::pair_extractor<index_value_type>;
	using shared_index_type    = tbox::spatial::shared::index::RTree<index_value_type, index_extractor_type>;

	// Get MPI communicator, etc.
	namespace mpi = boost::mpi;
	mpi::communicator world;
	auto my_rank   = world.rank();
	auto num_ranks = world.size();

	// ----------------------------------------------------------
	//      Build Indexer for each Local Coordinate Location
	// ----------------------------------------------------------
	pio::perr << "Tolerance = " << tolerance << std::endl;
	pio::perr << "Build Indexer for each Local Coordinate Location" << std::endl;

	// Build "index_value_type" for each local coordinate and place into a spatial index
	shared_index_type local_bound_indexer;
	std::vector<index_bound_type> local_bounds_by_id;
	local_bounds_by_id.reserve(xyz.size());
	for(auto i = 0; i < xyz.size(); ++i){

		// Build a slightly larger bounding box
		auto min_xyz = xyz[i];
		auto max_xyz = xyz[i];
		for(auto d = 0; d < min_xyz.size(); ++d){
			min_xyz[d] -= tolerance;
			max_xyz[d] += tolerance;
		}

		// Place into coordinate index tree
		local_bounds_by_id.emplace_back(index_bound_type(min_xyz,max_xyz));
		auto index_value = index_value_type(local_bounds_by_id.back(), index_key_type(i));
		local_bound_indexer.insert(index_value);
	}

	// ----------------------------------------------------------
	//          Search Indexer for each Local Matches
	// ----------------------------------------------------------
	pio::perr << "Search Indexer for each Local Matches" << std::endl;

	//
	// For each of our local bounds we need to find the
	// - Matching local bounds
	//
	auto& my_rank_recv_map = recv_map_[my_rank];
	//my_rank_recv_map.resize(local_bounds_by_id.size());
	std::vector<bool> local_already_mapped(local_bounds_by_id.size(),false);
	for(size_type id = 0; id < local_bounds_by_id.size(); ++id){

		// Search for local mapping if not already found
		if(not local_already_mapped[id]){
			std::vector<index_value_type> search_results;
			local_bound_indexer.query( tbox::spatial::shared::predicate::Intersects(local_bounds_by_id[id]), std::back_inserter(search_results) );

			// Build list of just the Local ID's
			std::set<size_type> search_ids;
			for(auto& found_pair: search_results){
				search_ids.insert(found_pair.second);
			}

			// The result "should" be the same for every location within results
			for(auto& match_id: search_ids){
				my_rank_recv_map[match_id] = search_ids;
				local_already_mapped[match_id] = true;
			}
		}
	}

	// ----------------------------------------------------------
	//  Gather/Scatter my Local Coordinate Others who Overlap
	// ----------------------------------------------------------
	pio::perr << "Gather/Scatter my Local Coordinate Others who Overlap" << std::endl;

	// Get Bounding Box for this Ranks Coordinates
	// Perform MPI_Allgather so everyone gets a bound region for each processor
	std::vector<index_bound_type> bounds_by_rank;
	mpi::all_gather(world, local_bound_indexer.bounds(), bounds_by_rank);

	// For each Ranks Bounding Region
	// - Build list of local indexed values to send to each overlapping rank
	// - Submit a receive request to get indexed values from the rank
	// - Submit a send request to send indexed values from this rank
	std::map<rank_type, mpi::request> send_requests;
	std::map<rank_type, mpi::request> recv_requests;
	std::map<rank_type, std::vector<index_value_type>> send_to_ranks;
	std::map<rank_type, std::vector<index_value_type>> recv_from_ranks;
	for(auto rank = 0; rank < num_ranks; ++rank){
		if( rank != my_rank ){ // Don't search myself (we know the answer is everything)

			std::vector<index_value_type> search_results;
			local_bound_indexer.query( tbox::spatial::shared::predicate::Intersects(bounds_by_rank[rank]), std::back_inserter(search_results) );

			// If found any bounds within the ranks bounding region
			// - Move into buffer map (i.e. no copy)
			// - Submit a receive request
			// - Submit a send request
			if( search_results.size() > 0 ){
				send_to_ranks.emplace(rank, search_results);

				// Tag = 0
				recv_requests[rank] = world.irecv(rank, 0, recv_from_ranks[rank]);
				send_requests[rank] = world.isend(rank, 0, send_to_ranks[rank]);
			}
		}
	}

	// ----------------------------------------------------------
	//     Search Local Indexer for each Global Match
	// ----------------------------------------------------------
	pio::perr << "Search Local Indexer for each Global Match" << std::endl;

	//
	// For each of our global bound we received we need to find the
	// - Matching local bounds
	//
	// We already required a Local Coordinate Indexer for local matches
	// so we'll re-use that Indexer for searching against each global
	// coordinate location we received.
	//
	for(auto& [rank, pairs_from_rank]: recv_from_ranks){
		recv_requests[rank].wait(); // Wait to receive data

		// Loop over each remote pair<bound, id>
		for(auto& [remote_bnd, remote_id]: pairs_from_rank){
			std::vector<index_value_type> search_results;
			local_bound_indexer.query( tbox::spatial::shared::predicate::Intersects(remote_bnd), std::back_inserter(search_results));

			for(auto& [local_bnd, local_id]: search_results){
				send_map_[rank].insert(local_id);
				recv_map_[rank][remote_id].insert(local_id);
			}

		}
	}
	for(auto& [rank, req]: send_requests){
		req.wait(); // Clear out the send requests
	}



	auto local_bounded_region = local_bound_indexer.bounds();
	pio::perr << "Report:" << std::endl;
	pio::perr << "Tolerance = " << tolerance << std::endl;
	for(auto d = 0; d < GDIM; ++d){
		pio::perr << "Range " << d << "  " << local_bounded_region.min(d) << " - " << local_bounded_region.max(d) << std::endl;
	}
	for(auto& [rank, vec] : send_to_ranks){
		pio::perr << "Sending " << vec.size() << " bounds to rank " << rank << std::endl;
	}
	for(auto& [rank, vec] : recv_from_ranks){
		pio::perr << "Received " << vec.size() << " bounds from rank " << rank << std::endl;
	}
	world.barrier();
	pio::perr << std::endl;

	for(auto rank = 0 ; rank < num_ranks; ++rank){
		if(my_rank == rank){

			for(auto& [srank, vec] : send_map_){
				pio::perr << "  Sending to rank " << srank << " nodes " << vec.size() << std::endl;;
			}
			for(auto& [rrank, vec] : recv_map_){
				pio::perr << "  Recving to rank " << rrank << " nodes " << vec.size() << std::endl;;
			}

			pio::perr << std::flush;
		}
		world.barrier();
	}
	world.barrier();
	pio::perr << std::endl;


	std::map<size_type,size_type> count;
	std::vector<size_type> imult(xyz.size());
	this->get_imult(imult);
	for(auto& im : imult){
		if(count.count(im) == 0){
			count[im] = 1;
		}
		else {
			count[im]++;
		}
	}
	for(auto rank = 0 ; rank < num_ranks; ++rank){
		if(my_rank == rank){
			for(auto& [mult, sz] : count){
				pio::perr << "imult " << mult << "  size = " << sz << std::endl;
			}
		}
		world.barrier();
	}


	world.barrier();
	exit(0);

	return true;
}



template<typename T>
template<typename ValueArray>
GBOOL
GGFX<T>::doOp(ValueArray& u,  GGFX_OP op){



	return true;
}

template<typename T>
GC_COMM
GGFX<T>::getComm() const {
	return boost::mpi::communicator();
}

template<typename T>
template<typename CountArray>
void
GGFX<T>::get_imult(CountArray& imult) const {

	// Zero whole array
	for(size_type i = 0; i < imult.size(); ++i){
		imult[i] = 0;
	}

	// Accumulate the counts in each index
	// - Loop over buffers received from each rank
	// - Loop over each value within the buffer
	// - Loop over each local index that value contributes to
	for (auto& [rank, matrix_map] : recv_map_) {
		for(auto& [remote_id, local_id_set] : matrix_map) {
			for(auto& id : local_id_set) {
				++(imult[id]);
			}
		}
	}
}

#endif

