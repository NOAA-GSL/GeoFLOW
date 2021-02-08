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

#include "tbox/assert.hpp"
#include "tbox/spatial.hpp"
#include "tbox/pio.hpp"
#include "tbox/tracer.hpp"

#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>



template<typename ValueType>
class GGFX {

public:

	struct Max{
		template<typename VectorType>
		typename VectorType::value_type
		operator()(const VectorType& vec) const {
			using value_type = typename VectorType::value_type;
			constexpr auto small = std::numeric_limits<value_type>::lowest();
			return std::accumulate(vec.begin(),vec.end(),small,[](auto& a, auto& b){
				return std::max(a,b);
			});
		}
	};

	struct Min{
		template<typename VectorType>
		typename VectorType::value_type
		operator()(const VectorType& vec) const {
			using value_type = typename VectorType::value_type;
			constexpr auto large = std::numeric_limits<value_type>::max();
			return std::accumulate(vec.begin(),vec.end(),large,[](auto& a, auto& b){
				return std::min(a,b);
			});
		}
	};

	struct Smooth{
		template<typename VectorType>
		typename VectorType::value_type
		operator()(const VectorType& vec) const {
			using value_type = typename VectorType::value_type;
			return std::accumulate(vec.begin(),vec.end(),value_type(0)) / vec.size();
		}
	};

	struct Sum{
		template<typename VectorType>
		typename VectorType::value_type
		operator()(const VectorType& vec) const {
			using value_type = typename VectorType::value_type;
			return std::accumulate(vec.begin(),vec.end(),value_type(0));
		}
	};




public:
	GGFX()                             = default;
	GGFX(const GGFX& other)            = default;
	GGFX(GGFX&& other)                 = default;
	~GGFX()                            = default;
	GGFX& operator=(const GGFX& other) = default;
	GGFX& operator=(GGFX&& other)      = default;

	template<typename Coordinates>
	GBOOL init(const ValueType tolerance, Coordinates& xyz);

	template<typename ValueArray, typename ReductionOp>
	GBOOL doOp(ValueArray& u,  ReductionOp op);

	GC_COMM  getComm() const;

	template<typename CountArray>
	void get_mult(CountArray& imult) const;

	template<typename CountArray>
	void get_imult(CountArray& imult) const;

	void display() const;

private:
	using rank_type  = int;
	using size_type  = std::size_t;
	using value_type = ValueType;
	static constexpr size_type MAX_DUPLICATES = std::pow(2,GDIM);

	std::map<rank_type, std::set<size_type>>                      send_map_; // [Rank][1:Nsend] = Local Index
	std::map<rank_type, std::map<size_type, std::set<size_type>>> recv_map_; // [Rank][1:Nrecv][1:Nshare] = Local Index

	std::map<rank_type, std::vector<value_type>>                  send_buffer_;      // [Rank][1:Nsend] = Value sending
	std::map<rank_type, std::vector<value_type>>                  recv_buffer_;      // [Rank][1:Nrecv] = Value received
	std::vector<std::vector<value_type>>                          reduction_buffer_; // [1:Nlocal][1:Nreduce] = Value to reduce
};



template<typename T>
void GGFX<T>::display() const {
	using namespace geoflow::tbox;
	namespace mpi = boost::mpi;
	mpi::communicator world;
	auto my_rank   = world.rank();
	auto num_ranks = world.size();
	for(auto rank = 0 ; rank < num_ranks; ++rank){
		if(my_rank == rank){
			pio::plog << "send_map_.size() = " << send_map_.size() <<std::endl;
			for(auto& [srank, vec] : send_map_){
				pio::plog << "  Sending to rank " << srank << " nodes " << vec.size() << std::endl;;
			}
			pio::plog << "recv_map_.size() = " << recv_map_.size() <<std::endl;
			for(auto& [rrank, vec] : recv_map_){
				pio::plog << "  Recving to rank " << rrank << " nodes " << vec.size() << std::endl;;
			}
			pio::plog << std::flush;
		}
		world.barrier();
	}
}

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
	GEOFLOW_TRACE();
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
//	pio::perr << "Tolerance = " << tolerance << std::endl;
//	pio::perr << "Build Indexer for each Local Coordinate Location" << std::endl;

	// Build "index_value_type" for each local coordinate and place into a spatial index
	GEOFLOW_TRACE_START("Build Local Index");
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
	GEOFLOW_TRACE_STOP();

	// ----------------------------------------------------------
	//          Search Indexer for each Local Matches
	// ----------------------------------------------------------
//	pio::perr << "Search Indexer for each Local Matches" << std::endl;

	//
	// For each of our local bounds we need to find the
	// - Matching local bounds
	//
	GEOFLOW_TRACE_START("Index Local Nodes");
	auto& my_rank_recv_map = recv_map_[my_rank];
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
	GEOFLOW_TRACE_STOP();

	// ----------------------------------------------------------
	//  Gather/Scatter my Local Coordinate Others who Overlap
	// ----------------------------------------------------------
//    pio::perr << "Gather/Scatter my Local Coordinate Others who Overlap" << std::endl;

	// Get Bounding Box for this Ranks Coordinates
	// Perform MPI_Allgather so everyone gets a bound region for each processor
	GEOFLOW_TRACE_START("AllGather Bounds");
	std::vector<index_bound_type> bounds_by_rank;
	const index_bound_type bnd = local_bound_indexer.bounds();
	mpi::all_gather(world, bnd, bounds_by_rank);
	GEOFLOW_TRACE_STOP();

	// For each Ranks Bounding Region
	// - Build list of local indexed values to send to each overlapping rank
	// - Submit a receive request to get indexed values from the rank
	// - Submit a send request to send indexed values from this rank
	GEOFLOW_TRACE_START("Gather/Scatter Coordinates");
	std::map<rank_type, mpi::request> send_requests;
	std::map<rank_type, mpi::request> recv_requests;
	std::map<rank_type, std::vector<index_value_type>> send_to_ranks;
	std::map<rank_type, std::vector<index_value_type>> recv_from_ranks;
	for(rank_type rank = 0; rank < num_ranks; ++rank){
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
	GEOFLOW_TRACE_STOP();

	// ----------------------------------------------------------
	//     Search Local Indexer for each Global Match
	// ----------------------------------------------------------
//	pio::perr << "Search Local Indexer for each Global Match" << std::endl;

	//
	// For each of our global bound we received we need to find the
	// - Matching local bounds
	//
	// We already required a Local Coordinate Indexer for local matches
	// so we'll re-use that Indexer for searching against each global
	// coordinate location we received.
	//
	GEOFLOW_TRACE_START("Build Global Maps");
	for(auto& [rank, req]: recv_requests){
		req.wait(); // Wait to receive data
		auto& pairs_from_rank = recv_from_ranks[rank];

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
	GEOFLOW_TRACE_STOP();

	GEOFLOW_TRACE_START("Allocate Buffers");
	// Allocate Buffers for Send/Receive in the future
	for(auto& [rank, send_ids]: send_map_){
		send_buffer_[rank].resize(send_ids.size());
	}
	for(auto& [rank, recv_id_pairs]: recv_map_){
		recv_buffer_[rank].resize(recv_id_pairs.size());
	}
	reduction_buffer_.resize(xyz.size());
	for(auto& buf : reduction_buffer_){
		buf.reserve(MAX_DUPLICATES);
	}
	GEOFLOW_TRACE_STOP();

	world.barrier(); // TODO: Remove

	return true;
}



template<typename T>
template<typename ValueArray, typename ReductionOp>
GBOOL
GGFX<T>::doOp(ValueArray& u,  ReductionOp oper){
	GEOFLOW_TRACE();
	using namespace geoflow::tbox;

	// Get size of data being used
	const auto N = u.size();

	// Get MPI communicator, etc.
	namespace mpi = boost::mpi;
	mpi::communicator world;
	auto my_rank   = world.rank();
	auto num_ranks = world.size();

	// Submit the Non-Blocking receive requests
	GEOFLOW_TRACE_START("Submit Receive Requests");
	std::vector<mpi::request> recv_requests;
	for(auto& [rank, mapping]: recv_map_){
		if( rank != my_rank ){
			auto tag = rank + my_rank;
			auto req = world.irecv(rank, tag, recv_buffer_[rank]);
			recv_requests.emplace_back(req);
		}
	}
	GEOFLOW_TRACE_STOP();

	// Copy values into Send Buffers & Non-Block Send
	std::vector<mpi::request> send_requests;
	GEOFLOW_TRACE_START("Pack Send Buffers");
	for (auto& [rank, local_ids_to_send] : send_map_) {

		// Pack into sending buffer
		size_type i = 0;
		auto& buffer_to_send = send_buffer_[rank];
		buffer_to_send.resize(local_ids_to_send.size());
		for (auto& id : local_ids_to_send){
			ASSERT(id < N);
			buffer_to_send[i++] = u[id];
		}

		// Send buffer to receiving processor
		auto tag = rank + my_rank;
		auto req = world.isend(rank, tag, buffer_to_send);
		send_requests.emplace_back(req);
	}
	GEOFLOW_TRACE_STOP();

//	// Prepare the buffer used for the reduction
//	// - This should only resize during first call
	reduction_buffer_.resize(N);
	for(auto& buf : reduction_buffer_){
		buf.resize(0);
		buf.reserve(MAX_DUPLICATES);
	}

	// Insert local data into the reduction buffer
	GEOFLOW_TRACE_START("Pack Local Data");
	for(auto& [local_id, local_id_set]: recv_map_[my_rank]){
		for(auto& id: local_id_set){
			reduction_buffer_[id].push_back( u[local_id] );
		}
	}
	GEOFLOW_TRACE_STOP();

	// Insert global data into the reduction buffer
	GEOFLOW_TRACE_START("Pack Global Data");
	size_type recv_count = 0;
	for (auto& [rank, remote_local_map] : recv_map_) {
		if( rank != my_rank ){
			recv_requests[recv_count++].wait();
			auto& buffer_for_rank = recv_buffer_[rank];

			// Loop over each value received from rank
			size_type n = 0;
			for(auto& [remote_id, local_id_set] : remote_local_map){

				// Loop over each local ID the value is used in
				for(auto& id : local_id_set) {
					ASSERT(reduction_buffer_.size() > id);
					ASSERT(reduction_buffer_[id].size() < MAX_DUPLICATES);
					reduction_buffer_[id].push_back( buffer_for_rank[n] );
				}
				++n;
			}
		}
	}
	GEOFLOW_TRACE_STOP();

	// Call the Reduction Operation on each
	GEOFLOW_TRACE_START("Reduce");
	size_type i = 0;
	for(auto& dup_values : reduction_buffer_) {
		u[i++] = oper(dup_values);
	}
	GEOFLOW_TRACE_STOP();

	// Clear all send requests
	mpi::wait_all(send_requests.begin(), send_requests.end());

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
GGFX<T>::get_mult(CountArray& mult) const {
	GEOFLOW_TRACE();

	// Zero whole array
	for(size_type i = 0; i < mult.size(); ++i){
		mult[i] = 0;
	}

	// Accumulate the counts in each index
	// - Loop over buffers received from each rank
	// - Loop over each value within the buffer
	// - Loop over each local index that value contributes to
	for (auto& [rank, matrix_map] : recv_map_) {
		for(auto& [remote_id, local_id_set] : matrix_map) {
			for(auto& id : local_id_set) {
				ASSERT(id < mult.size());
				++(mult[id]);
			}
		}
	}
}

template<typename T>
template<typename CountArray>
void
GGFX<T>::get_imult(CountArray& imult) const {
	GEOFLOW_TRACE();

	// Get multiplicity
	this->get_mult(imult);

	// Invert whole array
	for(size_type i = 0; i < imult.size(); ++i){
		imult[i] = static_cast<T>(1.0) / imult[i];
	}
}

#endif

