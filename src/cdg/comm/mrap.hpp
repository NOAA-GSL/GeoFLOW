//==================================================================================
// Module       : mrap.hpp
// Date         : 1/5/21 (BTF)
// Description  : Matrix Reductions Across Processors for exchanging values at
//                shared nodes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(MRAP_HPP)
#define MRAP_HPP

#include <map>
#include <vector>


#include <boost/mpi.hpp>

namespace mpi = boost::mpi;


template<typename ValueType>
class MRAP {

public:
	using rank_type  = int;
	using size_type  = std::size_t;
	using value_type = ValueType;


	MRAP()                              = default;
	MRAP(const MRAP& other)             = default;
	MRAP(MRAP&& other)                  = default;
	~MRAP()                             = default;
	MRAP& operator=(const MRAP& other)  = default;
	MRAP& operator=(MRAP&& other)       = default;


	void add_send_mapping(const rank_type rank, const size_type N, const size_type* local_ids) {
		send_map_.emplace(rank,std::vector<size_type>(local_ids, local_ids+N));
		send_buffer_[rank].resize(N);
	}

	void add_send_mapping(const rank_type rank, const std::vector<size_type> local_ids) {
		send_map_.emplace(rank,local_ids);
		send_buffer_[rank].resize(local_ids.size());
	}


	void add_recv_mapping(const rank_type rank, const std::vector<std::vector<size_type>> local_ids) {
		recv_map_.emplace(rank,local_ids);
		recv_buffer_[rank].resize(local_ids.size());
	}

	template<typename ReductionOper>
	void reduce(const size_type N, const value_type* local_values, ReductionOper& op) {
		constexpr size_type MAX_DUPLICATES = 8;

		// Get MPI Communicator
		mpi::communicator world;

		// Submit the Non-Blocking receive requests
		std::vector<mpi::request> recv_requests;
		for (auto& [rank, vector_of_values] : recv_buffer_) {
			auto tag = rank + world.rank();
			auto req = world.irecv(rank, tag, vector_of_values);
			recv_requests.emplace_back(req);
		}

		// Copy values into Send Buffers & Send
		std::vector<mpi::request> send_requests;
		for (auto& [rank, vector_of_ids] : send_map_) {

			// Pack into sending buffer
			size_type i = 0;
			auto& buffer_for_rank = send_buffer_[rank];
			for (auto& id : vector_of_ids){
				ASSERT(id < N);
				buffer_for_rank[i++] = local_values[id];
			}

			// Send buffer to receiving processor
			auto tag = rank + world.rank();
			auto req = world.isend(rank, tag, buffer_for_rank);
			send_requests.emplace_back(req);
		}

		// Prepare the buffer used for the reduction
		// - This should only resize after first call
		size_type n = 0;
		reduction_buffer_.resize(N);
		for(auto& buf : reduction_buffer_){
			buf.resize(0);
			buf.reserve(MAX_DUPLICATES);
			buf.push_back(local_values[n++]);
		}

		// Loop over receive map
		size_type i = 0;
		for (auto& [rank, matrix_map] : recv_map_) {
			recv_requests[i++].wait(); // wait for data to appear
			auto& buffer_for_rank = recv_buffer_[rank];

			// Loop over each value received from the rank
			size_type n = 0;
			for(auto& local_id_map : matrix_map){

				// Loop over each local ID the value is used in
				for(auto& id : local_id_map) {
					ASSERT(id < reduction_buffer_.size());
					reduction_buffer_[id].push_back( buffer_for_rank[n] );
					ASSERT(reduction_buffer_[id].size() < MAX_DUPLICATES);
				}
				++n;
			}
		}

		// Call the Reduction Operation on each
		n = 0;
		for(auto& dup_values : reduction_buffer_) {
			local_values[n++] = op(dup_values);
		}

		// Don't exit before all sends complete
		for(auto& req : send_requests){
			req.wait();
		}
	}

private:

	std::map<rank_type, std::vector<size_type>>              send_map_;         // [Rank][1:Nsend] = Local Index
	std::map<rank_type, std::vector<std::vector<size_type>>> recv_map_;         // [Rank][1:Nrecv][1:Nshare] = Local Index
	std::map<rank_type, std::vector<value_type>>             send_buffer_;      // [Rank][1:Nsend] = Value sending
	std::map<rank_type, std::vector<value_type>>             recv_buffer_;      // [Rank][1:Nrecv] = Value received
	std::vector<std::vector<size_type>>                      reduction_buffer_; // [1:Nlocal][1:Nreduce] = Value to reduce
};

#endif /* MRAP_HPP */
