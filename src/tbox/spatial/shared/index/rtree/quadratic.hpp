/*
 * quadratic.hpp
 *
 *  Created on: Sep 17, 2020
 *      Author: bflynt
 */

#ifndef ALGORITHM_QUADRATIC_HPP_
#define ALGORITHM_QUADRATIC_HPP_

#include "tbox/spatial/bound/concept.hpp"  // spatial::bound::IncreaseToHold
#include "tbox/spatial/shared/index/rtree/algorithm.hpp"


#include <cassert>
#include <iterator>
#include <limits>
#include <memory>
#include <utility>

namespace tbox {
namespace spatial {
namespace shared {
namespace index {
namespace rtree {


template<std::size_t MaxChildren, std::size_t MinChildren = MaxChildren/2>
struct Quadratic final : public Algorithm<Quadratic<MaxChildren,MinChildren>> {
	static_assert((MinChildren > 1) && (MinChildren <= MaxChildren/2));

	static constexpr std::size_t max_children = MaxChildren;
	static constexpr std::size_t min_children = MinChildren;

	/** Choose two seed children from provided Node pointer
	 *
	 *  Function selects two child Nodes from within the node_ptr
	 *  to used as seed Nodes for a the start of Node splitting.
	 *  Returns both Nodes inside a std::pair.  Provided parent
	 *  Node is not modified in any way.
	 *
	 *  @param[in] parent_ptr Pointer to the parent node to select a child seeds from
	 *
	 *	@returns Pair of Node pointers to the seeds
	 */
	template<typename NodePtr>
	static
	std::pair<NodePtr,NodePtr>
	pick_seeds(const NodePtr& parent_ptr) {
		using node_type  = typename NodePtr::element_type;
		using bound_type = typename node_type::bound_type;
		using value_type = typename bound_type::value_type;
		assert(parent_ptr);

		// Make a returns value
		auto ans = std::make_pair<NodePtr,NodePtr>(nullptr,nullptr);

		// SetUp some values
		auto max_wasted_area = std::numeric_limits<value_type>::lowest();

		// Loop over N*log(N) iterations
		for(auto it = parent_ptr->begin(); it != parent_ptr->end(); ++it){
			auto i_area = (*it)->area();
			for(auto jt = std::next(it);  jt != parent_ptr->end(); ++jt){

				// Calculate wasted area
				auto j_area = (*jt)->area();
				auto c_area = Union((*it)->getBound(), (*jt)->getBound()).area();
				auto wasted_area = c_area - i_area - j_area;

				// If the packing is worst hold pointers to Nodes
				if( wasted_area > max_wasted_area ){
					max_wasted_area = wasted_area;
					ans = std::pair<NodePtr,NodePtr>(*it,*jt);
				}
			}
		}
		assert(ans.first);
		assert(ans.second);
		return ans;
	}


	/** Find the next child to place and where to place it
	 *
	 *  Function selects a child Node from within the parent_ptr
	 *  to be placed within either A or B node.  Returns a pair of
	 *  Node points with the frist entry being the child node to place
	 *  and the second being the Node to place it into.
	 *
	 *  @param[in] parent_ptr Pointer to the parent node to select a child from
	 *  @param[in] a_node_ptr Pointer to Node candidate A
	 *  @param[in] b_node_ptr Pointer to Node candidate B
	 *
	 *	@returns Pair of Node pointers
	 */
	template<typename NodePtr>
	static
	std::pair<NodePtr,NodePtr>
	pick_next(const NodePtr& parent_ptr, const NodePtr& a_node_ptr, const NodePtr& b_node_ptr) {
		using node_type  = typename NodePtr::element_type;
		using bound_type = typename node_type::bound_type;
		using value_type = typename bound_type::value_type;
		using spatial::bound::IncreaseToHold;
		assert(parent_ptr);
		assert(a_node_ptr);
		assert(b_node_ptr);

		// Pointers which will hold the results
		NodePtr next_child(nullptr);
		NodePtr place_into(nullptr);

		// SetUp some values
		auto max_increase = std::numeric_limits<value_type>::lowest();

		// Loop over all children finding the next one to place
		for(const NodePtr& child : *parent_ptr){
			auto a_increased_area = IncreaseToHold(a_node_ptr->getBound(), child->getBound());
			auto b_increased_area = IncreaseToHold(b_node_ptr->getBound(), child->getBound());
			auto diff = std::abs(a_increased_area-b_increased_area);

			if( diff > max_increase ){
				max_increase = diff;
				next_child   = child;
				place_into   = (a_increased_area < b_increased_area) ? a_node_ptr: b_node_ptr;
			}
		}
		assert(next_child);
		assert(place_into);
		return std::make_pair<NodePtr,NodePtr>(std::move(next_child),std::move(place_into));
	}

};



} /* namespace rtree */
} /* namespace index */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */


#endif /* ALGORITHM_QUADRATIC_HPP_ */
