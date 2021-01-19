/*
 * linear.hpp
 *
 *  Created on: Sep 17, 2020
 *      Author: bflynt
 */

#ifndef ALGORITHM_LINEAR_HPP_
#define ALGORITHM_LINEAR_HPP_


#include "tbox/spatial/bound/concept.hpp"                          // spatial::bound::IncreaseToHold
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
struct Linear final : public Algorithm<Linear<MaxChildren,MinChildren>> {
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
		assert(parent_ptr->size() > 1);

		// Make a returns value
		auto ans = std::make_pair<NodePtr,NodePtr>(nullptr,nullptr);

		// SetUp some values
		auto max_scaled_length = std::numeric_limits<value_type>::lowest();

		// Loop over Each Dimension
		for(std::size_t i = 0; i < bound_type::ndim; ++i) {

			// Loop over each child
			NodePtr minimum_max_child(nullptr);
			NodePtr maximum_min_child(nullptr);
			auto minimum_max_value = std::numeric_limits<value_type>::max();
			auto maximum_min_value = std::numeric_limits<value_type>::lowest();
			for(const NodePtr& child : *parent_ptr){
				if( child->getBound().min(i) > maximum_min_value ){
					maximum_min_value = child->getBound().min(i);
					maximum_min_child = child;
				}
				if( child->getBound().max(i) < minimum_max_value ){
					minimum_max_value = child->getBound().max(i);
					minimum_max_child = child;
				}
			}

			auto scaled_length = std::abs(minimum_max_value - maximum_min_value)/parent_ptr->getBound().length(i);

			if( scaled_length > max_scaled_length ){
				max_scaled_length = scaled_length;
				ans = std::pair<NodePtr,NodePtr>(minimum_max_child,maximum_min_child);
			}
		}

		// Make sure they are not the same
		// - Seems like the algorithm from the original paper has a flaw
		// - I sometimes get first and second equal to each other when pages are Cartesian aligned cells
		if(ans.first == ans.second) {
			ans.second = (ans.first == parent_ptr->front()) ? parent_ptr->back(): parent_ptr->front();
		}

		assert(ans.first);
		assert(ans.second);
		assert(ans.first != ans.second);
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

		// Select first child as the next
		NodePtr next_child(parent_ptr->front());
		NodePtr place_into(nullptr);

		// Pick smallest increase to go into
		auto a_increased_area = IncreaseToHold(a_node_ptr->getBound(), next_child->getBound());
		auto b_increased_area = IncreaseToHold(b_node_ptr->getBound(), next_child->getBound());
		place_into = (a_increased_area < b_increased_area) ? a_node_ptr: b_node_ptr;

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

#endif /* ALGORITHM_LINEAR_HPP_ */
