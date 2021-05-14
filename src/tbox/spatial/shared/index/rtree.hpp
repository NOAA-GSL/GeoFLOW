/*
 * rtree.hpp
 *
 *  Created on: Sep 18, 2020
 *      Author: bflynt
 */

#ifndef RTREE_HPP_
#define RTREE_HPP_

#include "tbox/spatial/shared/index/rtree/algorithm.hpp"
#include "tbox/spatial/shared/index/rtree/leaf.hpp"
#include "tbox/spatial/shared/index/rtree/linear.hpp"
#include "tbox/spatial/shared/index/rtree/node.hpp"
#include "tbox/spatial/shared/index/rtree/page.hpp"
#include "tbox/spatial/shared/index/rtree/quadratic.hpp"

#include "tbox/spatial/common/truncated_multiset.hpp"
#include "tbox/spatial/shared/predicate/all.hpp"

#include <algorithm>  // std::remove_if
#include <functional> // std::equal_to
#include <memory>     // std::allocator

namespace tbox {
namespace spatial {
namespace shared {
namespace index {



template<typename Value,
		 typename BoundGetter,
         typename Parameters  = rtree::Quadratic<10,4>,
		 typename EqualOp     = std::equal_to<Value>,
		 typename Allocator   = std::allocator<Value>
>
class RTree final {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
private:
	using node_type            = rtree::Node<Value,BoundGetter>;
	using node_reference       = node_type&;
	using const_node_reference = node_type const&;
	using node_pointer         = std::shared_ptr<node_type>;
	using Algorithm            = Parameters;

public:

	using value_type       = Value;
	using size_type        = std::size_t;
	using reference        = value_type&;
	using const_reference  = value_type const&;
	using allocator_type   = Allocator;
	//using iterator         = typename list_type::iterator;
	//using const_iterator   = typename list_type::const_iterator;

	using parameters       = Parameters;
	using equal_operator   = EqualOp;
	using bound_extractor  = BoundGetter;
	using bound_type       = typename bound_extractor::bound_type;
	using bound_value_type = typename bound_type::value_type;



	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	RTree() = default;

	RTree(const RTree& other) = default;

	RTree(RTree&& other) = default;

	~RTree() = default;

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	RTree& operator=(const RTree& box) = default;

	RTree& operator=(RTree&& box) = default;


	//-------------------------------------------------------------------------
	// Modifiers
	//-------------------------------------------------------------------------

	void insert(value_type const& value) {
		auto new_leaf = std::make_shared<node_type>(value);
		root_node_ptr_ = Algorithm::insert(root_node_ptr_,new_leaf);
	}

	template<typename Iterator>
	void insert(Iterator first, Iterator last) {
		while(first != last) {
			this->insert(*first);
			++first;
		};
	}

	void remove(value_type const& value) {
		auto del_leaf = std::make_shared<node_type>(value);
		root_node_ptr_ = Algorithm::remove(root_node_ptr_,del_leaf);
	}

	template<typename Iterator>
	void remove(Iterator first, Iterator last) {
		while(first != last) {
			this->remove(*first);
			++first;
		};
	}

	void clear() const noexcept {
		root_node_ptr_.reset();
	}

	//-------------------------------------------------------------------------
	// Iterators
	//-------------------------------------------------------------------------


	//-------------------------------------------------------------------------
	// Capacity
	//-------------------------------------------------------------------------

	size_type size() const noexcept {
		return 0;
	}

	bool empty() const noexcept {
		return (root_node_ptr_);
	}

	constexpr size_type max_size() const noexcept {
		return std::numeric_limits<size_type>::max();
	}

	//-------------------------------------------------------------------------
	// Indexing
	//-------------------------------------------------------------------------

	bound_type const& bounds() const noexcept {
		return root_node_ptr_->getBound();
	}

	template<typename Predicates, typename OutIter>
	size_type query(Predicates const& pred, OutIter out_it) const noexcept {
		return query_dispatch(pred,out_it, predicate::is_distance_predicate<Predicates>());
	}


	void display() const {
		Algorithm::Diagnostics(root_node_ptr_);
	}


	//-------------------------------------------------------------------------
	// Data [Private]
	//-------------------------------------------------------------------------
private:
	node_pointer root_node_ptr_;


	//-------------------------------------------------------------------------
	// Internal Functions [Private]
	//-------------------------------------------------------------------------
private:

	template<typename Predicates, typename OutIter>
	size_type query_dispatch(Predicates const& pred, OutIter out_it, const std::false_type& /* is_distance_pred */) const noexcept {

		// Quick return if nothing to be done
		if(!root_node_ptr_) return 0;

		// Get stack for searching
		std::list<node_pointer> candidate_nodes;

		// Iterate over stack till all pages are processed
		size_type count = 0;
		candidate_nodes.push_back(root_node_ptr_);
		while(candidate_nodes.size() > 0){

			auto const& current_candidate = candidate_nodes.front();
			auto candidate_is_leaf = current_candidate->isLeaf();

			if( pred(current_candidate->getBound(), candidate_is_leaf) ) {
				if( candidate_is_leaf ){
					*out_it = current_candidate->getValue();
					++out_it;
					++count;
				}
				else {
					for(auto& child : *current_candidate){
						candidate_nodes.push_back(child);
					}
				}
			}
			candidate_nodes.pop_front();
		}
		return count;
	}


	template<typename Predicates, typename OutIter>
	size_type query_dispatch(Predicates const& pred, OutIter out_it, const std::true_type& /* is_distance_pred */) const noexcept {
		constexpr auto bound_from_leaf = std::true_type();
		constexpr auto bound_from_page = std::false_type();

		// Quick return if nothing to be done
		if(!root_node_ptr_) return 0;

		// Get 2 truncated stacks for searching
		using pair_less_op       = LessPair<0>;
		using dist_node_ptr_pair = std::pair<bound_value_type, node_pointer>;
		std::multiset<dist_node_ptr_pair,pair_less_op>     candidate_nodes;
		TruncatedMultiSet<dist_node_ptr_pair,pair_less_op> candidate_leafs(pred.count());

		// Insert the Root Node into the candidate_nodes
		auto distance_threshhold = std::numeric_limits<bound_value_type>::max();
		auto current_candidate   = node_pointer(nullptr);
		auto is_leaf             = root_node_ptr_->isLeaf();
		auto dist                = pred(root_node_ptr_->getBound(), is_leaf);
		candidate_nodes.emplace(dist,root_node_ptr_);

		// Iterate over stack till all pages are processed
		while(candidate_nodes.size() > 0){

			auto front        = candidate_nodes.begin();
			dist              = front->first;
			current_candidate = front->second;
			candidate_nodes.erase(front);

			// Ensure the current_candidate is still a possible candidate
			if( dist <= distance_threshhold ) {

				// If Leaf
				// - Place into the Leaf set
				// - If we have K leafs then update the new distance tolerance
				if( current_candidate->isLeaf() ) {
					candidate_leafs.emplace(dist,current_candidate);
					if(candidate_leafs.size() >= pred.count()) {
						distance_threshhold = candidate_leafs.crbegin()->first;
					}
				}
				else {
					// Loop over children into the candidates
					for(auto& child : *current_candidate){
						dist = pred(child->getBound(), child->isLeaf());
						candidate_nodes.emplace(dist,child);
					}
				}
			}
		}

		// Copy resulting values into provided output iterator
		std::transform(std::cbegin(candidate_leafs), std::cend(candidate_leafs), out_it, [](auto& value_pair){
			return value_pair.second->getValue();
		});

		return candidate_leafs.size();
	}
};

} /* namespace index */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */


#endif /* RTREE_HPP_ */
