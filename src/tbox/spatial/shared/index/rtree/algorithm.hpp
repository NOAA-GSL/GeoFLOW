/*
 * base.hpp
 *
 *  Created on: Sep 22, 2020
 *      Author: bflynt
 */

#ifndef RTREE_ALGORITHM_BASE_HPP_
#define RTREE_ALGORITHM_BASE_HPP_


#include "tbox/spatial/bound/concept.hpp"


#include <cassert>
#include <memory>
#include <utility>

namespace tbox {
namespace spatial {
namespace shared {
namespace index {
namespace rtree {


template<typename SplittingType>
struct Algorithm {
	static constexpr std::size_t min_children = SplittingType::min_children;
	static constexpr std::size_t max_children = SplittingType::max_children;


	template<typename NodePtr>
	static
	std::pair<NodePtr,NodePtr>
	pick_seeds(const NodePtr& parent_ptr) {
		return SplittingType::pick_seeds(parent_ptr);
	}

	template<typename NodePtr>
	static
	std::pair<NodePtr,NodePtr>
	pick_next(const NodePtr& parent_ptr, const NodePtr& a_node_ptr, const NodePtr& b_node_ptr) {
		return SplittingType::pick_next(parent_ptr,a_node_ptr,b_node_ptr);
	}


	/** Split a single Node into two minimum entry Nodes
	 *
	 *  Function splits the provided parent Node into two
	 *  "optimal" split nodes and returns them as a std::pair.
	 *  The original parent Node is returned minus the
	 *  children which were transfered to the to new Nodes.
	 *
	 *  @param[in] parent_ptr Pointer to the parent node to select a child seeds from
	 *
	 *	@returns Pair of Node pointers to the seeds
	 */
	template<typename NodePtr>
	static
	std::pair<NodePtr,NodePtr>
	split_node(NodePtr& parent_ptr) {
		using node_type = typename NodePtr::element_type;
		assert(parent_ptr);
		assert(parent_ptr->isPage());

		// Find the two optimal seed children
		auto children_seeds = pick_seeds(parent_ptr);

		// Create two new Nodes using the seeds as first entries
		auto a_node_ptr = std::make_shared<node_type>();
		a_node_ptr->insert(std::get<0>(children_seeds));
		parent_ptr->remove(std::get<0>(children_seeds), std::false_type());

		auto b_node_ptr = std::make_shared<node_type>();
		b_node_ptr->insert(std::get<1>(children_seeds));
		parent_ptr->remove(std::get<1>(children_seeds), std::false_type());

		// Loop over all parent children placing them
		while( (parent_ptr->size() > 0)                               and
			   (parent_ptr->size()+a_node_ptr->size() > min_children) and
			   (parent_ptr->size()+b_node_ptr->size() > min_children) ){
			auto node_pair = pick_next(parent_ptr, a_node_ptr, b_node_ptr);
			std::get<1>(node_pair)->insert(std::get<0>(node_pair));
			parent_ptr->remove(std::get<0>(node_pair), std::false_type());
		}

		// Make sure new nodes meet minimum size
		if(parent_ptr->size() > 0){
			if(a_node_ptr->size() < min_children) {
				do {
					a_node_ptr->insert(parent_ptr->front());
					parent_ptr->remove(parent_ptr->front(), std::false_type());
				} while (parent_ptr->size() > 0);
			}
			else {
				do {
					b_node_ptr->insert(parent_ptr->front());
					parent_ptr->remove(parent_ptr->front(), std::false_type());
				} while (parent_ptr->size() > 0);
			}
		}
		parent_ptr->clear();

		assert(a_node_ptr->size() >= min_children);
		assert(b_node_ptr->size() >= min_children);
		assert(a_node_ptr->size() <= max_children);
		assert(b_node_ptr->size() <= max_children);
		assert(parent_ptr->size() == 0);
		return std::make_pair<NodePtr,NodePtr>(std::move(a_node_ptr),std::move(b_node_ptr));
	}

	/**
	 *  Search the current node for best fit
	 *
	 *  Searches the provided node for the best location to insert an object
	 *  with the provided bounding object. An assert will be thrown if the
	 *  current node is invalid.  This function can return a leaf Node
	 *  which is used by the calling function to identify when the bottom of
	 *  the tree has been reached.
	 *
	 *  @param starting_node[in] Pointer to top level Node to start search from
	 *  @param bounding_box[in] Bounding object to use for search
	 *
	 *  @return Pointer to the Node (non-leaf) which is the best geometric fit
	 */
	template<typename BBox, typename NodePtr>
	static NodePtr find_best_fit_in_node(const BBox& bounding_box, const NodePtr& current_node) {
		using spatial::bound::IncreaseToHold;
		assert(current_node);
		assert(current_node->size() > 0);

		// Quick Return if children are leafs
		if( current_node->front()->isLeaf() ) {
			return current_node->front();
		}

		// Search for child with smallest increase to hold
		NodePtr best_node(nullptr);
		auto current_minimum = std::numeric_limits<typename BBox::value_type>::max();

		// Loop over all Children
		for(const NodePtr& child : *current_node){

			// Get area required to hold both
			auto area_to_enlarge = IncreaseToHold(child->getBound(), bounding_box);

			// If the area increase is smaller than accept
			if( area_to_enlarge < current_minimum ){
				current_minimum = area_to_enlarge;
				best_node = child;
			}
			else if (area_to_enlarge == current_minimum) {
				assert(best_node);
				if( (child->area() < best_node->area()) or
					(child->size() < best_node->size())){
					best_node = child;
				}
			}
		}
		assert(best_node);
		return best_node;
	}

	/**
	 *  Search the full tree from the starting node for best fit
	 *
	 *  Searches the down the tree from the provided node for the
	 *  best location to insert an object with the provided bounding object.
	 *  If the starting node is empty then a new node is created as the best
	 *  place to insert the object.
	 *
	 *  @param starting_node[in] Pointer to top level Node to start search from
	 *  @param bounding_box[in] Bounding object to use for search
	 *
	 *  @return Pointer to the Node (non-leaf) which is the best geometric fit
	 */
	template<typename BBox, typename NodePtr>
	static NodePtr find_best_fit_in_tree(const NodePtr& starting_node, const BBox& bounding_box){
		using node_type = typename NodePtr::element_type;

		// If the starting_node is NULL then create a new one for placement
		if(not starting_node) {
			return std::make_shared<node_type>();
		}

		// Search within the Node and descend down tree
		auto current_node = starting_node;
		while(current_node->isPage()){
			assert(current_node->size() > 0 );
			current_node = find_best_fit_in_node(bounding_box,current_node);
		}
		assert(current_node->isLeaf());
		return current_node->getParent();
	}

	/**
	 *  Search the full tree from the starting node for best fit
	 *
	 *  @param starting_node[in] Pointer to top level Node to start search from
	 *  @param bounding_box[in] Bounding object to use for search
	 *
	 *  @return Pointer to the Node (non-leaf) which is the best geometric fit
	 */
	template<typename NodePtr>
	static NodePtr expand_tree(const NodePtr& starting_node){
		using node_type = typename NodePtr::element_type;
		assert(starting_node);

		NodePtr current_node_ptr(starting_node);
		while(current_node_ptr->hasParent()){

			auto& parent_node = current_node_ptr->getParent();

			// Check if it needs to be split
			if( current_node_ptr->size() > max_children ){
				auto split_pair  = split_node(current_node_ptr);
				parent_node->remove(current_node_ptr, std::false_type());
				parent_node->insert(std::get<0>(split_pair));
				parent_node->insert(std::get<1>(split_pair));
			}
			//else {
			//	current_node_ptr->stretch(*starting_node);
			//}
			current_node_ptr = parent_node;
			current_node_ptr->stretch(*starting_node);
		}

		// At this point the current_node_ptr == root pointer
		if( current_node_ptr->size() > max_children ){
			assert(not current_node_ptr->hasParent());

			// Split the current_node_ptr
			auto split_pair = split_node(current_node_ptr);

			// Create new root Node using one of the split Nodes
			auto new_root_ptr = std::make_shared<node_type>();

			// Add in the split nodes
			new_root_ptr->insert(std::get<0>(split_pair));
			new_root_ptr->insert(std::get<1>(split_pair));

			// Assign the new root pointer to be returned
			current_node_ptr = new_root_ptr;
		}
		return current_node_ptr;
	}



	template<typename NodePtr>
	static NodePtr condense_tree(const NodePtr& starting_node){
		using node_type = typename NodePtr::element_type;
		using do_not_restretch = std::false_type;
		assert(starting_node);

		// Get a list for orphan nodes
		std::list<NodePtr> orphan_node_list;

		NodePtr current_node_ptr(starting_node);
		while(current_node_ptr->hasParent()){

			auto& parent_node = current_node_ptr->getParent();

			// Eliminate under sized nodes
			if( current_node_ptr->size() < min_children ){
				std::copy(current_node_ptr->begin(),
						  current_node_ptr->end(),
						  std::back_inserter(orphan_node_list));
				parent_node->remove(current_node_ptr, do_not_restretch());
			}

			current_node_ptr = parent_node;
			current_node_ptr->restretch();
		}

		// Re-Insert all Orphans
		for(NodePtr& orphan : orphan_node_list) {
			current_node_ptr = insert(current_node_ptr,orphan);
		}

		// If the root only has one child then make that the root
		if( (current_node_ptr->size() == 1) and (not current_node_ptr->front()->isLeaf()) ){
			current_node_ptr = current_node_ptr->front();
			current_node_ptr->setParent(nullptr);
		}

		return current_node_ptr;
	}

	/**
	 * Place a new node within the tree under starting_node
	 *
	 * @param starting_node[in] Pointer to the Node to place the inserted Node under
	 * @param place_node[in] Pointer to the Node that is being inserted
	 *
	 * @returns Pointer to new starting Node (if needed to expand)
	 */
	template<typename NodePtr>
	static NodePtr insert(const NodePtr& starting_node, NodePtr& place_node){

		// Find the node to insert into
		// - new pointer if starting_node == nullptr
		// - existing node if starting_node != nullptr
		auto best_node_ptr = find_best_fit_in_tree(starting_node,place_node->getBound());

		// Insert into the best node
		best_node_ptr->insert(place_node);

		// Non-Recursively ascend tree adjusting boxes
		// - assigns into root for case where new root was created
		return expand_tree(best_node_ptr);
	}

	/**
	 * Remove an existing Node within the tree under starting_node
	 *
	 * Removes all nodes which have the same properties (box & data)
	 * as the provided remove_node pointer. It does not attempt to match
	 * the location pointed at by the shared_ptr.
	 *
	 * @param starting_node[in] Pointer to the Node to place the inserted Node under
	 * @param remove_node[in] Pointer to the Node that is properties of Nodes to remove
	 *
	 * @returns Pointer to new starting Node (if needed to shrink)
	 */
	template<typename NodePtr>
	static NodePtr remove(const NodePtr& starting_node, const NodePtr& remove_node){

		/**
		NodePtr new_root(nullptr);

		// Get list of all nodes with a matching bounding_box
		auto matching_list = matches(remove_node->getBound(),starting_node);

		for(NodePtr& match : matching_list){
			if(remove_node->data() == match->data()) {
				match->getParent()->remove(remove_node);
				new_root = condense_tree(match);
			}
		}
		**/

		// Find the node where it should live
		// - new pointer if starting_node == nullptr
		// - existing node if starting_node != nullptr
		auto best_node_ptr = find_best_fit_in_tree(starting_node,remove_node->getBound());
		assert(best_node_ptr);

		// Search all children for matching properties and remove
		std::list<NodePtr> matching_node_list;
		for(NodePtr& child : *best_node_ptr){
			assert(child->isLeaf());
			if( (child->getBound()  == remove_node->getBound()) and
				(child->data() == remove_node->data())  ) {
				matching_node_list.push_front(child);
			}
		}

		// Remove all matching Nodes + ReStretch
		using restretch_node = std::false_type;
		for(NodePtr& match : matching_node_list){
			best_node_ptr->remove(match, restretch_node());
		}
		best_node_ptr->restretch();

		// Non-Recursively ascend tree adjusting boxes
		// - assigns into root for case where new root was eliminated
		return condense_tree(best_node_ptr);
	}


	template<typename NodePtr>
	static void Diagnostics(const NodePtr& starting_node){

		std::size_t leaf_count = 0;
		std::size_t page_count = 0;
		std::array<std::size_t,max_children+1> child_count;
		child_count.fill(0);

		std::list<NodePtr> node_list;
		node_list.push_back(starting_node);
		while(node_list.size() > 0){
			auto child = node_list.front();
			node_list.pop_front();

			if( child->isLeaf() ){
				++leaf_count;
			}
			else {
				++page_count;
				child_count[child->size()] += 1;
				for(auto& grand_child : *child ){
					node_list.push_back(grand_child);
				}
			}
		}

		std::cout << "Diagnostics"            << std::endl;
		std::cout << "Pages = " << page_count << std::endl;
		std::cout << "Leafs = " << leaf_count << std::endl;
		for(auto i = 0; i < child_count.size(); ++i){
			std::cout << "Pages with " << i << " Children = " << child_count[i] << std::endl;
		}

	}





};


} /* namespace rtree */
} /* namespace index */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */


#endif /* RTREE_ALGORITHM_BASE_HPP_ */
