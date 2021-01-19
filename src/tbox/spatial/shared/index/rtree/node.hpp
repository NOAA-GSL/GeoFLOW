/*
 * node.hpp
 *
 *  Created on: Sep 22, 2020
 *      Author: bflynt
 */

#ifndef RTREE_NODE_HPP_
#define RTREE_NODE_HPP_

#include "tbox/spatial/shared/index/rtree/page.hpp"
#include "tbox/spatial/shared/index/rtree/leaf.hpp"

#include <cassert>
#include <memory>
#include <variant>

namespace tbox {
namespace spatial {
namespace shared {
namespace index {
namespace rtree {


template<typename Value, typename BoundExtractor>
class Node final : public std::enable_shared_from_this<Node<Value,BoundExtractor>> {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
private:
	using self_type       = Node<Value,BoundExtractor>;
	using leaf_type       = Leaf<Value,BoundExtractor>;
	using page_type       = Page<Value,BoundExtractor>;
	using variant_type    = std::variant<page_type,leaf_type>;
	using bound_extractor = BoundExtractor;

public:
	using value_type           = Value;
	using size_type            = typename page_type::size_type;
	using bound_type           = typename bound_extractor::bound_type;
	using bound_value_type     = typename bound_type::value_type;
	using node_pointer         = std::shared_ptr<self_type>;
	using child_iterator       = typename page_type::iterator;
	using const_child_iterator = typename page_type::const_iterator;


	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	//Node() = default;
	Node() : parent_ptr_(nullptr) {
		std::get<page_type>(data_).restretch();
	}

	Node(const Node& other) = default;

	Node(Node&& other) = default;

	~Node() = default;

	//Node(node_pointer& child_a_ptr, node_pointer& child_b_ptr) : parent_ptr_(nullptr) {
	//	std::get<page_type>(data_).getBound().reset();
	//	this->insert(child_a_ptr);
	//	this->insert(child_b_ptr);
	//}

	Node(value_type const& value) : parent_ptr_(nullptr), data_(value) {
	}


	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	Node& operator=(const Node& box) = default;

	Node& operator=(Node&& box) = default;

	//-------------------------------------------------------------------------
	// Modifiers
	//-------------------------------------------------------------------------

	void clear() {
		assert(this->isPage());
		std::get<page_type>(data_).clear();
	}

	/**
	 * Assign a parent node to this
	 *
	 * @param[in] parent The Node to assign as the parent of this Node
	 */
	void setParent(node_pointer const& parent) noexcept {
		parent_ptr_ = parent;
	}

	/**
	 * Add a child Node to this Node
	 */
	void insert(node_pointer& child_ptr) {
		assert(this->isPage());

		// Set myself as child parent
		child_ptr->setParent(this->shared_from_this());

		// Insert this Child to my page
		std::get<page_type>(data_).insert(child_ptr);
	}

	/**
	 * Remove a child Node to from this Node
	 */
	void remove(node_pointer& child_ptr, const bool re_stretch = true) {
		assert(this->isPage());
		if(re_stretch) {
			std::get<page_type>(data_).remove(child_ptr, std::true_type());
		}
		std::get<page_type>(data_).remove(child_ptr, std::false_type());
	}

	void remove(node_pointer& child_ptr, std::true_type /* re-stretch */) {
		assert(this->isPage());
		std::get<page_type>(data_).remove(child_ptr, std::true_type());
	}

	void remove(node_pointer& child_ptr, std::false_type /* re-stretch */) {
		assert(this->isPage());
		std::get<page_type>(data_).remove(child_ptr, std::false_type());
	}


	void stretch(self_type const& other) noexcept {
		assert(this->isPage());
		std::get<page_type>(data_).stretch(other.getBound());
	}


	void restretch() noexcept {
		assert(this->isPage());
		std::get<page_type>(data_).restretch();
	}

	//-------------------------------------------------------------------------
	// Element Access
	//-------------------------------------------------------------------------

	/** Test if Node has a Parent
	 */
	bool hasParent() const noexcept {
		return (parent_ptr_.get());
	}

	/**
	 *	Get the parent of this Node
	 */
	const node_pointer getParent() const noexcept {
		return parent_ptr_;
	}

	value_type const& getValue() const noexcept {
		assert(this->isLeaf());
		return std::get<leaf_type>(data_).getValue();
	}

	node_pointer& front() {
		assert(this->isPage());
		return std::get<page_type>(data_).front();
	}

	const node_pointer& front() const {
		assert(this->isPage());
		return std::get<page_type>(data_).front();
	}

	node_pointer& back() {
		assert(this->isPage());
		return std::get<page_type>(data_).back();
	}

	const node_pointer& back() const {
		assert(this->isPage());
		return std::get<page_type>(data_).back();
	}

	//-------------------------------------------------------------------------
	// Iterators
	//-------------------------------------------------------------------------


	child_iterator begin() {
		assert(this->isPage());
		return std::get<page_type>(data_).begin();
	}

	const_child_iterator begin() const {
		assert(this->isPage());
		return std::get<page_type>(data_).begin();
	}

	const_child_iterator cbegin() const {
		assert(this->isPage());
		return std::get<page_type>(data_).cbegin();
	}

	child_iterator end() {
		assert(this->isPage());
		return std::get<page_type>(data_).end();
	}

	const_child_iterator end() const {
		assert(this->isPage());
		return std::get<page_type>(data_).end();
	}

	const_child_iterator cend() const {
		assert(this->isPage());
		return std::get<page_type>(data_).cend();
	}


	//-------------------------------------------------------------------------
	// Capacity
	//-------------------------------------------------------------------------

	bool isLeaf() const {
		return std::holds_alternative<leaf_type>(data_);
	}

	bool isPage() const {
		return std::holds_alternative<page_type>(data_);
	}

	bool empty() const noexcept {
		assert(this->isPage());
		return std::get<page_type>(data_).empty();
	}

	size_type size() const noexcept {
		//assert(this->isPage());
		//return std::get<page_type>(data_).size();

		if(this->isLeaf()){
			return 0;
		}
		return std::get<page_type>(data_).size();
	}

	size_type max_size() const noexcept {
		assert(this->isPage());
		return std::get<page_type>(data_).max_size();
	}

	bound_value_type area() const noexcept {
		return this->getBound().area();
	}

	//-------------------------------------------------------------------------
	// Indexing
	//-------------------------------------------------------------------------

	const bound_type& getBound() const noexcept {
		if(this->isLeaf()){
			return std::get<leaf_type>(data_).getBound();
		}
		return std::get<page_type>(data_).getBound();
	}

	//-------------------------------------------------------------------------
	// Data [Private]
	//-------------------------------------------------------------------------
private:
	variant_type data_;
	node_pointer parent_ptr_;

	//-------------------------------------------------------------------------
	// Internal Functions [Private]
	//-------------------------------------------------------------------------
private:



};




} /* namespace rtree */
} /* namespace index */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */



#endif /* RTREE_NODE_HPP_ */
