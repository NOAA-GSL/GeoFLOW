/*
 * page.hpp
 *
 *  Created on: Sep 22, 2020
 *      Author: bflynt
 */

#ifndef RTREE_PAGE_HPP_
#define RTREE_PAGE_HPP_


#include <list>
#include <memory>

namespace tbox {
namespace spatial {
namespace shared {
namespace index {
namespace rtree {

//
// Declare Node since a Page holds a list of them
//
template<typename Value, typename BoundExtractor>
class Node;





template<typename Value, typename BoundExtractor>
class Page {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
private:
	using self_type        = Page<Value,BoundExtractor>;
	using node_type        = Node<Value,BoundExtractor>;
	using node_pointer     = std::shared_ptr<node_type>;
	using list_type        = std::list<node_pointer>;

public:
	using value_type       = Value;
	using size_type        = typename list_type::size_type;
	using iterator         = typename list_type::iterator;
	using const_iterator   = typename list_type::const_iterator;
	using reference        = typename list_type::reference;
	using const_reference  = typename list_type::const_reference;

	using bound_extractor  = BoundExtractor;
	using bound_type       = typename bound_extractor::bound_type;
	using bound_value_type = typename bound_type::value_type;


	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	Page() = default;

	Page(const Page& other) = default;

	Page(Page&& other) = default;

	~Page() = default;

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	Page& operator=(const Page& box) = default;

	Page& operator=(Page&& box) = default;

	//-------------------------------------------------------------------------
	// Modifiers
	//-------------------------------------------------------------------------

	void clear() noexcept {
		this->nodes_.clear();
	}

	void stretch(bound_type const& other_bound) noexcept {
		this->bound_.stretch(other_bound);
	}

	void restretch() {
		this->bound_.reset();
		for(auto& node_ptr : this->nodes_) {
			this->stretch(node_ptr->getBound());
		}
	}

	/**
	 * Add a child Node to this Node
	 */
	void insert(node_pointer& child_ptr) {
		this->nodes_.push_back(child_ptr);
		this->stretch(child_ptr->getBound());
	}

	/**
	 * Remove a child Node to from this Node
	 */
	void remove(node_pointer& child_ptr, const bool re_stretch = true) {
		if(re_stretch) {
			this->remove(child_ptr, std::true_type());
		}
		this->remove(child_ptr, std::false_type());
	}

	void remove(node_pointer& child_ptr, std::true_type /* re_stretch */) {
		this->nodes_.remove(child_ptr);
		this->restretch();
	}

	void remove(node_pointer& child_ptr, std::false_type /* re_stretch */) {
		this->nodes_.remove(child_ptr);
	}


	//-------------------------------------------------------------------------
	// Element Access
	//-------------------------------------------------------------------------

	reference front() {
		return nodes_.front();
	}

	const_reference front() const {
		return nodes_.front();
	}

	reference back() {
		return nodes_.back();
	}

	const_reference back() const {
		return nodes_.back();
	}

	//-------------------------------------------------------------------------
	// Iterators
	//-------------------------------------------------------------------------

	iterator begin() noexcept {
		return nodes_.begin();
	}

	const_iterator begin() const noexcept {
		return nodes_.begin();
	}

	const_iterator cbegin() const noexcept {
		return nodes_.cbegin();
	}

	iterator end() noexcept {
		return nodes_.end();
	}

	const_iterator end() const noexcept {
		return nodes_.end();
	}

	const_iterator cend() const noexcept {
		return nodes_.cend();
	}

	iterator rbegin() noexcept {
		return nodes_.rbegin();
	}

	const_iterator rbegin() const noexcept {
		return nodes_.rbegin();
	}

	const_iterator crbegin() const noexcept {
		return nodes_.crbegin();
	}

	iterator rend() noexcept {
		return nodes_.rend();
	}

	const_iterator rend() const noexcept {
		return nodes_.rend();
	}

	const_iterator crend() const noexcept {
		return nodes_.crend();
	}


	//-------------------------------------------------------------------------
	// Capacity
	//-------------------------------------------------------------------------

	bool empty() const noexcept {
		return this->nodes_.empty();
	}

	size_type size() const noexcept {
		return this->nodes_.size();
	}

	size_type max_size() const noexcept {
		return this->nodes_.max_size();
	}

	//-------------------------------------------------------------------------
	// Indexing
	//-------------------------------------------------------------------------


	const bound_type& getBound() const noexcept {
		return bound_;
	}

	//-------------------------------------------------------------------------
	// Data [Private]
	//-------------------------------------------------------------------------
private:
	list_type  nodes_;
	bound_type bound_;

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




#endif /* RTREE_PAGE_HPP_ */
