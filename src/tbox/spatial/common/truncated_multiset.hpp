/*
 * trucated_multiset.hpp
 *
 *  Created on: Sep 21, 2020
 *      Author: bflynt
 */

#ifndef COMMON_TRUCATED_SET_HPP_
#define COMMON_TRUCATED_SET_HPP_


#include <iterator>
#include <memory>
#include <set>


template<std::size_t N>
struct LessPair {
	template<typename T1, typename T2>
	bool operator()(const std::pair<T1,T2>& p1, const std::pair<T1,T2>& p2) const {
		return std::get<N>(p1) < std::get<N>(p2);
	}
};




template<typename Key,
		 typename Compare   = std::less<Key>,
		 typename Allocator = std::allocator<Key>>
class TruncatedMultiSet {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
private:
	using set_type       = std::multiset<Key,Compare,Allocator>;

public:
	using value_type              = typename set_type::value_type;
	using size_type               = typename set_type::size_type;
	using difference_type         = typename set_type::difference_type;
	using reference               = typename set_type::reference;
	using const_reference         = typename set_type::const_reference;
	using allocator_type          = typename set_type::allocator_type;
	using iterator                = typename set_type::iterator;
	using const_iterator          = typename set_type::const_iterator;
	using reverse_iterator        = typename set_type::reverse_iterator;
	using const_reverse_iterator  = typename set_type::const_reverse_iterator;

	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	TruncatedMultiSet() : _count(std::numeric_limits<difference_type>::max()){
	}

	TruncatedMultiSet(TruncatedMultiSet const& other) = default;
	TruncatedMultiSet(TruncatedMultiSet&& other)      = default;

	TruncatedMultiSet(const size_type count) : _count(count){
	}

	~TruncatedMultiSet() = default;

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	TruncatedMultiSet& operator=(const TruncatedMultiSet& box) = default;

	TruncatedMultiSet& operator=(TruncatedMultiSet&& box) = default;

	//-------------------------------------------------------------------------
	// Iterators
	//-------------------------------------------------------------------------

	iterator begin() noexcept {
		return _data_set.begin();
	}

	const_iterator begin() const noexcept {
		return _data_set.begin();
	}

	const_iterator cbegin() const noexcept {
		return _data_set.cbegin();
	}

	iterator end() noexcept {
		return _data_set.end();
	}

	const_iterator end() const noexcept {
		return _data_set.end();
	}

	const_iterator cend() const noexcept {
		return _data_set.cend();
	}

	reverse_iterator rbegin() noexcept {
		return _data_set.rbegin();
	}

	const_reverse_iterator rbegin() const noexcept {
		return _data_set.rbegin();
	}

	const_reverse_iterator crbegin() const noexcept {
		return _data_set.crbegin();
	}

	reverse_iterator rend() noexcept {
		return _data_set.rend();
	}

	const_reverse_iterator rend() const noexcept {
		return _data_set.rend();
	}

	const_reverse_iterator crend() const noexcept {
		return _data_set.crend();
	}

	//-------------------------------------------------------------------------
	// Capacity
	//-------------------------------------------------------------------------

	bool empty() const noexcept {
		return this->_data_set.empty();
	}

	size_type size() const noexcept {
		return this->_data_set.size();
	}

	size_type max_size() const noexcept {
		return this->_data_set.max_size();
	}

	//-------------------------------------------------------------------------
	// Modifiers
	//-------------------------------------------------------------------------

	void clear() const noexcept {
		this->_data_set.clear();
	}

	void insert(value_type&& value) {
		this->_data_set.insert(std::move(value));
		this->truncate_();
	}

	void insert(const value_type& value){
		this->_data_set.insert(value);
		this->truncate_();
	}

	template<typename InputIt>
	void insert(InputIt first, InputIt last) {
		this->_data_set.insert(first,last);
		this->truncate_();
	}

	template< class... Args >
	void emplace(Args&&... args) {
		this->_data_set.emplace(std::move(args)...);
		this->truncate_();
	}

	iterator erase(const_iterator pos) {
		return this->_data_set.erase(pos);
	}

	iterator erase(const_iterator first, const_iterator last){
		return this->_data_set.erase(first,last);
	}

	size_type erase(const value_type& key) {
		return this->_data_set.erase(key);
	}

private:
	set_type  _data_set;
	size_type _count;


	void truncate_() {
		if(_data_set.size() > _count){
			this->_data_set.erase(std::next(_data_set.begin(),_count),_data_set.end());
		}
	}
};


#endif /* COMMON_TRUCATED_SET_HPP_ */
