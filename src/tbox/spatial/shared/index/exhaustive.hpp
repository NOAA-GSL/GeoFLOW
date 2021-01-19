/*
 * exhaustive.hpp
 *
 *  Created on: Sep 18, 2020
 *      Author: bflynt
 */

#ifndef EXHAUSTIVE_HPP_
#define EXHAUSTIVE_HPP_


#include "tbox/spatial/common/truncated_multiset.hpp"
#include "tbox/spatial/shared/predicate/all.hpp"


#include <algorithm>  // std::remove_if
#include <functional> // std::equal_to
#include <list>       // std::list
#include <memory>     // std::allocator

namespace tbox {
namespace spatial {
namespace shared {
namespace index {


template<typename Value,
		 typename BoundExtractor,
		 typename EqualOp     = std::equal_to<Value>,
		 typename Allocator   = std::allocator<Value>
>
class Exhaustive final {


	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
private:
	using list_type       = std::list<Value,Allocator>;

public:

	using value_type      = typename list_type::value_type;
	using size_type       = typename list_type::size_type;
	using reference       = typename list_type::reference;
	using const_reference = typename list_type::const_reference;
	using allocator_type  = typename list_type::allocator_type;
	using iterator        = typename list_type::iterator;
	using const_iterator  = typename list_type::const_iterator;

	using equal_operator  = EqualOp;
	using bound_extractor = BoundExtractor;
	using bound_type      = typename bound_extractor::bound_type;
	using bound_value_type = typename bound_type::value_type;


	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	Exhaustive() = default;

	Exhaustive(const Exhaustive& other) = default;

	Exhaustive(Exhaustive&& other) = default;

	~Exhaustive() = default;

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	Exhaustive& operator=(const Exhaustive& box) = default;

	Exhaustive& operator=(Exhaustive&& box) = default;

	//-------------------------------------------------------------------------
	// Modifiers
	//-------------------------------------------------------------------------

	void insert(value_type const& value) {
		this->_stretch(value);
		this->_values.push_back(value);
	}

	template<typename Iterator>
	void insert(Iterator first, Iterator last) {
		while(first != last) {
			this->insert(*first);
			++first;
		}
	}

	void remove(value_type const& value) {
		this->_remove(value);
		this->_restretch();
	}

	template<typename Iterator>
	void remove(Iterator first, Iterator last) {
		this->_bulk_remove(first,last);
		this->_restretch();
	}

	void clear() const noexcept {
		this->_values.clear();
		this->_bound.reset();
	}

	//-------------------------------------------------------------------------
	// Iterators
	//-------------------------------------------------------------------------

	iterator begin() noexcept {
		return this->_values.begin();
	}

	const_iterator begin() const noexcept {
		return this->_values.begin();
	}

	iterator end() noexcept {
		return this->_values.end();
	}

	const_iterator end() const noexcept {
		return this->_values.end();
	}

	//-------------------------------------------------------------------------
	// Capacity
	//-------------------------------------------------------------------------

	size_type size() const noexcept {
		return this->_values.size();
	}

	bool empty() const noexcept {
		return this->_values.empty();
	}

	size_type max_size() const noexcept {
		return this->_values.max_size();
	}

	//-------------------------------------------------------------------------
	// Indexing
	//-------------------------------------------------------------------------

	bound_type const& bounds() const noexcept {
		return this->_bound;
	}


	template<typename Predicates, typename OutIter>
	size_type query(Predicates const& pred, OutIter out_it) const noexcept {
		return _query_dispatch(pred,out_it, tbox::spatial::shared::predicate::is_distance_predicate<Predicates>());
	}


	//-------------------------------------------------------------------------
	// Data [Private]
	//-------------------------------------------------------------------------
private:
	list_type  _values;  // List of value stored within the Index
	bound_type _bound;   // Bound of Index


	//-------------------------------------------------------------------------
	// Internal Functions [Private]
	//-------------------------------------------------------------------------
private:

	bound_type const& _extract_bound(value_type const& value) const noexcept {
		static const bound_extractor extract;
		return extract(value);
	}

	template<typename Iterator>
	void _remove(value_type const& value) {
		std::remove_if(_values.begin(), _values.end(), [&](auto& v){
			return equal_operator()(value,v);
		});
	}

	template<typename Iterator>
	void _bulk_remove(Iterator first, Iterator last) {
		while(first != last) {
			this->_remove(*first);
			++first;
		}
	}

	void _stretch(value_type const& value) {
		_bound.stretch( _extract_bound(value) );
	}

	void _restretch() {
		this->_bound.reset();
		for(auto it = _values.cbegin(); it != _values.cend(); ++it){
			this->_bound.stretch( _extract_bound(*it));
		}
	}

	template<typename Predicates, typename OutIter>
	size_type _query_dispatch(Predicates const& pred, OutIter out_it, const std::false_type& /* is_distance_pred */) const noexcept {
		constexpr auto bound_from_leaf = std::true_type();

		// Loop over every single value and count if predicate passes
		size_type count = 0;
		for(auto& value : this->_values){
			if( pred(_extract_bound(value),bound_from_leaf)) {
				*out_it = value;
				++out_it;
				++count;
			}
		}
		return count;
	}

	template<typename Predicates, typename OutIter>
	size_type _query_dispatch(Predicates const& pred, OutIter out_it, const std::true_type& /* is_distance_pred */) const noexcept {
		constexpr auto bound_from_leaf = std::true_type();

		// Get a multiset which ranks and automatically trims excess
		using dist_value_pair = std::pair<bound_value_type, value_type>;
		using pair_less_op    = LessPair<0>;
		TruncatedMultiSet<dist_value_pair,pair_less_op> min_set(pred.count());

		// Place all evaluated predicates into MultSet which ranks and limits size
		for(auto& value : this->_values){
			min_set.emplace(pred(_extract_bound(value),bound_from_leaf),value);
		}

		// Copy results into provided output iterator
		std::transform(std::cbegin(min_set), std::cend(min_set), out_it, [](auto& value_pair){
			return value_pair.second;
		});

		return pred.count(); //min_set.size();
	}


};



} /* namespace index */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */

#endif /* EXHAUSTIVE_HPP_ */
