/*
 * nearest.hpp
 *
 *  Created on: Sep 21, 2020
 *      Author: bflynt
 */

#ifndef PREDICATE_NEAREST_HPP_
#define PREDICATE_NEAREST_HPP_

#include "tbox/spatial/shared/predicate/dispatch.hpp"

#include <type_traits>

namespace tbox {
namespace spatial {
namespace shared {
namespace predicate {

/**
 * Resulting type of predicate operation to perform
 */
template<typename BoundType, typename NodeOpTag, typename LeafOpTag>
class distance_predicate final {

public:
	using size_type  = std::size_t;
	using value_type = typename BoundType::value_type;

	distance_predicate() = delete;
	distance_predicate(const distance_predicate& other) = default;
	distance_predicate(distance_predicate&& other) = default;
	~distance_predicate() = default;
	distance_predicate& operator=(const distance_predicate& other) = default;
	distance_predicate& operator=(distance_predicate&& other) = default;

	distance_predicate(BoundType const& bound, size_type const n) : _bound(bound), _count(n) {
	}

	value_type operator()(const BoundType& bound, const std::true_type /* is_leaf */) const noexcept {
		return detail::dispatch<LeafOpTag>::apply(bound,_bound);
	}

	value_type operator()(const BoundType& bound, const std::false_type /* is_leaf */) const noexcept {
		return detail::dispatch<NodeOpTag>::apply(bound,_bound);
	}

	value_type operator()(const BoundType& bound, const bool is_leaf) const noexcept {
		if( is_leaf ) {
			return this->operator()(bound,std::true_type());
		}
		return this->operator()(bound,std::false_type());
	}

	size_type count() const noexcept {
		return _count;
	}

protected:
	BoundType   _bound;
	std::size_t _count;
};


template<typename Predicate>
struct is_distance_predicate : public std::false_type {
};

template<typename BoundType, typename NodeOpTag, typename LeafOpTag>
struct is_distance_predicate<distance_predicate<BoundType,NodeOpTag,LeafOpTag>> : public std::true_type {
};

} /* namespace predicate */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */


#endif /* PREDICATE_NEAREST_HPP_ */
