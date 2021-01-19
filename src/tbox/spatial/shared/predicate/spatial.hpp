/*
 * spatial.hpp
 *
 *  Created on: Sep 18, 2020
 *      Author: bflynt
 */

#ifndef PREDICATE_SPATIAL_HPP_
#define PREDICATE_SPATIAL_HPP_

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
class spatial_predicate final {


public:
	spatial_predicate() = delete;
	spatial_predicate(const spatial_predicate& other) = default;
	spatial_predicate(spatial_predicate&& other) = default;
	~spatial_predicate() = default;
	spatial_predicate& operator=(const spatial_predicate& other) = default;
	spatial_predicate& operator=(spatial_predicate&& other) = default;

	spatial_predicate(BoundType const& bound) : _bound(bound) {
	}

	bool operator()(const BoundType& bound, const std::true_type /* is_leaf */) const noexcept {
		return detail::dispatch<LeafOpTag>::apply(bound,_bound);
	}

	bool operator()(const BoundType& bound, const std::false_type /* is_leaf */) const noexcept {
		return detail::dispatch<NodeOpTag>::apply(bound,_bound);
	}

	bool operator()(const BoundType& bound, const bool is_leaf) const noexcept {
		if( is_leaf ) {
			return this->operator()(bound,std::true_type());
		}
		return this->operator()(bound,std::false_type());
	}


protected:
	BoundType _bound;
};


template<typename Predicate>
struct is_spatial_predicate : public std::false_type {
};

template<typename BoundType, typename NodeOpTag, typename LeafOpTag>
struct is_spatial_predicate<spatial_predicate<BoundType,NodeOpTag,LeafOpTag>> : public std::true_type {
};


} /* namespace predicate */
} /* namespace shared */
} /* namespace spatial */
} /* namespace tbox */


#endif /* PREDICATE_SPATIAL_HPP_ */
