/*
 * box.hpp
 *
 *  Created on: Sep 15, 2020
 *      Author: bflynt
 */

#ifndef BOUNDING_BOX_HPP_
#define BOUNDING_BOX_HPP_


#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>

namespace tbox {
namespace spatial {
namespace bound {

//
// Declare the following function
//
template<typename RANGE_TYPE, std::size_t NDIM>
class Box;

template<typename T, std::size_t N>
bool Disjoint(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
bool Intersects(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
bool Overlaps(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
bool Contains(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
bool Covers(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
bool Equals(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
typename Box<T,N>::value_type
Nearest(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
typename Box<T,N>::value_type
Centroid(Box<T,N> const& a, Box<T,N> const& b);

template<typename T, std::size_t N>
typename Box<T,N>::value_type
Furthest(Box<T,N> const& a, Box<T,N> const& b);









template<typename RANGE_TYPE, std::size_t NDIM>
class Box final {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
public:
	using self_type  = Box<RANGE_TYPE,NDIM>;
	using array_type = std::array<RANGE_TYPE,NDIM>;
	using value_type = typename array_type::value_type;
	using size_type  = typename array_type::size_type;
	static constexpr size_type ndim = NDIM;


	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------
public:

	Box() = default;

	Box(const Box& box) = default;

	Box(Box&& box) = default;

	Box(const array_type& min_corner,
		const array_type& max_corner)
	: min_(min_corner),
	  max_(max_corner){
	}

	~Box() = default;

	//-------------------------------------------------------------------------
	// Access Operators
	//-------------------------------------------------------------------------

	value_type min(const size_type dim) const noexcept {
			return min_[dim];
	}

	value_type max(const size_type dim) const noexcept {
			return max_[dim];
	}

	value_type center(const size_type dim) const noexcept {
		return 0.5*(max_[dim] + min_[dim]);
	}

	value_type length(const size_type dim) const noexcept {
		return (max_[dim] - min_[dim]);
	}
//
//	array_type const& min() const noexcept {
//			return min_;
//	}
//
//	array_type const& max() const noexcept {
//			return max_;
//	}

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	Box& operator=(const Box& box) = default;

	Box& operator=(Box&& box) = default;

	//-------------------------------------------------------------------------
	// Properties
	//-------------------------------------------------------------------------

	value_type area() const noexcept {
		value_type s = max_[0] - min_[0];
		for(size_type i = 1; i < min_.size(); ++i){
			s *= max_[i] - min_[i];
		}
		return s;
	}

	//-------------------------------------------------------------------------
	// Operations
	//-------------------------------------------------------------------------

	/**
	 * Set this BBox
	 *
	 * Set the bounding box to provided values
	 */
	void set(const array_type& min_corner,
	         const array_type& max_corner) noexcept {
		min_ = min_corner;
		max_ = max_corner;
	}

	/**
	 * Reset this BBox to default
	 *
	 * Reset the bounding box to the most negative volume
	 * where a stretch by any valid box will result in an
	 * size increase.
	 */
	void reset() noexcept {
		for(size_type i = 0; i < min_.size(); ++i){
			min_[i] = std::numeric_limits<value_type>::max();
			max_[i] = std::numeric_limits<value_type>::lowest();
		}
	}

	/**
	 * Stretch this BBox to fit provided box
	 *
	 * Stretch this bounding box to enclose the provided box
	 */
	void stretch(const Box& box) noexcept {
		using std::min;
		using std::max;
		for(size_type i = 0; i < min_.size(); ++i){
			min_[i] = min(min_[i],box.min_[i]);
			max_[i] = max(max_[i],box.max_[i]);
		}
	}

	//-------------------------------------------------------------------------
	// Boolean Operators
	//-------------------------------------------------------------------------

	/**
	 * Test if two Boxes are equal
	 *
	 * Both boxes are equal if the coordinates for each dimensions are equal
	 */
	bool
	operator==(self_type const& other) const {
		return ((this->min_ == other.min_) and (this->max_ == other.max_));
	}

	/**
	 * Test if two Boxes are not equal
	 *
	 * Both boxes are equal if the coordinates for each dimensions are equal
	 */
	bool
	operator!=(self_type const& other) const {
		return (not (*this == other));
	}

	//-------------------------------------------------------------------------
	// Tests [Friends]
	//-------------------------------------------------------------------------
	friend bool Disjoint<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend bool Intersects<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend bool Overlaps<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend bool Contains<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend bool Covers<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend bool Equals<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);

	friend value_type Nearest<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend value_type Centroid<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);
	friend value_type Furthest<RANGE_TYPE,NDIM>(self_type const& a, self_type const& b);

	//-------------------------------------------------------------------------
	// Data [Private]
	//-------------------------------------------------------------------------
protected:
	array_type min_;
	array_type max_;



	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & min_;
        ar & max_;
    }
};


/**
 * Test if Box A and Box B are disjoint
 *
 * Test if the Box A and Box B do not touch at any
 * location.
 */
template<typename T, std::size_t N>
bool Disjoint(Box<T,N> const& a, Box<T,N> const& b) {
	// return (not Intersects(a,b);
	for(std::size_t i = 0; i < N; ++i){
		if((a.max_[i] < b.min_[i]) or (b.max_[i] < a.min_[i])) {
			return true;
		}
	}
	return false;
}

/**
 * Test if Box A intersects Box B
 *
 * Test if the intersection of Box A with Box B
 * would result in a Box with an area equal to or greater
 * than zero. (ie. True if they touch anywhere)
 */
template<typename T, std::size_t N>
bool Intersects(Box<T,N> const& a, Box<T,N> const& b) {
	//return ((a.min_ <= b.max_) and (a.max_ >= b.min_));
	return std::equal(a.min_.cbegin(), a.min_.cend(), b.max_.cbegin(), std::less_equal<T>()) and
		   std::equal(a.max_.cbegin(), a.max_.cend(), b.min_.cbegin(), std::greater_equal<T>());
}

/**
 * Test if Box A overlaps Box B
 *
 * Test if the intersection of this Box with Box B
 * would result in a Box with an area greater than zero.
 * Also, can be thought of as an intersection where the
 * boxes more than just touch they "overlap" some amount.
 */
template<typename T, std::size_t N>
bool Overlaps(Box<T,N> const& a, Box<T,N> const& b) {
	//return ((a.min_ < b.max_) and (a.max_ > b.min_));
	return std::equal(a.min_.cbegin(), a.min_.cend(), b.max_.cbegin(), std::less<T>()) and
		   std::equal(a.max_.cbegin(), a.max_.cend(), b.min_.cbegin(), std::greater<T>());
}

/**
 * Test if Box A fully Contains Box B
 *
 * Test if Box A extents are further or equal to
 * Box B extents.  They can be touching.
 */
template<typename T, std::size_t N>
bool Contains(Box<T,N> const& a, Box<T,N> const& b) {
	//return ((a.min_ <= b.min_) and (a.max_ >= b.max_));
	return std::equal(a.min_.cbegin(), a.min_.cend(), b.min_.cbegin(), std::less_equal<T>()) and
		   std::equal(a.max_.cbegin(), a.max_.cend(), b.max_.cbegin(), std::greater_equal<T>());
}

/**
 * Test if Box A fully Covers Box B
 *
 * Test if Box A extents are further than Box B in every
 * direction.  They can not be touching on any side.
 */
template<typename T, std::size_t N>
bool Covers(Box<T,N> const& a, Box<T,N> const& b) {
	//return ((a.min_ < b.min_) and (a.max_ > b.max_));
	return std::equal(a.min_.cbegin(), a.min_.cend(), b.min_.cbegin(), std::less<T>()) and
		   std::equal(a.max_.cbegin(), a.max_.cend(), b.max_.cbegin(), std::greater<T>());
}

/**
 * Test if Box A and Box A are equal
 *
 * Test if the coordinates for each dimensions are equal
 */
template<typename T, std::size_t N>
bool Equals(Box<T,N> const& a, Box<T,N> const& b) {
	//return ((a.min_ == b.min_) and (a.max_ == b.max_));
	return std::equal(a.min_.cbegin(), a.min_.cend(), b.min_.cbegin()) and
		   std::equal(a.max_.cbegin(), a.max_.cend(), b.max_.cbegin());
}


/**
 * Nearest distance metric between Box A and Box B
 *
 * Calculate a measure of the nearest distance between
 * Box A and Box B. The distance may not be the true
 * distance but rather a distance metric which can be
 * compared against another measure.  In the case of
 * the Box it is the Euclidean Distance Squared.
 *
 * Notice:
 * If the boxes overlap in any coordinate direction that
 * contribution is zero. Therefore, boxes which touch will
 * have a nearest distance of Zero.
 */
template<typename T, std::size_t N>
typename Box<T,N>::value_type
Nearest(Box<T,N> const& a, Box<T,N> const& b) {
	using std::max;
	using std::pow;
	using size_type  = typename Box<T,N>::size_type;
	using value_type = typename Box<T,N>::value_type;
	constexpr value_type zero = 0;
	value_type dist_sq = 0;
	for(size_type i = 0; i < N; ++i){
		auto b_bigger  = max(zero, b.min_[i]-a.max_[i]);
		auto b_smaller = max(zero, a.min_[i]-b.max_[i]);
		auto b_dist    = max(b_bigger,b_smaller);
		dist_sq += pow(b_dist,2);
	}
	return dist_sq;
}

/**
 * Center distance metric between Box A and Box B
 *
 * Calculate a measure of the center distance between
 * Box A and Box B. The distance may not be the true
 * distance but rather a distance metric which can be
 * compared against another measure.  In the case of
 * the Box it is the Euclidean Distance Squared.
 */
template<typename T, std::size_t N>
typename Box<T,N>::value_type
Centroid(Box<T,N> const& a, Box<T,N> const& b) {
	using std::pow;
	using size_type  = typename Box<T,N>::size_type;
	using value_type = typename Box<T,N>::value_type;
	value_type dist_sq = 0;
	for(size_type i = 0; i < N; ++i){
		//auto ac = 0.5 * (a.max_[i]+a.min_[i]);
		//auto bc = 0.5 * (b.max_[i]+b.min_[i]);
		//auto dc = std::pow(ac-bc,2);
		dist_sq += pow(0.5*(a.max_[i]+a.min_[i]-b.max_[i]-b.min_[i]),2);
	}
	return dist_sq;
}

/**
 * Furthest distance metric between Box A and Box B
 *
 * Calculate a measure of the furthest distance between
 * any two locations on Box A and Box B. The distance
 * may not be the true distance but rather a distance
 * metric which can be compared against another measure.
 * In the case of the Box it is the Euclidean Distance
 * Squared.
 *
 * Notice:
 * If the boxes are inside of each other in any coordinate
 * direction that contribution is zero. Therefore, boxes
 * which contain another box have a furthest distance of Zero.
 */
template<typename T, std::size_t N>
typename Box<T,N>::value_type
Furthest(Box<T,N> const& a, Box<T,N> const& b) {
	using std::max;
	using std::pow;
	using size_type  = typename Box<T,N>::size_type;
	using value_type = typename Box<T,N>::value_type;
	value_type dist_sq = 0;
	for(size_type i = 0; i < N; ++i){
		if( (a.max_[i] < b.max_[i]) xor (b.min_[i] < a.min_[i]) ){
			auto b_bigger  = pow(b.max_[i]-a.min_[i],2);
			auto b_smaller = pow(b.min_[i]-a.max_[i],2);
			dist_sq += max(b_bigger,b_smaller);
		}
	}
	return dist_sq;
}


} /* namespace bound */
} /* namespace spatial */
} /* namespace tbox */


//-------------------------------------------------------------------------
// Stream Operators
//-------------------------------------------------------------------------
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const tbox::spatial::bound::Box<T,N>& a){
	os << "min(";
	for(auto i = 0; i < N; ++i){
		os << " " << a.min(i);
	}
	os << ") max(";
	for(auto i = 0; i < N; ++i){
		os << " " << a.max(i);
	}
	os << ")";
	return os;
}


#endif /* BOUNDING_BOX_HPP_ */
