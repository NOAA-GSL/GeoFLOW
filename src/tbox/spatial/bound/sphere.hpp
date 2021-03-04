/*
 * sphere.hpp
 *
 *  Created on: Sep 16, 2020
 *      Author: bflynt
 */

#ifndef BOUNDING_SPHERE_HPP_
#define BOUNDING_SPHERE_HPP_


#include <algorithm>
#include <array>
#include <cmath>          // std::pow(
#include <functional>
#include <limits>
#include <numeric>
#include <ostream>

namespace tbox {
namespace spatial {
namespace bound {

namespace detail {

/**
 * Constant to calculate the volume of a sphere in N dimensions
 *
 * Overloaded struct to hold constants for calculating the volume
 * of a sphere in NDIM dimensions.  The formula followed is
 * volume = sphere_multipler * radius^NDIM.
 *
 * @tparam NDIM Number of dimensions
 */
template<std::size_t NDIM>
struct sphere_multipler {
};

template<>
struct sphere_multipler<1> {
	static constexpr long double value = 1;
};
template<>
struct sphere_multipler<2> {
	static constexpr long double value = M_PI;
};
template<>
struct sphere_multipler<3> {
	static constexpr long double value = 4*M_PI/3;
};
template<>
struct sphere_multipler<4> {
	static constexpr long double value = std::pow(M_PI,2)/2;
};
template<>
struct sphere_multipler<5> {
	static constexpr long double value = 8*std::pow(M_PI,2)/15;
};
template<>
struct sphere_multipler<6> {
	static constexpr long double value = std::pow(M_PI,3)/6;
};
template<>
struct sphere_multipler<7> {
	static constexpr long double value = 16*std::pow(M_PI,3)/105;
};
template<>
struct sphere_multipler<8> {
	static constexpr long double value = std::pow(M_PI,4)/24;
};
template<>
struct sphere_multipler<9> {
	static constexpr long double value = 32*std::pow(M_PI,4)/945;
};
template<>
struct sphere_multipler<10> {
	static constexpr long double value = std::pow(M_PI,5)/120;
};
template<>
struct sphere_multipler<11> {
	static constexpr long double value = 64*std::pow(M_PI,5)/10395;
};
template<>
struct sphere_multipler<12> {
	static constexpr long double value = std::pow(M_PI,6)/720;
};


} /* namespace detail */

//
// Declare the following function
//
template<typename RANGE_TYPE, std::size_t NDIM>
class Sphere;

template<typename T, std::size_t N>
bool Disjoint(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
bool Intersects(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
bool Overlaps(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
bool Contains(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
bool Covers(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
bool Equals(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Nearest(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Centroid(Sphere<T,N> const& a, Sphere<T,N> const& b);

template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Furthest(Sphere<T,N> const& a, Sphere<T,N> const& b);




template<typename RANGE_TYPE, std::size_t NDIM>
class Sphere final {

	//-------------------------------------------------------------------------
	// Types & Constants
	//-------------------------------------------------------------------------
public:
	using self_type  = Sphere<RANGE_TYPE,NDIM>;
	using array_type = std::array<RANGE_TYPE,NDIM>;
	using value_type = typename array_type::value_type;
	using size_type  = typename array_type::size_type;
	static constexpr size_type ndim = NDIM;

	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------
public:

	Sphere() = default;

	Sphere(const Sphere& other) = default;

	Sphere(Sphere&& other) = default;

	Sphere(const array_type& center,
		   const value_type& radius)
	: center_(center),
	  radius_(radius){
	}

	~Sphere() = default;

	//-------------------------------------------------------------------------
	// Access Operators
	//-------------------------------------------------------------------------

	//-------------------------------------------------------------------------
	// Access Operators
	//-------------------------------------------------------------------------

	value_type min(const size_type dim) const noexcept {
			return (center_[dim] - radius_);
	}

	value_type max(const size_type dim) const noexcept {
			return (center_[dim] + radius_);
	}

	value_type length(const size_type dim) const noexcept {
			return (2 * radius_);
	}

	value_type center(const size_type dim) const noexcept {
		return center_[dim];
	}

	value_type radius() const noexcept {
		return radius_;
	}

	//-------------------------------------------------------------------------
	// Assignment Operators
	//-------------------------------------------------------------------------

	Sphere& operator=(const Sphere& other) = default;

	Sphere& operator=(Sphere&& other) = default;

	//-------------------------------------------------------------------------
	// Properties
	//-------------------------------------------------------------------------

	value_type area() const noexcept {
		return detail::sphere_multipler<NDIM>::value * std::pow(radius_,NDIM);
	}

	//-------------------------------------------------------------------------
	// Operations
	//-------------------------------------------------------------------------

	/**
	 * Set this BSphere
	 *
	 * Set the bounding Sphere to provided values
	 */
	void set(const array_type& center,
	         const value_type& radius) noexcept {
		center_ = center;
		radius_ = radius;
	}

	/**
	 * Reset this Sphere to default
	 *
	 * Reset the sphere to a negative volume
	 * where a stretch by any valid sphere will result in an
	 * size increase.
	 */
	void reset() noexcept {
		center_.fill(0);
		radius_ = std::numeric_limits<value_type>::lowest();
	}

	/**
	 * Stretch this Sphere to fit provided
	 *
	 * Reposition this bounding sphere to enclose the provided sphere
	 */
	void stretch(const Sphere& other) noexcept {

		// Calculate the distance between centers
		auto distance = std::sqrt(distance_squared_(other));

		// If this Sphere encloses the other Sphere do nothing
		if(distance+other.radius_ <= radius_){
			((void)0); // No Operation
		}

		// If this Sphere is enclosed by other then become other Sphere
		else if( distance+radius_ <= other.radius_ ) {
			radius_ = other.radius_;
			center_ = other.center_;
		}

		// We must adjust ourself to enclose both
		else {

			auto new_radius = 0.5*(radius_ + other.radius_ + distance);
			auto dterm = (new_radius - radius_)/distance;
			for(size_type i = 0; i < center_.size(); ++i){
				center_[i] += (other.center_[i] - center_[i]) * dterm;
			}
			radius_ = new_radius;
		}
	}

	/**
	 * Scale this Sphere radius
	 *
	 * Scale this Sphere radius by s
	 */
	void scale(const value_type& s) noexcept {
		radius_ *= s;
	}

	//-------------------------------------------------------------------------
	// Boolean Operators
	//-------------------------------------------------------------------------

	/**
	 * Test if two Spheres are equal
	 *
	 * Both Spheres are equal if radius and center is same
	 */
	bool
	operator==(self_type const& other) const {
		return ((this->radius_ == other.radius_) and
			    (this->center_ == other.center_));
	}

	/**
	 * Test if two Spheres are not equal
	 *
	 * Both Spheres are not equal if radius or center is different
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
	array_type center_;
	value_type radius_;

	/**
	 * Get the distance squared between this and another sphere
	 */
	value_type distance_squared_(const Sphere& other) const {
		value_type ans = 0;
		for(size_type i = 0; i < center_.size(); ++i){
			ans += std::pow(center_[i] - other.center_[i], 2);
		}
		return ans;
	}
};

/**
 * Test if Sphere A and Sphere B are disjoint
 *
 * Test if the Sphere A and Sphere B do not touch at any
 * location.
 */
template<typename T, std::size_t N>
bool Disjoint(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	return (not Intersects(a,b));
}

/**
 * Test if Sphere A intersects Sphere B
 *
 * Test if the intersection of Sphere A with Sphere B
 * would result in a Sphere with an area equal to or greater
 * than zero. (ie. True if they touch anywhere)
 */
template<typename T, std::size_t N>
bool Intersects(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	auto dist_sq = Centroid(a,b);
	auto rad_sq  = pow(a.radius_+b.radius_, 2);
	return (dist_sq <= rad_sq);
}

/**
 * Test if Sphere A overlaps Sphere B
 *
 * Test if the intersection of this Sphere with Sphere B
 * would result in a Sphere with an area greater than zero.
 * Also, can be thought of as an intersection where the
 * Spherees more than just touch they "overlap" some amount.
 */
template<typename T, std::size_t N>
bool Overlaps(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	auto dist_sq = Centroid(a,b);
	auto rad_sq  = pow(a.radius_+b.radius_, 2);
	return (dist_sq < rad_sq);
}

/**
 * Test if Sphere A fully Contains Sphere B
 *
 * Test if Sphere A extents are further or equal to
 * Sphere B extents.  They can be touching.
 */
template<typename T, std::size_t N>
bool Contains(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	if(b.radius_ > a.radius_)
		return false;
	auto dist_sq = Centroid(a,b);
	auto rad_sq  = pow(a.radius_-b.radius_, 2);
	return (dist_sq <= rad_sq);
}

/**
 * Test if Sphere A fully Covers Sphere B
 *
 * Test if Sphere A extents are further than Sphere B in every
 * direction.  They can not be touching on any side.
 */
template<typename T, std::size_t N>
bool Covers(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	if(b.radius_ >= a.radius_)
		return false;
	auto dist_sq = Centroid(a,b);
	auto rad_sq  = pow(a.radius_-b.radius_, 2);
	return (dist_sq < rad_sq);
}

/**
 * Test if Sphere A and Sphere A are equal
 *
 * Test if the coordinates for each dimensions are equal
 */
template<typename T, std::size_t N>
bool Equals(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	return (a.radius_ == b.radius_) and (a.center_ == b.center_);
}


/**
 * Nearest distance metric between Sphere A and Sphere B
 *
 * Calculate a measure of the nearest distance between
 * Sphere A and Sphere B. The distance may not be the true
 * distance but rather a distance metric which can be
 * compared against another measure.  In the case of
 * the Sphere it is the Euclidean Distance Squared.
 *
 * Notice:
 * Spheres which touch will have a nearest distance of Zero.
 */
template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Nearest(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	using std::max;
	using value_type = typename Sphere<T,N>::value_type;
	constexpr value_type zero = 0;
	auto dist_sq  = Centroid(a,b);
	auto rad_a_sq = pow(a.radius_,2);
	auto rad_b_sq = pow(b.radius_,2);
	return max(zero, dist_sq-rad_a_sq-rad_b_sq);
}

/**
 * Center distance metric between Sphere A and Sphere B
 *
 * Calculate a measure of the center distance between
 * Sphere A and Sphere B. The distance may not be the true
 * distance but rather a distance metric which can be
 * compared against another measure.  In the case of
 * the Sphere it is the Euclidean Distance Squared.
 */
template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Centroid(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	using size_type  = typename Sphere<T,N>::size_type;
	using value_type = typename Sphere<T,N>::value_type;
	value_type dist_sq = 0;
	for(size_type i = 0; i < N; ++i){
		dist_sq += pow(a.center_[i] - b.center_[i], 2);
	}
	return dist_sq;
}

/**
 * Furthest distance metric between Sphere A and Sphere B
 *
 * Calculate a measure of the furthest distance between
 * Sphere A and Sphere B. The distance may not be the true
 * distance but rather a distance metric which can be
 * compared against another measure.  In the case of
 * the Sphere it is the Euclidean Distance Squared.
 *
 * Notice:
 * Spheres inside each other will have a furthest distance of Zero.
 */
template<typename T, std::size_t N>
typename Sphere<T,N>::value_type
Furthest(Sphere<T,N> const& a, Sphere<T,N> const& b) {
	using std::pow;
	using std::sqrt;
	auto dist_sq  = Centroid(a,b);
	auto d_rad_sq = pow(a.radius_-b.radius_,2);

	// Sphere contains Sphere
	if( dist_sq <= d_rad_sq ){
		return 0;
	}

	// Spheres do/do not touch
	auto dist = sqrt(dist_sq);
	auto rads = a.radius_ + b.radius_;
	return pow(rads + dist,2);
}


} /* namespace bound */
} /* namespace spatial */
} /* namespace tbox */

//-------------------------------------------------------------------------
// Stream Operators
//-------------------------------------------------------------------------
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, tbox::spatial::bound::Sphere<T,N> const& a){
	os << "cen(";
	for(auto i = 0; i < N; ++i){
		os << " " << a.center(i);
	}
	os << ") rad(";
	os << " " << a.radius();
	os << ")";
	return os;
}



#endif /* BOUNDING_SPHERE_HPP_ */
