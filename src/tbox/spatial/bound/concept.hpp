/*
 * concept.hpp
 *
 *  Created on: Sep 16, 2020
 *      Author: bflynt
 */

#ifndef BOUNDING_CONCEPT_HPP_
#define BOUNDING_CONCEPT_HPP_

//#include <concepts>
#include <type_traits>

namespace tbox {
namespace spatial {
namespace bound {

/**
 * Bound Concept
 *
 * Concept defines all functions a bounding region
 * must support to be used by the index algorithms.
 */
/**
template<typename T>
concept BoundConcept = requires(T a, T b) {
	typename T::value_type;
	typename T::size_type;
	{T::ndim}                            -> std::same_as<typename T::size_type>;
	std::is_default_constructible<T>();
	std::is_copy_constructible<T>();
	std::is_copy_assignable<T>();
	{a.min((typename T::size_type)0)}    -> std::same_as<typename T::value_type>;
	{a.max((typename T::size_type)0)}    -> std::same_as<typename T::value_type>;
	{a.length((typename T::size_type)0)} -> std::same_as<typename T::value_type>;
	{a.area()}                           -> std::same_as<typename T::value_type>;
	{a.reset()}                          -> std::same_as<void>;
	{a.stretch(b)}                       -> std::same_as<void>;
	{Disjoint(a,b)}                      -> std::same_as<bool>;
	{Intersects(a,b)}                    -> std::same_as<bool>;
	{Overlaps(a,b)}                      -> std::same_as<bool>;
	{Contains(a,b)}                      -> std::same_as<bool>;
	{Covers(a,b)}                        -> std::same_as<bool>;
	{Equals(a,b)}                        -> std::same_as<bool>;
	{Nearest(a,b)}                       -> std::same_as<typename T::value_type>;
	{Centroid(a,b)}                      -> std::same_as<typename T::value_type>;
	{Furthest(a,b)}                      -> std::same_as<typename T::value_type>;
};
*/

/**
 * Union of two boxes
 * Returns a box large enough to fit both boxes
 */
template<typename T>
T
Union(const T& a, const T& b){
	T ans(a);
	ans.stretch(b);
	return ans;
}

/**
 * Get increased area of Box A needed to contain Box B
 *
 * Returns the increase in area that Box A would be required to
 * undergo so that Box B can fit inside.
 */
template<typename T>
typename T::value_type
IncreaseToHold(const T& a, const T& b){
	return (Union(a,b).area() - a.area());
}


} /* namespace bound */
} /* namespace spatial */
} /* namespace tbox */

#endif /* BOUNDING_CONCEPT_HPP_ */
