/*
 * property_tree.ipp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_PROPERTY_TREE_IPP_
#define SRC_GEOFLOW_TBOX_PROPERTY_TREE_IPP_

#include "tbox/assert.hpp"
#include "tbox/error_handler.hpp"

#include "boost/archive/basic_archive.hpp"
#include "boost/property_tree/ptree_serialization.hpp"

namespace geoflow {
namespace tbox {


template<class Archive>
void
PropertyTree::serialize(Archive & ar, const unsigned int version){
	ar & this->node_;
}


template<typename T>
bool
PropertyTree::isValue(const std::string& key) const {
	return (this->key_exists_(key)            &&
		   (not this->key_has_children_(key)) &&
		   (this->value_is_type_<T>(key)) );
}

template<typename T>
bool
PropertyTree::isArray(const std::string& key) const {
	return (this->key_exists_(key) &&
			this->key_has_children_(key) &&
		    (not this->key_children_have_names_(key)) &&
			this->array_is_type_<T>(key)  );
}

template<typename T>
bool
PropertyTree::isArray2D(const std::string& key) const {
	return (this->key_exists_(key)                    &&
			this->key_has_children_(key)              &&
		    (not this->key_children_have_names_(key)) &&
			this->array_is_array_type_<T>(key) );
}

template<typename T>
T
PropertyTree::getValue(const std::string& key) const {
	if( not key_exists_(key) ){
		EH_ERROR("Key '"<<key<<"' does not exist");
	}
	if( key_has_children_(key) ){
		EH_ERROR("Expected single value in key '" << key << "'");
	}
	if( not value_is_type_<T>(key) ) {
		EH_ERROR("Value within key '"<<key<<"' could not be converted to requested type");
	}
	return node_.get<T>(key);
}

template<typename T>
T
PropertyTree::getValue(const std::string& key, const T& dval) const {
	if( not this->isValue<T>(key) ) {
		return dval;
	}
	return node_.get<T>(key);
}

template<typename T>
std::vector<T>
PropertyTree::getArray(const std::string& key) const {
	if( not key_exists_(key) ){
		EH_ERROR("Key '"<<key<<"' does not exist");
	}
	if( not key_has_children_(key) ){
		EH_ERROR("Expected array in key '" << key << "'");
	}
	if( key_children_have_names_(key) ){
		EH_ERROR("Expected array in key '" << key << "'");
	}
	if( not array_is_type_<T>(key) ) {
		EH_ERROR("Value within key '"<<key<<"' could not be converted to requested array type");
	}
	return this->get_array_impl_<T>(key);
}

template<typename T>
std::vector<T>
PropertyTree::getArray(const std::string& key, const std::vector<T>& dval) const {
	if( not this->isArray<T>(key) ) {
		return dval;
	}
	return this->get_array_impl_<T>(key);
}

template<typename T>
std::vector<std::vector<T>>
PropertyTree::getArray2D(const std::string& key) const {
	if( not key_exists_(key) ){
		EH_ERROR("Key '"<<key<<"' does not exist");
	}
	if( not key_has_children_(key) ){
		EH_ERROR("Expected 2D array in key '" << key << "'");
	}
	if( key_children_have_names_(key) ){
		EH_ERROR("Expected 2D array in key '" << key << "'");
	}
	if( not array_is_array_type_<T>(key) ) {
		EH_ERROR("Value within key '"<<key<<"' could not be converted to requested 2D array type");
	}
	return this->get_array_of_array_impl_<T>(key);
}

template<typename T>
std::vector<std::vector<T>>
PropertyTree::getArray2D(const std::string& key,
		                 const std::vector<std::vector<T>>& dval) const {
	if( not this->isArray2D<T>(key) ) {
		return dval;
	}
	return this->get_array_of_array_impl_<T>(key);
}

template<typename T>
void
PropertyTree::setValue(const std::string& key, const T& val) {
	node_.put(key,val);
}


template<typename T>
void
PropertyTree::setArray(const std::string& key, const std::vector<T>& vec) {
	  boost::property_tree::ptree array;
	  const std::string empty_string;
	  for(auto const& val : vec){
	    boost::property_tree::ptree array_entry;
	    array_entry.put(empty_string,val);
	    array.push_back(std::make_pair(empty_string,array_entry));
	  }
	  this->node_.put_child(key,array);
}


template<typename T>
bool
PropertyTree::value_is_type_(const std::string& key) const {
	ASSERT( key_exists_(key) );
	ASSERT( not key_has_children_(key) );
	return node_.get_optional<T>(key).is_initialized(); // Boost 1.65 only has this call
	//return node_.get_optional<T>(key).has_value();    // Boost 1.68 added this call
}

template<typename T>
bool
PropertyTree::array_is_type_(const std::string& key) const {
	ASSERT( key_exists_(key) );
	ASSERT( key_has_children_(key) );
	ASSERT(not key_children_have_names_(key));
	auto key_node = node_.get_child(key);
	for(auto const& child_pair : key_node){ // Loop over child nodes
		if( not child_pair.second.get_value_optional<T>() ) {
			return false;
		}
	}
	return true;
}

template<typename T>
bool
PropertyTree::array_is_array_type_(const std::string& key) const {
	ASSERT( key_exists_(key) );
	ASSERT( key_has_children_(key) );
	ASSERT(not key_children_have_names_(key));
	auto key_node = node_.get_child(key);
	// Loop over nodes of child with key
	for(auto const& child_pair : node_.get_child(key)){
	    // Child must have children
		if( child_pair.second.empty() ){
			return false;
		}
		// Loop over child nodes
		for(auto const& sub_child_pair : child_pair.second){
			// Test if sub_child key is empty std::string
			if( not sub_child_pair.first.empty() ){
				return false;
			}
			// Test if sub_child value can be converted to T
			if( not sub_child_pair.second.get_value_optional<T>() ) {
				return false;
			}
		}
	}
	return true;
}

// Implementation of getArray (No error checks)
template<typename T>
std::vector<T>
PropertyTree::get_array_impl_(const std::string& key) const {
	ASSERT(key_exists_(key));
	ASSERT(key_has_children_(key));
	ASSERT(not key_children_have_names_(key));
	ASSERT(array_is_type_<T>(key));
	std::vector<T> result;
	auto key_node = node_.get_child(key);
	for(auto const& child_pair : key_node){
		result.push_back(child_pair.second.get_value<T>());
	}
	return result;
}

// Implementation of getArray (No error checks)
template<typename T>
std::vector<std::vector<T>>
PropertyTree::get_array_of_array_impl_(const std::string& key) const {
	ASSERT(key_exists_(key));
	ASSERT(key_has_children_(key));
	ASSERT(not key_children_have_names_(key));
	ASSERT(array_is_array_type_<T>(key));

	std::vector<std::vector<T>> result;
	auto key_node = node_.get_child(key);
	for(auto const& child_pair : key_node){                  // Loop over child nodes
		std::vector<T> row;
		for(auto const& sub_child_pair : child_pair.second){ // Loop over sub-child nodes
			row.push_back(sub_child_pair.second.get_value<T>());
		}
		result.push_back(row);
	}
	return result;
}




} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_TBOX_PROPERTY_TREE_IPP_ */
