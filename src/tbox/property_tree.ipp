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
void PropertyTree::serialize(Archive & ar, const unsigned int version){
	//boost::property_tree::serialize(ar,this->node_,version);
	ar & this->node_;
}

template<typename T>
bool PropertyTree::isValue(const std::string& key) const{
	bool result = (keyExists(key) && key_is_terminal(key));
	if( result ){
		try{
			const T b = this->get_value_impl<T>(key);
		}
		catch(...){
			result = false;
		}
	}
	return result;
}

template<typename T>
bool PropertyTree::isArray(const std::string& key) const{
	bool result = (keyExists(key) && key_is_array(key));
	if( result ){
		try{
			std::vector<T> tmp = this->get_array_impl<T>(key);
		}
		catch(...){
			result = false;
		}
	}
	return result;
}



template<typename T>
T
PropertyTree::getValue(const std::string& key) const{
	if( !isValue<T>(key) ){
		EH_ERROR("PropertyTree::getValue() key = " << key);
	}
	return this->get_value_impl<T>(key);
}

template<typename T>
T
PropertyTree::getValue(const std::string& key, const T& dval) const{
	if( !isValue<T>(key) ){
		return dval;
	}
	return this->get_value_impl<T>(key);
}

template<typename T>
std::vector<T>
PropertyTree::getArray(const std::string& key) const{
	if( !isArray<T>(key) ){
		EH_ERROR("PropertyTree::getArray() key = " << key);
	}
	return this->get_array_impl<T>(key);
}

template<typename T>
std::vector<T>
PropertyTree::getArray(const std::string& key, const std::vector<T>& dval) const{
	if( !isArray<T>(key) ){
		return dval;
	}
	return this->get_array_impl<T>(key);
}


template<typename T>
void
PropertyTree::setValue(const std::string& key, const T& val){
	node_.put(key,val);
}

template<typename T>
void
PropertyTree::setArray(const std::string& key, const std::vector<T>& vec){
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
T
PropertyTree::get_value_impl(const std::string& key) const{
	ASSERT(key_is_terminal(key));
	return node_.get<T>(key);
}

template<typename T>
std::vector<T>
PropertyTree::get_array_impl(const std::string& key) const{
	ASSERT(key_is_array(key));
	std::vector<T> result;
	result.reserve(node_.get_child(key).size());
	for(auto const& kv_pair : node_.get_child(key)){
		result.emplace_back(kv_pair.second.get<T>(""));
	}
	return result;
}






} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_TBOX_PROPERTY_TREE_IPP_ */
