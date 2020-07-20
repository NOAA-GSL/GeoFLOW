/*
 * property_tree.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "tbox/property_tree.hpp"

#include "boost/property_tree/json_parser.hpp"

namespace geoflow {
namespace tbox {


PropertyTree::PropertyTree(){
}


PropertyTree::PropertyTree(const PropertyTree &pt) :
		node_(pt.node_){
}

PropertyTree::~PropertyTree(){
}


PropertyTree& PropertyTree::operator=(const PropertyTree &pt){
	node_ = pt.node_;
	return *this;
}

void PropertyTree::load_file(const std::string& file_name){
	try{
		boost::property_tree::read_json(file_name,this->node_);
	} catch(std::exception &e){
		EH_ERROR("PropertyTree::load_file(...) -> " <<  e.what());
	}
}


void PropertyTree::save_file(const std::string& file_name) const{
	try{
		boost::property_tree::write_json(file_name, node_);
	} catch(std::exception &e){
		EH_ERROR("PropertyTree::save_file(...) -> " <<  e.what());
	}
}

void PropertyTree::load_string(const std::string& content){
	try{
		std::stringstream ss;
		ss << content;
		boost::property_tree::read_json(ss,this->node_);
	} catch(std::exception &e){
		EH_ERROR("PropertyTree::load_string(...) -> " <<  e.what());
	}
}

bool
PropertyTree::keyExists(const std::string& key) const {
	return this->key_exists_(key);
}

std::vector<std::string>
PropertyTree::getKeys() const {
	std::vector<std::string> result;
	for(auto kv : node_){
		result.push_back(kv.first);
	}
	return result;
}

bool
PropertyTree::isPropertyTree(const std::string& key) const {
	return (this->key_exists_(key) &&
			this->key_has_children_(key) &&
			this->key_children_have_names_(key));
}

PropertyTree
PropertyTree::getPropertyTree(const std::string& key) const {
	if( not key_exists_(key) ){
		EH_ERROR("PropertyTree key '"<<key<<"' does not exist");
	}
	if( not key_has_children_(key) ){
		EH_ERROR("PropertyTree key '"<<key<<"' contains single value");
	}
	if( not key_children_have_names_(key) ){
		EH_ERROR("PropertyTree key '"<<key<<"' contains no child names");
	}
	return this->get_tree_impl_(key);
}

/**
 * Get the PropertyTree at key returned as type T
 *
 * Get the PropertyTree stored at key.  If the key
 * doesn't exist the provided default value is
 * returned.
 */
PropertyTree
PropertyTree::getPropertyTree(const std::string& key, const PropertyTree& dval) const {
	if( this->isPropertyTree(key) ){
		return dval;
	}
	return this->get_tree_impl_(key);
}

void
PropertyTree::setPropertyTree(const std::string& key, const PropertyTree& val){
	this->node_.put_child(key, val.node_);
}

bool
PropertyTree::key_exists_(const std::string& key) const {
	return (node_.count(key) > 0);
}

bool
PropertyTree::key_has_children_(const std::string& key) const {
	ASSERT( key_exists_(key) );
	return (node_.get_child(key).size() > 0);
}

bool
PropertyTree::key_children_have_names_(const std::string& key) const {
	ASSERT( key_exists_(key) );
	ASSERT( key_has_children_(key) );
	auto key_node =  node_.get_child(key);
	for(auto const& child_pair : key_node){ // Loop over child nodes
		auto child_key = child_pair.first;  // Get std::string name of child
		if( not child_key.empty() ){        // Test if std::string is empty
			return true;
		}
	}
	return false;
}


PropertyTree
PropertyTree::get_tree_impl_(const std::string& key) const {
	ASSERT(key_exists_(key));
	ASSERT(key_has_children_(key));
	ASSERT(key_children_have_names_(key));
	PropertyTree pt;
	pt.node_ = node_.get_child(key);
	return pt;
}

} // namespace tbox
} // namespace geoflow
