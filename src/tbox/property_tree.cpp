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
PropertyTree::keyExists(const std::string& key) const{
	return (node_.count(key) > 0);
}

std::vector<std::string>
PropertyTree::getKeys() const{
	std::vector<std::string> result;
	result.reserve(node_.size());
	for(auto kv : node_){
		result.emplace_back(kv.first);
	}
	return result;
}

bool PropertyTree::isPropertyTree(const std::string& key) const{
	return ( this->keyExists(key)     &&
			 this->key_is_list(key)   &&
			 this->list_has_keys(key) );
}

PropertyTree
PropertyTree::getPropertyTree(const std::string& key) const{
	if( !isPropertyTree(key) ){
		EH_ERROR("PropertyTree::getPropertyTree() key = " << key);
	}
	return this->get_ptree_impl(key);
}

PropertyTree
PropertyTree::getPropertyTree(const std::string& key, const PropertyTree& dval) const{
	if( !isPropertyTree(key) ){
		return dval;
	}
	return this->get_ptree_impl(key);
}


void
PropertyTree::setPropertyTree(const std::string& key, const PropertyTree& val){
	this->node_.put_child(key, val.node_);
}



// Key leads to a terminal value (not array or ptree)
bool PropertyTree::key_is_terminal(const std::string& key) const {
	ASSERT(this->keyExists(key));
	boost::property_tree::ptree const child = node_.get_child(key);
	return (child.empty() && !child.data().empty());
}

// Key leads to a list (array or ptree)
bool PropertyTree::key_is_list(const std::string& key) const {
	ASSERT(this->keyExists(key));
	boost::property_tree::ptree const child = node_.get_child(key);
	return  (!child.empty() && child.data().empty());
}

// Key leads to an array
bool PropertyTree::key_is_array(const std::string& key) const{
	return ( this->keyExists(key)     &&
			 this->key_is_list(key)   &&
			!this->list_has_keys(key) );
}

// Key leads to a list with key names
bool PropertyTree::list_has_keys(const std::string& key) const {
	ASSERT(this->key_is_list(key));
	for(auto const& kv_pair : node_.get_child(key)){
		if( !kv_pair.first.empty() ){
			return true;
		}
	}
	return false;
}


PropertyTree
PropertyTree::get_ptree_impl(const std::string& key) const{
	ASSERT(this->key_is_list(key) && !this->key_is_array(key));
	PropertyTree pt;
	pt.node_ = node_.get_child(key);
	return pt;
}




} // namespace tbox
} // namespace geoflow
