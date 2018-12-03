/*
 * property_tree.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_PROPERTY_TREE_HPP_
#define SRC_GEOFLOW_TBOX_PROPERTY_TREE_HPP_


#include <string>
#include <vector>

/**
 * Define the library used under the hood
 */
#include "boost/property_tree/ptree.hpp"

namespace geoflow {
namespace tbox {


class PropertyTree;
PropertyTree parse_json(const std::string& file_name);


/**
 * @brief
 * Class to read and write information to a file
 *
 *
 *
 */
class PropertyTree {

public:

	// *********************************************************************
	//                       Constructors / Destructors
	// *********************************************************************

	/**
	 * Default constructor for a PropertyTree
	 */
	PropertyTree();

	/**
	 * Copy constructor for a PropertyTree
	 */
	PropertyTree(const PropertyTree &DB);

	/**
	 * Default destructor for a PropertyTree
	 */
	virtual ~PropertyTree();

	// *********************************************************************
	//                                Copy
	// *********************************************************************

	/**
	 * Copy a PropertyTree
	 */
	PropertyTree& operator=(const PropertyTree &DB);


	// *********************************************************************
	//                                File I/O
	// *********************************************************************


	void load_file(const std::string& file_name);


	void save_file(const std::string& file_name) const;


	void load_string(const std::string& content);

	// *********************************************************************
	//                           Serialization
	// *********************************************************************

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	// *********************************************************************
	//                          Inquiry Functions
	// *********************************************************************

	/**
	 * Return true if the specified key exists in the database and false
	 * otherwise.
	 *
	 * @param key Key name to lookup.
	 */
	bool keyExists(const std::string& key) const;

	/**
	 * Return all keys in the database.
	 *
	 * Get all the keys at this level of the database.
	 */
	std::vector<std::string> getKeys() const;

	// *********************************************************************
	//                        Check Type Functions
	// *********************************************************************

	/**
	 * Check if value of type T exists at key
	 *
	 * Checks the key value to see if it exists and can
	 * be converted to type T.
	 */
	template<typename T>
	bool isValue(const std::string& key) const;

	/**
	 * Check if array of type T exists at key
	 *
	 * Checks the key value to see if it exists and can
	 * be converted to type array of T.
	 */
	template<typename T>
	bool isArray(const std::string& key) const;

	/**
	 * Check if a PropertyTree exists at key
	 *
	 * Checks the key value to see if it exists and is of
	 * type PropertyTree
	 */
	bool isPropertyTree(const std::string& key) const;

	// *********************************************************************
	//                          Get Functions
	// *********************************************************************

	/**
	 * Get the value at key returned as type T
	 *
	 * Get the value stored at key converted to type T.
	 * If the key doesn't exist or cannot be converted an error
	 * occurs.
	 */
	template<typename T>
	T getValue(const std::string& key) const;

	/**
	 * Return the array at key returned as type T
	 *
	 * Get the array stored at key converted to type std::vector<T>.
	 * If the key doesn't exist or cannot be converted an error
	 * occurs.
	 */
	template<typename T>
	std::vector<T> getArray(const std::string& key) const;


	/**
	 * Return the PropertyTree at key returned as type T
	 *
	 * Get the PropertyTree stored at key.  If the key
	 * doesn't exist or cannot be converted an error
	 * occurs.
	 */
	PropertyTree getPropertyTree(const std::string& key) const;



	// *********************************************************************
	//                          Get Functions
	// *********************************************************************

	/**
	 * Get the value at key returned as type T
	 *
	 * Get the value stored at key converted to type T.
	 * If the key doesn't exist the provided default value is
	 * returned.
	 */
	template<typename T>
	T getValue(const std::string& key, const T& dval) const;

	/**
	 * Get the array at key returned as type T
	 *
	 * Get the array stored at key converted to type std::vector<T>.
	 * If the key doesn't exist the provided default value is
	 * returned.
	 */
	template<typename T>
	std::vector<T> getArray(const std::string& key, const std::vector<T>& dval) const;

	/**
	 * Get the PropertyTree at key returned as type T
	 *
	 * Get the PropertyTree stored at key.  If the key
	 * doesn't exist the provided default value is
	 * returned.
	 */
	PropertyTree getPropertyTree(const std::string& key, const PropertyTree& dval) const;

	// *********************************************************************
	//                          Set Functions
	// *********************************************************************

	/**
	 * Place value of type T into PropertyTree.
	 *
	 * Insert a value of type T into the PropertyTree
	 * using the provided key name.  If this key is already
	 * associated with another value it will be overwritten.
	 */
	template<typename T>
	void setValue(const std::string& key, const T& val);

	/**
	 * Place array of type T into PropertyTree.
	 *
	 * Insert an array of type T into the PropertyTree
	 * using the provided key name.  If this key is already
	 * associated with another value it will be overwritten.
	 */
	template<typename T>
	void setArray(const std::string& key, const std::vector<T>& val);

	/**
	 * Place PropertyTree into PropertyTree.
	 *
	 * Insert a PropertyTree T into the PropertyTree
	 * using the provided key name.  If this key is already
	 * associated with another value it will be overwritten.
	 */
	void setPropertyTree(const std::string& key, const PropertyTree& val);


private:
	friend class boost::serialization::access;
	boost::property_tree::ptree node_;

	// *********************************************************************
	//              [private] Helper Functions
	// *********************************************************************

	// Determine if key is not array or ptree (No error checks)
	bool key_is_terminal(const std::string& key) const;

	// Determine if key is array or ptree (No error checks)
	bool key_is_list(const std::string& key) const;

	// Determine if key is an array (No error checks)
	bool key_is_array(const std::string& key) const;

	// Determine if list has none empty keys (No error checks)
	bool list_has_keys(const std::string& key) const;

	// Implementation of getValue (No error checks)
	template<typename T>
	T get_value_impl(const std::string& key) const;

	// Implementation of getArray (No error checks)
	template<typename T>
	std::vector<T> get_array_impl(const std::string& key) const;

	// Implementation of getPropertyTree (No error checks)
	PropertyTree get_ptree_impl(const std::string& key) const;
};


} // namespace tbox
} // namespace geoflow

#include "tbox/property_tree.ipp"

#endif /* SRC_GEOFLOW_TBOX_PROPERTY_TREE_HPP_ */
