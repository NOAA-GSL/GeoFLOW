/*
 * command_line.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: bryan.flynt
 */
#ifndef SRC_XSTD_COMMAND_LINE_HPP_
#define SRC_XSTD_COMMAND_LINE_HPP_

#include <map>
#include <string>

namespace geoflow {
namespace tbox {

/// CommandLine parses simple input line commands
/**
 * The class parses the command line arguments into string tokens.
 * It assumes all keys start with a "-" symbol.
 *
 * Keys that start with a single "-" will be treated as multiple
 * key arguments with the following value (if provided) applying
 * to each character in the key.
 *
 * Keys that start with a double "--" will be treated as a single
 * key argument with the following value (if provided) applying
 * only to that full key word.
 *
 * Usage:
 * \code{.cpp}
 * int main(int argc, char* argv[]){
 *
 *    // Parse the command line
 *    CommandLine cmd_line(argc,argv);
 *
 *    // Get the "-i" or "--i" argument (with a default)
 *    std::string file_name = cmd_line.get("i", "default.txt");
 *
 *    // Get the "--verbose" argument (with a default)
 *    std::string verbose = cmd_line.get("verbose", "0");
 *
 *    // Determine if a "z" was provided after a single "-"
 *    // "-z" or "-xyz" or "-xzvf" would all qualify
 *    bool z_flag = cmd_line.exists("z");
 *
 *    return 0;
 * }
 * \endcode
 */
class CommandLine final {

private:
	using key_value_type = std::map<std::string,std::string>;

public:

	/// Empty data structure
	CommandLine() = default;

	/// Copy constructor
	CommandLine(const CommandLine& CL) = default;

	/// Move constructor
	CommandLine(CommandLine&& CL) = default;

	/// Construct with command line arguments
	CommandLine(int argc, char* argv[]);

	/// Destroy the class
	~CommandLine() = default;

	/// Assign the class
	CommandLine& operator=(const CommandLine& CL) = default;

	/// Move the class
	CommandLine& operator=(CommandLine&& CL) = default;

	void process(int argc, char* argv[]);

	/// Insert a key-value pair into the data
	void insert(const std::string& key, const std::string& value = "");

	/// Return true if the key exists within the command line arguments
	bool exists(const std::string& key) const;

    /// Return true if the key exists within the command line arguments
    bool exists(const std::string& short_key, const std::string& long_key) const;

	/// Returns the value following the search key
	/*
	 * Returns the value that followed the key within the argument list.
	 * If the key existed within the arguments but had no value the method
	 * will return an empty result even if a default is provided.  The default
	 * is returned only when they key did not exist on the command line.
	 */
	std::string get(const std::string& key, const std::string& def = "");

	std::string get(const std::string& short_key, const std::string& long_key, const std::string& def = "");

	template<typename OStreamType>
	OStreamType& display(OStreamType& sout);

private:
	key_value_type key_value_pairs_;

	// Test if key exists
	bool is_key_(const std::string& str) const;

	// Test if string is a single "-" key
	bool is_short_key_(const std::string& str) const;

	// Strip leading "-" from keys
	std::string get_key_(const std::string& str) const;

	// Generate key-value pairs from a tagged "-" or "--" key
	key_value_type generate_pairs_(const std::string& tagged_key,
			                       const std::string& value) const;
};



template<typename OStreamType>
OStreamType& CommandLine::display(OStreamType& sout){
	sout << "Parsed Command Line:\n";
	for(auto keyval : key_value_pairs_){
		sout << std::get<0>(keyval) << "  " << std::get<1>(keyval) << "\n";
	}
	return sout;
}



} // namespace tbox
} // namespace geoflow


#endif /* SRC_XSTD_COMMAND_LINE_HPP_ */
