/*
 * command_line.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: bryan.flynt
 */

#ifndef SRC_TBOX_COMMAND_LINE_HPP_
#define SRC_TBOX_COMMAND_LINE_HPP_

#include <string>
#include <vector>

namespace geoflow {
namespace tbox {

/**
 * \brief
 * Class CommandLine parses simple input line commands from Linux terminals.
 *
 * \details
 * The class parses the command line arguments into string tokens.  When a user
 * requests the value of that tag it searches the tokens and returns the next
 * value within the list of tokens.  This simple parser does NOT accept multiple
 * option tags within the same option.
 *
 * To get a single value following "-i" use:
 * @code{.cpp}
 * int main(int argc, char* argv[]){
 *
 *     std::string file_name = "default.txt";
 *	   CommandLine cline(argc,argv);
 *	   if( cline.tagExists("-i") ){
 *	       file_name = cline.getValue("-i");
 *	   }
 *
 *	   // Do something with file
 *
 *	   return 0;
 *	}
 * @endcode
 *
 */
class CommandLine {
public:

	/**
	 * Can not construct a CommandLine without arguments
	 */
	CommandLine() = delete;

	/**
	 * Construct with command line arguments.
	 */
	CommandLine(int argc, char* argv[]);

	/**
	 * Copy the class
	 */
	CommandLine(const CommandLine& CL);

	/**
	 * Destroy the class
	 */
	~CommandLine();

	/**
	 * Copy and swap the class
	 */
	CommandLine& operator=(CommandLine CL);

	/**
	 * Search if the provided tag exists within the command line arguments
	 */
	bool tagExists(const std::string& tag) const;

	/**
	 * Get the value following the provided tag.
	 *
	 * Returns the single value following the requested search tag.
	 * If the tag is not found or the tag exists but no value follows it
	 * then an empty string is returned.
	 */
	std::string getValue(const std::string& tag) const;


private:
        std::vector<std::string> tokens_;
};


} // namespace tbox
} // namespace geoflow


#endif /* SRC_TBOX_COMMAND_LINE_HPP_ */
