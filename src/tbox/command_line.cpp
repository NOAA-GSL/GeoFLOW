/*
 * command_line.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: bryan.flynt
 */

#include <algorithm>
#include "tbox/assert.hpp"
#include "tbox/command_line.hpp"

#include <algorithm>

namespace geoflow {
namespace tbox {


CommandLine::~CommandLine(){
}

CommandLine::CommandLine(int argc, char* argv[]){
	for (int i = 1; i < argc; ++i) {
		this->tokens_.push_back(std::string(argv[i]));
    }
}

CommandLine::CommandLine(const CommandLine& CL){
	tokens_ = CL.tokens_;
}

CommandLine& CommandLine::operator=(CommandLine CL){
	tokens_.swap(CL.tokens_);
	return *this;
}

bool
CommandLine::tagExists(const std::string& tag) const{
	bool ans = false;
	auto loc = std::find(tokens_.begin(), tokens_.end(), tag);
	if( loc != tokens_.end() ){
		ans = true;
	}
	return ans;
}

std::string
CommandLine::getValue(const std::string& tag) const{
	std::string ans("");
	auto loc = std::find(tokens_.begin(), tokens_.end(), tag);
	if (loc != tokens_.end() && ++loc != tokens_.end()){
		ans = *loc;
    }
	return ans;
}


} // namespace tbox
} // namespace geoflow
