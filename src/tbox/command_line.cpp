/*
 * command_line.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: bryan.flynt
 */



#include "tbox/command_line.hpp"


namespace geoflow {
namespace tbox {


CommandLine::CommandLine(int argc, char* argv[]){
    this->process(argc,argv);
}

void
CommandLine::process(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        if( this->is_key_(std::string(argv[i])) ){
            std::string tagged_key = std::string(argv[i]);

            std::string value;
            if( ((i+1) < argc) && !this->is_key_(std::string(argv[i+1])) ){
                value = std::string(argv[i+1]);
                ++i;
            }
            auto kvpairs = generate_pairs_(tagged_key,value);
            key_value_pairs_.insert(kvpairs.begin(),kvpairs.end());
        }
    }
}

void
CommandLine::insert(const std::string& key, const std::string& value){
    key_value_pairs_[key] = value;
}


bool
CommandLine::exists(const std::string& key) const{
    return (key_value_pairs_.count(key) > 0);
}

bool
CommandLine::exists(const std::string& short_key, const std::string& long_key) const{
    return ( this->exists(short_key) || this->exists(long_key) );
}

std::string
CommandLine::get(const std::string& key, const std::string& def) {
    std::string ans(def);
    if( this->exists(key) ){
        ans = key_value_pairs_[key];
    }
    return ans;
}

std::string
CommandLine::get(const std::string& short_key,
                 const std::string& long_key,
                 const std::string& def) {
    std::string ans(def);
    if( this->exists(short_key) ){
        ans = key_value_pairs_[short_key];
    }
    else if( this->exists(long_key) ){
        ans = key_value_pairs_[long_key];
    }
    return ans;
}

bool
CommandLine::is_key_(const std::string& str) const {
    return (str.find("-") == 0);
}

bool
CommandLine::is_short_key_(const std::string& str) const {
    return (str.rfind("-") == 0);
}

std::string
CommandLine::get_key_(const std::string& str) const{
    return str.substr(str.rfind("-")+1);
}

typename CommandLine::key_value_type
CommandLine::generate_pairs_(const std::string& tagged_key,const std::string& value)const{
    key_value_type ans;
    auto key = get_key_(tagged_key);
    if( is_short_key_(tagged_key) ){
        for(auto it = key.begin(); it != key.end(); ++it){
            ans[std::string(1,*it)] = value;
        }
    }
    else {
        ans[key] = value;
    }
    return ans;
}

} // namespace tbox
} // namespace geoflow
