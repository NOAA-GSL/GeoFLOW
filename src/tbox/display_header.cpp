/*
 * display_header.cpp
 *
 *  Created on: Oct 3, 2019
 *      Author: bryan.flynt
 */



#include "configure.hpp"
#include "tbox/display_header.hpp"
#include "tbox/pio.hpp"

#include <algorithm>

namespace geoflow {
namespace tbox {

void
display_header(const std::string& code_name){
    const int field_width = 72;
    const int name_width = code_name.size();
    const int npad = std::max(0, (field_width - name_width) / 2);

    // Define some common string types
    const std::string bar(field_width, '=');
    const std::string pad(npad, ' ');

    pio::pout << "\n";
    pio::pout << bar << "\n";
    pio::pout << pad << code_name << "\n";
    pio::pout << bar << "\n";
    pio::pout << "Version:      " << PROJECT_VERSION << "\n";
    pio::pout << "GIT Hash:     " << GIT_LONG_HASH   << "\n";
    pio::pout << "GIT Time:     " << GIT_COMMIT_TIME << "\n";
    pio::pout << "Compile Time: " << __DATE__ << " " << __TIME__ << "\n";
    pio::pout << bar << std::endl;
    pio::pout << "\n";
}


} // namespace tbox
} // namespace geoflow


