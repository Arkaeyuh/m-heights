#include "utils.h"
#include <sstream>

// explicit instantiations you need (otherwise linker errors)
template std::string vec_to_string<int>(const std::vector<int>&);
template std::string vec_to_string<double>(const std::vector<double>&);
