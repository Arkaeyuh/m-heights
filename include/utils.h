#pragma once
#include <string>
#include <vector>

// stringify a std::vector<T>
template<typename T>
std::string vec_to_string(const std::vector<T>& v);
