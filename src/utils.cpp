#include "utils.h"
#include <sstream>

// general case: T must have operator<<
template<typename T>
std::string vec_to_string(const std::vector<T>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ", ";
        os << v[i];
    }
    os << "]";
    return os.str();
}

// specialization for vector<int> elements (i.e. for std::vector<std::vector<int>>)
template<>
std::string vec_to_string<std::vector<int>>(const std::vector<std::vector<int>>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ", ";
        // recursive call to the T=int version
        os << vec_to_string<int>(v[i]);
    }
    os << "]";
    return os.str();
}

// now explicitly instantiate only those two
template std::string vec_to_string<int>(const std::vector<int>&);
template std::string vec_to_string<std::vector<int>>(
    const std::vector<std::vector<int>>&);
