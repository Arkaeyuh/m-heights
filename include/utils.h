#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <sstream>
#include <type_traits>

// Forward-declare so recursive call works
template<typename T>
std::string vec_to_string(const std::vector<T>& v);

// arithmetic -> to_string
template<typename U>
std::enable_if_t<std::is_arithmetic<U>::value, std::string>
element_to_string(const U &x) {
    return std::to_string(x);
}

// string stays as-is
inline std::string element_to_string(const std::string &s) {
    return s;
}

// nested vector -> recursion
template<typename U>
std::string element_to_string(const std::vector<U> &x) {
    return vec_to_string(x);
}

// finally the generic vector printer
template<typename T>
std::string vec_to_string(const std::vector<T>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ", ";
        os << element_to_string(v[i]);
    }
    os << "]";
    return os.str();
}

#endif // UTILS_H
