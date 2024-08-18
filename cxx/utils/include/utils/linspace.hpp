#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

namespace utils {

std::vector<double> linspace(double start, double end, size_t num);
std::vector<double> linspace(std::tuple<double, double> interval, size_t num);

} // namespace utils
