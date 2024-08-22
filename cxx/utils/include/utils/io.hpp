#pragma once

#include <string>
#include <vector>

namespace utils {

void write1d(const std::vector<double> &u, double timestep);
void write1d(const std::vector<double> &u, const std::string &fpath);

} // namespace utils
