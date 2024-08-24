#pragma once

#include <armadillo>
#include <string>
#include <vector>

namespace utils {

void write1d(const std::vector<double> &u, double timestep);
void write1d(const std::vector<double> &u, const std::string &fpath);

void write2d(const arma::mat &u, double time);

} // namespace utils
