#include "utils/io.hpp"
#include "utils/linspace.hpp"

#include <algorithm> // std::transform
#include <cstdlib>   // size_t
#include <fstream>
#include <iomanip> // std::setprecision
#include <iostream>
#include <tuple> // std::tuple, std::make_tuple
#include <vector>

void initialize(std::vector<double> &u, double u_min, double u_max,
                std::tuple<double, double> interval);
void linearConvection(std::vector<double> &u, double c, double dt, double dx, unsigned int nt);
void updateBoundaryCondition(std::vector<double> &u, double c, double dt, double dx,
                             unsigned int nt);

int main()
{
    constexpr auto x_domain_interval = std::make_tuple(0.0, 2.0);
    constexpr auto nx = 41UL;
    constexpr auto dx = (std::get<1>(x_domain_interval) - std::get<0>(x_domain_interval)) /
                        (static_cast<double>(nx) - 1.0);
    constexpr auto nt = 25UL;
    constexpr auto dt = 0.025;
    constexpr auto c = 1.0;

    std::vector<double> x = utils::linspace(x_domain_interval, nx);

    // Initialize scalar field
    std::vector<double> u(nx, 1.0);
    initialize(u, 1.0, 2.0, x_domain_interval);

    linearConvection(u, c, dt, dx, nt);
    updateBoundaryCondition(u, c, dt, dx, nt);

    return EXIT_SUCCESS;
}

void initialize(std::vector<double> &u, double u_min, double u_max,
                std::tuple<double, double> interval)
{
    size_t nx = u.size();
    const auto [x_min, x_max] = interval;
    const double dx = (x_max - x_min) / (static_cast<double>(nx) - 1.0);

    for (int i = 0; i < nx; ++i)
    {
        const double x = x_min + i * dx;
        u[i] = (x >= 0.5 && x <= 1.0) ? u_max : u_min;
    }
}

void linearConvection(std::vector<double> &u, double c, double dt, double dx, unsigned int nt)
{
    const auto nx = u.size();

    std::vector<double> un(nx);

    // Temporal iteration
    for (size_t t = 0; t < nt; ++t)
    {
        utils::write1d(u, static_cast<double>(t));
        std::copy(u.begin(), u.end(), un.begin());

        // Spatial iteration over internal nodes
#ifdef SIMPLE
        for (size_t i = 1; i < nx - 1; ++i)
        {
            u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1]);
        }
#else
        std::transform(u.begin() + 1, u.end() - 1, u.begin() + 1,
                       [&un, c, dt, dx, index = 1](double) mutable {
                           double result = un[index] - c * dt / dx * (un[index] - un[index - 1]);
                           ++index;
                           return result;
                       });
#endif
    }
}

void updateBoundaryCondition(std::vector<double> &u, double c, double dt, double dx,
                             unsigned int nt)
{
}
