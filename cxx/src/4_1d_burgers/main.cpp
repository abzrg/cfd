#include "utils/io.hpp"
#include "utils/linspace.hpp"
#include "utils/mathematicalConstants.hpp"

#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <vector>

void initialize(std::vector<double> &u, const std::vector<double> &x, const double t,
                const double nu);
void burgers(std::vector<double> &u, double nu, double dt, double dx, unsigned int nt);

int main()
{
    constexpr auto pi = utils::constants::mathematical::pi;
    constexpr auto x_domain_interval = std::make_tuple(0.0, 2.0 * pi);
    constexpr auto nx = 101UL;
    constexpr auto dx = (std::get<1>(x_domain_interval) - std::get<0>(x_domain_interval)) /
                        (static_cast<double>(nx) - 1.0);
    constexpr auto nt = 100UL;
    constexpr auto nu = 0.07;
    constexpr auto dt = dx * nu;
    const auto t = 0.0;

    std::vector<double> x = utils::linspace(x_domain_interval, nx);

    // Initialize scalar field
    std::vector<double> u(nx, 1.0);
    initialize(u, x, t, nu);
    utils::write1d(u, "U0");

    burgers(u, nu, dt, dx, nt);

    return EXIT_SUCCESS;
}

void initialize(std::vector<double> &u, const std::vector<double> &x, const double t,
                const double nu)
{
    size_t nx = u.size();
    const auto x_max = x.back();
    const auto x_min = x.front();
    const double dx = (x_max - x_min) / (static_cast<double>(nx) - 1.0);
    constexpr auto pi = utils::constants::mathematical::pi;

    for (size_t i = 0; i < nx; ++i)
    {
        // clang-format off
        u[i] = -2 * nu *
               (
                   -(-8*t + 2*x[i]) * exp(-pow((-4*t + x[i]), 2) / (4*nu * (t + 1))) / (4*nu * (t + 1))
                   -
                   (-8*t + 2*x[i] - 4*pi) * exp(-pow((-4*t + x[i] - 2*pi), 2) / (4*nu * (t + 1))) / (4*nu * (t + 1))
               )
               /
               (
                   exp(-pow((-4*t + x[i] - 2*pi), 2) / (4*nu * (t + 1)))
                   +
                   exp(-pow((-4*t + x[i]), 2) / (4*nu * (t + 1)))
               )
               + 4;
        // clang-format on
    }
}

void burgers(std::vector<double> &u, double nu, double dt, double dx, unsigned int nt)
{
    const auto nx = u.size();

    std::vector<double> un(nx);

    // Temporal iteration
    for (size_t t = 0; t < nt; ++t)
    {
        utils::write1d(u, static_cast<double>(t));
        std::copy(u.begin(), u.end(), un.begin());

        // Spatial iteration
#ifdef SIMPLE
        for (size_t i = 1; i < nx - 1; ++i)
        {
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) +
                   nu * dt / dx / dx * (un[i + 1] - 2 * un[i] + un[i - 1]);
        }
#else
        std::transform(u.begin() + 1, u.end() - 1, u.begin() + 1,
                       [&un, nu, dt, dx, i = 1](double) mutable {
                           double res = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) +
                                        nu * dt / dx / dx * (un[i + 1] - 2 * un[i] + un[i - 1]);
                           ++i;
                           return res;
                       });
#endif
        // Update boundary conditions
        u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 2]) +
               nu * dt / dx / dx * (un[1] - 2 * un[0] + un[nx - 2]);
        u[nx - 1] = u[0];
    }
}
