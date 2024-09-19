#include "utils/io.hpp"

#include <armadillo>

void initialize(arma::mat &u, int nx, int ny);
void linearConvection(arma::mat &u, double nt, double dt, double nx, double dx, double ny,
                      double dy, double cx, double cy);

int main()
{
    constexpr auto nx = 81;
    constexpr auto dx = 2.0 / (nx - 1);

    constexpr auto ny = 81;
    constexpr auto dy = 2.0 / (ny - 1);

    constexpr auto cx = 1.0;
    constexpr auto cy = 1.0;

    constexpr auto sigma = 0.2;

    constexpr auto nt = 100;
    constexpr auto dt = sigma * dx;

    arma::mat u = arma::ones<arma::mat>(nx, ny);

    initialize(u, nx, ny);

    linearConvection(u, nt, dt, nx, dx, ny, dy, cx, cy);

    return 0;
}

void initialize(arma::mat &u, int nx, int ny)
{
    for (int i = nx / 4; i <= nx / 2; ++i)
    {
        for (int j = ny / 4; j <= ny / 2; ++j)
        {
            u(i, j) = 2.0;
        }
    }
}

void linearConvection(arma::mat &u, double nt, double dt, double nx, double dx, double ny,
                      double dy, double cx, double cy)
{
    arma::mat u_old{};
    for (int n = 0; n < nt; ++n)
    {
        u_old = u;

        for (size_t i = 1; static_cast<double>(i) < nx - 1; ++i)
        {
            for (size_t j = 1; static_cast<double>(j) < ny - 1; ++j)
            {
                u(i, j) = u_old(i, j) - cx * (dt / dx) * (u_old(i, j) - u_old(i - 1, j)) -
                          cy * (dt / dy) * (u_old(i, j) - u_old(i, j - 1));
            }
        }

        const auto current_time = n * dt;
        utils::write2d(u, current_time);
    }
}
