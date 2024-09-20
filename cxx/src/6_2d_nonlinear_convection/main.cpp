#include "utils/io.hpp"

#include <armadillo>

void initialize(arma::mat &field);
void nonlinearConvection(arma::mat &u, arma::mat &v, double nt, double dt, double dx, double dy);

int main()
{
    constexpr auto nx = 81;
    constexpr auto dx = 2.0 / (nx - 1);

    constexpr auto ny = 81;
    constexpr auto dy = 2.0 / (ny - 1);

    constexpr auto sigma = 0.2;

    constexpr auto nt = 100;
    constexpr auto dt = sigma * dx;

    arma::mat u = arma::ones<arma::mat>(nx, ny);
    arma::mat v{u};

    initialize(u);
    initialize(v);

    nonlinearConvection(u, v, nt, dt, dx, dy);

    return 0;
}

void initialize(arma::mat &field)
{
    size_t ny = field.n_rows;
    size_t nx = field.n_cols;

    for (size_t i = nx / 4; i <= nx / 2; ++i)
    {
        for (size_t j = ny / 4; j <= ny / 2; ++j)
        {
            field(i, j) = 2.0;
        }
    }
}

void nonlinearConvection(arma::mat &u, arma::mat &v, double nt, double dt, double dx, double dy)
{
    size_t ny = u.n_rows;
    size_t nx = u.n_cols;

    arma::mat u_old{};
    arma::mat v_old{};

    for (int n = 0; n < nt; ++n)
    {
        u_old = u;
        v_old = v;

        for (size_t i = 1; i < nx - 1; ++i)
        {
            for (size_t j = 1; j < ny - 1; ++j)
            {
                u(i, j) = u_old(i, j) - u_old(i, j) * (dt / dx) * (u_old(i, j) - u_old(i - 1, j)) -
                          v_old(i, j) * (dt / dy) * (u_old(i, j) - u_old(i, j - 1));

                v(i, j) = v_old(i, j) - u_old(i, j) * (dt / dx) * (v_old(i, j) - v_old(i - 1, j)) -
                          v_old(i, j) * (dt / dy) * (v_old(i, j) - v_old(i, j - 1));
            }
        }

        const auto current_time = n * dt;
        utils::write2d(u, current_time);
    }
}
