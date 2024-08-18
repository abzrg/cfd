#include "utils/linspace.hpp"

#include <algorithm> // std::generate
#include <cstddef>
#include <tuple>
#include <vector>

#ifdef SIMPLE

std::vector<double> utils::linspace(double start, double end, int num)
{
    std::vector<double> values(num);
    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i)
    {
        values[i] = start + i * step;
    }

    return values;
}
#else
// Function to create a vector with evenly spaced values
std::vector<double> utils::linspace(double start, double end, size_t num)
{
    std::vector<double> values(num);
    const auto step = (end - start) / (static_cast<double>(num) - 1);

    // Generate values using a lambda function
    std::generate(values.begin(), values.end(),
                  [start, step, i = 0]() mutable { return start + i++ * step; });

    return values;
}

// Override: Gets a interval tuple as input
std::vector<double> utils::linspace(std::tuple<double, double> interval, size_t num)
{
    auto [start, end] = interval;
    std::vector<double> values(num);
    const auto step = (end - start) / (static_cast<double>(num) - 1);

    std::generate(values.begin(), values.end(),
                  [start, step, i = 0]() mutable { return start + i++ * step; });

    return values;
}
#endif
