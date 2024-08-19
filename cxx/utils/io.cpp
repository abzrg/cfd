#include "utils/io.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

void utils::write1d(const std::vector<double> &u, double timestep)
{
    std::ofstream outfile("timestep_" + std::to_string(timestep) + ".txt");

    if (outfile.is_open())
    {
        for (const auto &value : u)
        {
            outfile << std::fixed << std::setprecision(2) << value << "\n";
        }
        outfile.close();
    }
    else
    {
        std::cerr << "Error opening file for timestep " << timestep << std::endl;
    }
}
