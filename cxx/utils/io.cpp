#include "utils/io.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

void utils::write1d(const std::vector<double> &u, double time)
{
    std::ofstream outfile("timestep_" + std::to_string(time) + ".txt");

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
        std::cerr << "Error opening file for timestep " << time << std::endl;
    }
}

void utils::write1d(const std::vector<double> &u, const std::string &fpath)
{
    std::ofstream outfile(fpath);

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
        std::cerr << "Error opening file '" << fpath << "'." << std::endl;
    }
}

void utils::write2d(const arma::mat &u, double time)
{
    std::stringstream filename;
    filename << "timestep_" << std::fixed << std::setprecision(2) << time << ".txt";

    std::ofstream file(filename.str());

    if (file.is_open())
    {
        u.save(file, arma::raw_ascii);
        file.close();
        std::cout << "Written timestep at time " << time << " to " << filename.str() << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file " << filename.str() << std::endl;
    }
}
