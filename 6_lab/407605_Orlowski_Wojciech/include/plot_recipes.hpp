#ifndef PLOT_RECIPES_HPP
#define PLOT_RECIPES_HPP

#include <vector>
#include <fstream>
#include <iostream>

void plot_2D(std::vector<double>, std::vector<double>, std::string);
void plot_D(std::vector<double>& , std::vector<double>& , std::vector<double>& , std::vector<double>&, std::string);
void plot_2D_closed(std::vector<double>& x, std::vector<double>& y, std::string png_name, double Ra);
void plot_n(std::vector<int>&, std::string png_name);
void plot_2D_different_times(std::vector<double>&,
                             std::vector<double>&,
                             std::vector<double>&,
                             std::vector<double>&,
                             std::vector<double>&,
                             std::vector<double>&,
                             std::string);

#endif