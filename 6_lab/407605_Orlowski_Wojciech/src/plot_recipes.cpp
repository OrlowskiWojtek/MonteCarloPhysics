#include "plot_recipes.hpp"

void plot_2D(std::vector<double> x, std::vector<double> y, std::string png_name)
{
    if(x.size() != y.size())
        return;

    std::string filename("tmp_filename");
    std::ofstream tmp;
    tmp.open(filename);

    for(uint64_t i = 0; i < x.size(); i++){
        tmp << x[i] << "\t" << y[i] << "\n"; 
    }

    tmp.close();


    if(system(("gnuplot -e 'outputName = \""+png_name+"\"' plot_scripts/2d.gnu").c_str()))
        std::cout << "Cannot plot 2d file" << std::endl;
    
    if(system("rm tmp_filename"))
        std::cout << "Cannot remove tmp_filename" << std::endl;
}

void plot_D(std::vector<double>& time, std::vector<double>& Dxx, std::vector<double>& Dyy, std::vector<double>& Dxy, std::string png_name)
{
    std::string filename("tmp_filename");
    std::ofstream tmp;
    tmp.open(filename);

    for(uint64_t i = 0; i < time.size(); i++){
        tmp << time[i] << "\t" << Dxx[i] << "\t" << Dyy[i] << "\t" << Dxy[i] << "\n"; 
    }

    tmp.close();

    if(system(("gnuplot -e 'outputName = \""+png_name+"\"' plot_scripts/D.gnu").c_str()))
        std::cout << "Cannot plot 2d file" << std::endl;
    
    if(system("rm tmp_filename"))
        std::cout << "Cannot remove tmp_filename" << std::endl;
}

void plot_2D_closed(std::vector<double>& x, std::vector<double>& y, std::string png_name, double Ra)
{
    if(x.size() != y.size())
        return;

    std::string filename("tmp_filename");
    std::ofstream tmp;
    tmp.open(filename);

    for(uint64_t i = 0; i < x.size(); i++){
        tmp << x[i] << "\t" << y[i] << "\n"; 
    }

    tmp.close();


    if(system(("gnuplot -e 'Ra = "+std::to_string(Ra)+";outputName = \""+png_name+"\"' plot_scripts/2d_closed.gnu").c_str()))
        std::cout << "Cannot plot 2d file" << std::endl;
    
    if(system("rm tmp_filename"))
        std::cout << "Cannot remove tmp_filename" << std::endl;
}

void plot_n(std::vector<int>& x, std::string png_name)
{

    std::string filename("tmp_filename");
    std::ofstream tmp;
    tmp.open(filename);

    for(uint64_t i = 0; i < x.size(); i++){
        tmp << x[i] << "\n"; 
    }

    tmp.close();

    if(system(("gnuplot -e 'outputName = \""+png_name+"\"' plot_scripts/n.gnu").c_str()))
        std::cout << "Cannot plot 2d file" << std::endl;
    
    if(system("rm tmp_filename"))
        std::cout << "Cannot remove tmp_filename" << std::endl;
}

void plot_2D_different_times(std::vector<double>& x1,
                             std::vector<double>& y1,
                             std::vector<double>& x2,
                             std::vector<double>& y2,
                             std::vector<double>& x3,
                             std::vector<double>& y3,
                             std::string png_name)
{
    std::string filename("tmp_filename");
    std::ofstream tmp;
    tmp.open(filename);

    for(uint64_t i = 0; i < x1.size(); i++){
        tmp << x1[i] << "\t" << y1[i] << "\t" <<x2[i] << "\t" << y2[i]<< "\t" << x3[i] << "\t" << y3[i] << "\n"; 
    }

    tmp.close();


    if(system(("gnuplot -e 'outputName = \""+png_name+"\"' plot_scripts/2d_times.gnu").c_str()))
        std::cout << "Cannot plot 2d file" << std::endl;
    
    if(system("rm tmp_filename"))
        std::cout << "Cannot remove tmp_filename" << std::endl;
}