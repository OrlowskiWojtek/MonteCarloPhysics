#include "tasks.hpp"
#include <chrono>
#include <iostream>

int main()
{
    srand(time(NULL));
    // all in tasks

    if(system("mkdir plots") || system("mkdir plots/task_2"))
        std::cout << "Can't make directories for plots" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    task_1();
    task_2();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

    std::cout << "Total time used: " << duration.count()/1000000. << " seconds" << std::endl;

    return 0;
}

