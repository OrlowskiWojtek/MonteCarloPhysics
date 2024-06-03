#include <iostream>
#include "dsmc.hpp"
#include <vector>


int main (){

	for(int i = 5 ; i <= 8; i++){
		std::string out_dir =  "tasks/task_"+std::to_string(i);
		std::string input = "inputs/input_"+std::to_string(i)+".in";

    DSMC_2D ob ;

    ob.read (input.c_str()); // wczytujemydane zpliku wejściowego    
    ob.init (); // automatyczna inicjalizacja położeń i prędkości
    ob.hist_velocity_all ((out_dir+"/hist_begin.dat").c_str() ,5.0 ,50); // zapis histogramu prędkosci do pliku
    
    //ob.write_position_velocity ((out_dir+"/rv_start.dat").c_str()); // zapis ustawień początkowych
    //ob.nthreads = 2;
    //ob.icol =1; // cząstki zderzają się
    //ob.evolution (0.0 , 500); 
    //ob.write_position_velocity((out_dir+"/rv.dat").c_str()); // zapis położeń i prędkości końcowych do pliku

		/*
		if(system(("mv wyniki "+out_dir).c_str())){
			std::cout << "wyniki couldnt be moved" << std::endl;
		}
		if(system(("mv xy_1.dat "+out_dir).c_str())){
			std::cout << "xy_1.dat couldnt be moved" << std::endl;
		}
		*/
	}
  return 0;
}
