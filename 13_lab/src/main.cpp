#include <iostream>
#include "dsmc.hpp"
#include <vector>


int main (){
	/*
	std::string out_dir =  "data/";
	std::string input = "inputs/input_left.in";

	DSMC_2D ob ;

	ob.read (input.c_str()); // wczytujemydane zpliku wejściowego    
	ob.init (); // automatyczna inicjalizacja położeń i prędkości
	ob.write_position_velocity ((out_dir+"left.dat").c_str()); // zapis ustawień początkowych

	// input = "inputs/input_right.in";

	// DSMC_2D ob_right;

	// ob_right.read(input.c_str());
	// ob_right.init();
	// ob_right.write_position_velocity ((out_dir+"right.dat").c_str());
	*/
	 
	std::string out_dir = "data/";
	std::string input = "inputs/input.in";

	DSMC_2D ob;

	ob.read(input.c_str());
	ob.init();
	ob.write_position_velocity((out_dir+"initial_pos.dat").c_str());

	ob.nthreads = 4;
	ob.icol =1; // cząstki zderzają się
	ob.evolution (0.0 , 2000);

	ob.write_position_velocity((out_dir+"rv_end.dat").c_str()); // zapis położeń i prędkości końcowych do pliku


	/*
	if(system(("mv wyniki "+out_dir).c_str())){
		std::cout << "wyniki couldnt be moved" << std::endl;
	}
	if(system(("mv xy_1.dat "+out_dir).c_str())){
		std::cout << "xy_1.dat couldnt be moved" << std::endl;
	}
	*/
	
  return 0;
}
