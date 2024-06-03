#ifndef PHOTON_DIFFUSION_H
#define PHOTON_DIFFUSION_H

#include<vector>
#include<string>
using namespace std;


class PHOTON_DIFFUSION_2D
{
	
	private: class BEAM{
			public:
				double w;
				double x,y; 		//polozenia
				double x_new,y_new; 	//nowe proponowanepolozenia
				double rx,ry; 		//kierunek wiazki
				int layer; 			//nr warstwy
				bool alive=true;		//
				int length;			//flaga: length=1 sprawdzamy warunki i przemieszczamy wiazke  
				vector<double> path;    //rejestrujemy sciezke 
			
	};
	
	
	public: 	
		//zmienne
			int nlayers;
			int nx,ny;
			double xmax,ymax,dx,dy;
			double x_source,dx_source;
			double x_detect,dx_detect;
			double p_min;
			double w_min;
			double rx0,ry0;
			
			
			int write_all_paths;
			int write_source_detection_paths;
			
		//wiazka
			BEAM beam;
			
			
			double abs_specular;
			
		// tablice
			vector<vector<double>> absorption;
			vector<double> reflectance;
			vector<double> transmittance;
			vector<vector<double>> layers_data;
			
		//funkcje	
			PHOTON_DIFFUSION_2D();
			void init();
			double uniform();
			void single_path();
			void roulette();
			void calculate_new_position();
			void scatter_in_layer(); //
			void scatter_up_down_boundary();
			double sign(double);
			void segment_intersection(double x1,double y1, double x2,double y2,double x3,double y3, double x4,double y4,
								 double & x_cross,double & y_cross, int & icross);
			void write_paths_to_file();
            void fill_m_layer(int m, double a, double s, double d, double g, double n);
			void export_to_files(string name);
};


#endif