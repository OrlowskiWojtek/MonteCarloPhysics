#ifndef DSMC_HPP
#define DSMC_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <time.h>
#include <vector>
#include <iomanip> 
#include <string>
#include <cstdlib>
#include <filesystem>
#include <omp.h>


using namespace std;
/* 
 * DSMC_2D: Direct Simulation Monte Carlo  w 2D 
 * 
 * klasa PARTICLE: tu trzymamy dane pojedynczej czastki
 * (x,y), (vx,vy), v=sqrt(vx^2+vy^2), mc-masa,rc-promien efektywny
 * ic-rodzaj czastki, (ix,iy) - polozenie komorki, irun-flaga dla wykonanego ruchu
 * 
 * 
 * zmienne:
 * xmin,xmax,ymin,ymax - wymiary pudla obliczeniowego
 * nx,ny,nxy=nx*ny - liczba komorek w kierunkach x/y
 * delta_x,delta_y - wymiary pojedynczej komorki
 * temp - temperatura
 * ntot - calkowita liczba czastek
 * nc[],mc[],rc[] - tablice zawierajace: liczbe czastek,mase i promien efektywny danego rodzaju czastek
 * dt - aktualny krok czasowy
 * time_sum - czas trwania symulacji 
 * 
 * par[ntot] - tablica obiektow typu PARTICLE
 * 
 * edge_out[] - tablica wierzcholkow pudla obliczeniowego (zapisane zgodnie z ruchem wskazowek zegara)
 * nodes_out - liczba wierzcholkow pudla
 * edge[] - tablica wierzcholkow barier wewnatrz pudla (zapisane przeciwnie do ruchu wskazowek zegara)
 * nodes - liczba wierzcholkow dla bariery wewnetrznej
 * nthreads - liczba uzywanych watkow OPENMP
 * 
 * indx,indx0 - tablice do sortowania czastek wedlug numerow komorek
 * 
 */

class DSMC_2D
{
	
	
private: class PARTICLE{
		public:
			double x,y,vx,vy,v,mc,rc;
			double x0,y0,vx0,vy0;
			double vxt=0.;
			double vyt=0.;
			double cv;
			int ic,ix,iy,irun;
			int indx=-1;
			int ncol=0;
			int nbound_col=0;
			double path=0.;
		};

	
public:
		double kb=1.38E-23;
		double xmin,xmax,ymin,ymax;
		int nx,ny,nxy,n_mix;
		double temp;
		double tempi[20];
		
		vector<vector<double>> temp_cell;
		vector<vector<double>> dens_cell;
		vector<vector<double>> vx_cell;
		vector<vector<double>> vy_cell;
		vector<vector<double>> press_cell;
		vector<vector<vector<double>>> press_tens;
		vector<double> velocity_x;
		vector<double> velocity_y;
		vector<double> press_x;
		vector<double> press_y;
		vector<int> indx;
		vector<vector<int>> indx0;
		vector<PARTICLE> par;
		
		int nc[5];
		int ntot;
		double mc[5],rc[5];
		int init_dist;
		double dt,time_sum;
		double delta_x,delta_y;
		

		int nodes_out=4;
		int nodes;
		double edge_out[100][2];
		double edge[100][2];	
		int icol=1;
		int nthreads=1;
		
		double l_mfp_num,l_mfp_teo,l_mfp_enskog, gamma_col;
		double cv_coeff=0.;
		double diffusion_coeff=0.;
		
	
		
		unsigned seed=12345;
		default_random_engine generator; //typ generatora 
		uniform_real_distribution<double> rng_U= uniform_real_distribution<double>(0.,1.);
		normal_distribution<double> rng_N=normal_distribution<double>(0.,1.);

		
/******* funkcje *****************************/		
		
	void read(const char * );
	void init();
	void hist_velocity_all(const char *,double, int);
	void write_nptv(const char *,int);
	void evolution(double,int);	
	void write_position_velocity(const char * );
	void write_cell_bounds(const char * );		
	double temperature_eff();
	void temperature_cell();
	void density_cell();
	void velocity_cell();
	void mean_free_path();
	void autocorrelation();
	void cell_ix_iy(double x,double y, int *ix,int *iy);
	void compute_pressure();
	int  global_index_cell(int,int);
	
	void step();	
	void sort();			
	void init_distribution();	
	int  intersection(double,double,double,double,double,double,double,double,double *,double *,double *,double *,int *);
	int  barrier_reflection(PARTICLE &,double *);
	void free_flight(PARTICLE &,double *);
	int  particle_collisions(int,int);
	void particle_particle_scattering_balls(PARTICLE &, PARTICLE &,int);
	void maxwell_distribution(double *, double *, double,double);
	void compute_pressure_tensor_interaction(PARTICLE &, PARTICLE &,PARTICLE &, PARTICLE &);
	void compute_pressure_tensor_kinetic();
	void test();
	
};

#endif