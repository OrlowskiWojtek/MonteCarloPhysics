#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<vector>
#include<complex>

#include "jlcxx/jlcxx.hpp"

using namespace std;

/*****************************************************************************
 *  liczymy rozwiazanie dokladne dla impulsu generowanego przez zrodlo
 *  oznaczenia:
 * 			x -polozenie, t - aktualny czas, t0 - maksimum sygnalu zrodla, Omega=2*Pi*f - (f to czestotliwosc zrodla)
 *                sigma - szerokosc impulsu,  R - opor linii, G - konduktancja linii, C-pojemnosc linii, Rg - opor zrodla
 *                Rl-opor cewki w x=l,length-dlugosc linii, number_nodes - liczba wezlow w calkowaniu po czestosci (omedze),
 *			n_sum_terms - maksymalna liczba wyrazow uwzgledniana w sumowaniu wkladow ui (i=0,1,...,n_sum_terms)
 * 
 * 
 *****************************************************************************/


double u_xt_exact(double x, double t, double t0,  double freq, double sigma, double R, double G, double L, double C, 
			    double Rg, double Rl, double length, int number_nodes, int n_sum_terms){
	
	double omega, domega, omega_min,omega_max,c_speed,aj,Omega;
	double sigma_w=1./sigma;
	double uxt=0.;
	complex<double> vg_m, I, z0, zl, ksi_omega, Gamma_l, Gamma_g,k_omega, u0_x_omega,u_sum_omega, multiplier,single_term, uxt_cmplx;
	int zakres;
	
	Omega=2.*M_PI*freq;
	zakres=4.;
	
	c_speed=1/sqrt(L*C); // predkosc swiatla
	I=0.0+1.0i;
	
	uxt_cmplx=0.+0.i;
	
	//wklady scentrowane w:  +/- Omega
	for(int m=-1;m<=1;m+=2){
		
		omega_min=m*Omega-zakres*sigma_w;
		omega_max=m*Omega+zakres*sigma_w;
		domega=(omega_max-omega_min)/number_nodes;
		
		
		//calkowanie po czestosci -(mala omega)
		for(int j=0;j<=number_nodes;j++){
			omega=omega_min+domega*j;
			vg_m=(-1.0*m)*I*sigma*sqrt(M_PI/2.)*exp(-pow((omega-m*Omega)/sigma_w,2)/2-I*(omega-m*Omega)*t0);
			z0=sqrt((R+I*omega*L)/(G+I*omega*C));
			ksi_omega=z0/(z0+Rg);
			
			zl=Rl+I*omega*0.; //tylko opor na wyjsciu
			
			Gamma_l=(zl-z0)/(zl+z0);
			Gamma_g=(Rg-z0)/(Rg+z0);
			k_omega=-I/c_speed*sqrt((I*omega+R/L)*(I*omega+G/C));
			u0_x_omega=ksi_omega*vg_m*( exp(-I*k_omega*x)+Gamma_l*exp(-I*k_omega*(2*length-x)) );
			
			//sumowanie wkladow: ui_x_omega
			u_sum_omega=0.+0.i;
			multiplier=Gamma_l*Gamma_g*exp(-2.*I*k_omega*length);
			single_term=1.;
			for(int ii=0;ii<=n_sum_terms;ii++){
				if(ii>0)single_term*=multiplier;
				u_sum_omega=u_sum_omega+u0_x_omega*single_term;
			}
			
			if(j==0 || j==number_nodes){
				aj=0.5;
			}else{
				aj=1.0;
			}	
			uxt_cmplx+=u_sum_omega/2./M_PI*exp(I*omega*t)*domega*aj;
			
		}//j
	}//m
	
	
	uxt=real(uxt_cmplx);
	
	
	return uxt;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("u_xt_exact",&u_xt_exact);
}


