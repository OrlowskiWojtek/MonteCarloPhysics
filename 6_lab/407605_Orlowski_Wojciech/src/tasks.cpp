#include "tasks.hpp"
#include "rand.hpp"
#include "plot_recipes.hpp"
#include "part_trans.hpp"

#include <cmath>

void task_1()
{
    double D = 1.0;
    std::vector<int> N_max{ (int)pow(10,2), (int)pow(10,3), (int)pow(10,4), (int)pow(10,5)};
    double x_start = 0., y_start = 0.;
    double dt = 0.1;
    double t_max = 100;

    int N_it = t_max/dt;
    int N_t = N_it;

    double sigma_dt = sqrt(2.*D*dt);

    for(auto &N: N_max){
        std::vector<double> x_i_t_1(N, x_start);
        std::vector<double> y_i_t_1(N, y_start);
        std::vector<double> x_i_t_5(N, x_start);
        std::vector<double> y_i_t_5(N, y_start);
        std::vector<double> x_i_t_20(N, x_start);
        std::vector<double> y_i_t_20(N, y_start);

        std::vector<double> x_i(N,x_start);
        std::vector<double> y_i(N,y_start);
        
        std::vector<double> D_xx(N_it);
        std::vector<double> D_yy(N_it);
        std::vector<double> D_xy(N_it);
        std::vector<double> time(N_it);

        std::cout << "starting for N = "+std::to_string(N) << std::endl;
        
        double x_exp = 0.0;
        double x_exp_sq = 0.0;
        double y_exp = 0.0;
        double y_exp_sq = 0.0;
        double x_y_exp = 0.0;

        double D_xx_exp = 0.0;
        double D_yy_exp = 0.0;
        double D_xy_exp = 0.0;
        double D_xx_exp_sq = 0.0;
        double D_yy_exp_sq = 0.0;
        double D_xy_exp_sq = 0.0;

        for(int it = 0; it < N_it; it++){
            for(int i =0; i < N; i++){
                std::tuple<double,double> move = Normal(0.,sigma_dt); 
                x_i[i] += std::get<0>(move);
                y_i[i] += std::get<1>(move);
            }
            
            if(it == 10){
                x_i_t_1 = x_i;
                y_i_t_1 = y_i;
            }else if(it == 50){
                x_i_t_5 = x_i;
                y_i_t_5 = y_i;
            }else if(it == 200){
                x_i_t_20 = x_i;
                y_i_t_20 = y_i;
            }

            x_exp = std::accumulate(x_i.begin(),x_i.end(),0.0)/static_cast<double>(N);
            y_exp = std::accumulate(y_i.begin(),y_i.end(),0.0)/static_cast<double>(N);
            x_exp_sq = std::inner_product(x_i.begin(), x_i.end(), x_i.begin(), 0.0)/static_cast<double>(N);
            y_exp_sq = std::inner_product(y_i.begin(), y_i.end(), y_i.begin(), 0.0)/static_cast<double>(N);
            x_y_exp  = std::inner_product(x_i.begin(), x_i.end(), y_i.begin(), 0.0)/static_cast<double>(N);

            
            D_xx[it] = (x_exp_sq - pow(x_exp,2))/(2.*(it+1.)*(t_max/N_it));
            D_yy[it] = (y_exp_sq - pow(y_exp,2))/(2.*(it+1.)*(t_max/N_it));
            D_xy[it] = (x_y_exp - x_exp*y_exp)/(2.*(it+1.)*(t_max/N_it));
            
            time[it] = it*(t_max/N_it);
        }

        D_xx_exp = std::accumulate(D_xx.begin(),D_xx.end(),0.0)/static_cast<double>(N_t);
        D_yy_exp = std::accumulate(D_yy.begin(),D_yy.end(),0.0)/static_cast<double>(N_t);
        D_xy_exp = std::accumulate(D_xy.begin(),D_xy.end(),0.0)/static_cast<double>(N_t);
        
        D_xx_exp_sq = std::inner_product(D_xx.begin(), D_xx.end(), D_xx.begin(), 0.0)/static_cast<double>(N_t);
        D_yy_exp_sq = std::inner_product(D_yy.begin(), D_yy.end(), D_yy.begin(), 0.0)/static_cast<double>(N_t);
        D_xy_exp_sq = std::inner_product(D_xy.begin(), D_xy.end(), D_xy.begin(), 0.0)/static_cast<double>(N_t);

        std::cout << "ending for N = "+std::to_string(N)    << "\n";
        std::cout << "====Diffusion coefficients values===" << "\n";
        std::cout << "D_xx: " << D_xx_exp << " std. err.: " << sqrt((D_xx_exp_sq - pow(D_xx_exp,2.))/static_cast<double>(N_t)) << "\n";
        std::cout << "D_yy: " << D_yy_exp << " std. err.: " << sqrt((D_yy_exp_sq - pow(D_yy_exp,2.))/static_cast<double>(N_t)) << "\n";
        std::cout << "D_xy: " << D_xy_exp << " std. err.: " << sqrt((D_xy_exp_sq - pow(D_xy_exp,2.))/static_cast<double>(N_t)) << "\n";
        std::cout << "====================================" << "\n" << std::endl;

        plot_2D_different_times(x_i_t_1,y_i_t_1,x_i_t_5,y_i_t_5,x_i_t_20,y_i_t_20,"diff_times_N="+std::to_string(N));
        plot_2D(x_i,y_i,"N="+std::to_string(N));
        plot_D(time, D_xx, D_yy, D_xy, "D_N="+std::to_string(N));
    }
}

void task_2()
{
    // setting parameters  
    double D = 1.0;
    double N_max = pow(10,4);
    double dt = 0.1;
    double t_max = pow(10,3);
    
    int N = t_max/dt;
    
    std::vector<double> omega{10.0,50.0,100.0};
    std::vector<double> dn{1.0,5.0,10.0};
    double xr = 0.0; double yr = 0.0; double Rr = 5.0;
    double xa = 3.0; double ya = 0.0; std::vector<double> Ra{0.1,0.5};
    double xs =-4.5; double ys = 0.0;

    std::vector<int> theta(N_max, 0);
    std::vector<double> x(N_max, 0.);
    std::vector<double> y(N_max, 0.);
    std::vector<double> plot_x;
    std::vector<double> plot_y;

    double sigma_dt = sqrt(2.*D*dt);

    // vectors for saving outputs;
    std::vector<int> n_vec(N);

    // variables used in simulation
    double x1, x2, y1, y2, dx, dy, length;

    // test for dn = 0.1
    for(auto &om: omega){
        for(auto &R_a: Ra){
            std::vector<int> theta(N_max, 0);

            for(int it = 0; it < N; it++){
                int n = 0;
                int n_new = 0;
                for(int i = 0; i < N_max; i++){
                    if(theta[i] == 0 && n_new < om*dt){
                        theta[i] = 1;
                        x[i] = xs;
                        y[i] = ys;
                        n_new++;
                    }

                    if(theta[i] == 1){
                        x1 = x[i];
                        y1 = y[i];

                        std::tuple<double,double> move = Normal(0.,sigma_dt);
                        dx = std::get<0>(move);
                        dy = std::get<1>(move);
                        x2 = x1 + dx;
                        y2 = y1 + dy;
                        do{
                            particle_translation(x1,y1,x2,y2,xr,yr,Rr,xa,ya,R_a,theta[i],length);   
                        }while( length > 1e-6);
                        if(theta[i] == 1){
                            n++;
                            x[i] = x2;
                            y[i] = y2; 
                        }
                    }
                }
                if(it == 1 || it == 10 || it == 100 || it == 5000){
                    plot_x.clear();
                    plot_y.clear();
                    for(uint64_t k = 0; k < theta.size(); k++){
                        if(theta[k]){
                            plot_x.push_back(x[k]);
                            plot_y.push_back(y[k]);
                        }
                    }
                    plot_2D_closed(plot_x,plot_y,"task_2/it="+std::to_string((int)it)+"om="+std::to_string((int)om)+"ra="+std::to_string(R_a),R_a);
                }
                n_vec[it] = n;
            }
            plot_n(n_vec, "task_2/n_for_om="+std::to_string((int)om)+"ra="+std::to_string(R_a));
        }
    }
}

