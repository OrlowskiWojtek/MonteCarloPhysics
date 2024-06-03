#include "utils.hpp"
#include "photon_difussion.hpp"

void task_2()
{
    for(int i = 0; i < 4; i++){
        PHOTON_DIFFUSION_2D photon;

        photon.nlayers = 3;
        photon.xmax = 0.2;
        photon.x_source = 0.1;
        photon.dx_source = 0.0;
        photon.x_detect = 0.15;
        photon.dx_detect = 0.01;
        photon.nx = 100;
        photon.ny = 100;
        photon.rx0 = 0.8;
        photon.ry0 = 0.6;
        
        // liczba wiązek fotonowych
        int N = 2e5;

        photon.write_all_paths = 0;
        photon.write_source_detection_paths = 0;
        
        double n1 = 1.3, n2 = 1.0, n3 = 1.0;
        double mu1_s = 10., mu2_s = 190., mu3_s = 90.; 
        
        switch(i){
            case 0: n2 = 1.5; break;
            case 1: n2 = 2.5; break;
            case 2: n1 = 1.0; n2 = 1.5; break;  
            case 3: n1 = 1.0; n2 = 1.5; mu2_s = 10.; break;
            default: break;
        }

        photon.fill_m_layer(1,1.,mu1_s,0.02,0.75,n1);
        photon.fill_m_layer(2,1.,mu2_s,0.02,0.075,n2);
        photon.fill_m_layer(3,10.,mu3_s,0.02,0.95,n3);

        photon.init();

        for(int k = 0; k < N; k++){
            photon.single_path();
        }

        photon.export_to_files("task2_"+to_string(i)+".dat");
    }    
}

void task_3()
{
  for(int i = 0; i < 4; i++){
        PHOTON_DIFFUSION_2D photon;

        photon.nlayers = 3;
        photon.xmax = 0.2;
        photon.x_source = 0.1;
        photon.dx_source = 0.0;
        photon.x_detect = 0.15;
        photon.dx_detect = 0.01;
        photon.nx = 100;
        photon.ny = 100;
        photon.rx0 = 0.0;
        photon.ry0 = 1.0;
        
        // liczba wiązek fotonowych
        int N = 2e5;

        photon.write_all_paths = 0;
        photon.write_source_detection_paths = 0;
        
        double n1 = 1.3, n2 = 1.0, n3 = 1.0;
        double mu1_s = 10., mu2_s = 190., mu3_s = 90.; 
        double mu1_a = 1., mu2_a = 1., mu3_a = 10.; 
        double g2 = 0.075; 
        switch(i){
            case 0: break;
            case 1: n1 = 1.0; mu2_a = 10.; mu2_s = 210.; n2 = 1.5; break;
            case 2: n1 = 1.0; mu2_a = 1. ; mu2_s = 210.; n2 = 1.5; break;  
            case 3: n1 = 1.0; mu2_a = 10.; mu2_s = 210.; n2 = 1.5; g2 = 0.75; break;
            default: break;
        }

        photon.fill_m_layer(1,mu1_a,mu1_s,0.02,0.75,n1);
        photon.fill_m_layer(2,mu2_a,mu2_s,0.02,g2,n2);
        photon.fill_m_layer(3,mu3_a,mu3_s,0.02,0.95,n3);

        photon.init();

        for(int k = 0; k < N; k++){
            photon.single_path();
        }

        photon.export_to_files("task3_"+to_string(i)+".dat");
    }    
}

