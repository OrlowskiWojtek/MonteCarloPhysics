#include <filesystem>
#include "dsmc.hpp"

using namespace std;


//=========================================================================================
// definicje funkcji
//=========================================================================================


/****************** INITIALIZATION ********************************************************
 * 
 * tworzymy tablice, generujemy rozklad poczatkowy, ustalamy brzegi pudla obliczeniowego
 * 
 */
void DSMC_2D::init(){
	
	generator.seed(seed);
	time_sum=0;
	
	
	nxy=nx*ny;
	delta_x=(xmax-xmin)/nx;
	delta_y=(ymax-ymin)/ny;
	
	ntot=0;
	for(int i=0;i<n_mix;i++) ntot+=nc[i];
		
	par.resize(ntot);
	indx0.resize(nx*ny,vector<int>(10));
	indx.resize(ntot);
	temp_cell.resize(nx,vector<double>(ny,temp));
	dens_cell.resize(nx,vector<double>(ny,0.0));
	vx_cell.resize(nx,vector<double>(ny));
	vy_cell.resize(nx,vector<double>(ny));
	velocity_x.resize(nx);
	velocity_y.resize(ny);
	press_x.resize(nx);
	press_y.resize(ny);
	press_cell.resize(nx,vector<double>(ny));
	press_tens.resize(nx,vector<vector<double>>(ny,vector<double>(4)));
	
/* ROZKLAD POCZATKOWY */
	init_distribution();
	
/* pierwsze sortowanie czstek w komorkach*/
	sort();
	
/*  wezly/krawedzie brzegu zewnetrznego */			
	  edge_out[0][0]=xmin;
	  edge_out[0][1]=ymin;
	  edge_out[1][0]=xmin;
	  edge_out[1][1]=ymax;
	  edge_out[2][0]=xmax;
	  edge_out[2][1]=ymax;
	  edge_out[3][0]=xmax;
	  edge_out[3][1]=ymin;
	  edge_out[4][0]=xmin;
	  edge_out[4][1]=ymin;
	  
	return;
}







/*
 ************************   EVOLUTION   ********************************************
 * 
 * dt=0.999*delta_min/vmax - czastki ze srodka nie dotra do brzegu zewnetrznego
 * 
 */

void DSMC_2D::evolution(double time_max,int nit){
	int i,it=0;
	time_t czas1,czas2;
	double czas; //czas trwania obliczen na komputerze [s]
	string nazwa,katalog;
	
	katalog="wyniki";
	if(!std::filesystem::exists(katalog)){
		std::filesystem::create_directory(katalog);
	}
	

	FILE *fp;
	fp=fopen("xy_1.dat","w");
		
	time(&czas1);
	printf("it\tdt\ttime\tczas_obl\ttemp_eff\n");
	while(time_sum<time_max || it<nit){
		
		step();
		it++;
		time(&czas2);
		czas=difftime(czas2,czas1);
		
		if((it%100) == 0){
			hist_velocity_all(("data/hist"+std::to_string(it)+".dat").c_str(),5.0,50);
		}
		
		//if((it%20)==1){	
		//	write_position_velocity(("data/data_evol/data_evolution_it"+to_string(it)).c_str());
		//}

		mean_free_path();
		autocorrelation();
		diffusion_coeff+=cv_coeff*dt/3;
		
		
		printf("%5d  %10.1E  %10.2E %10.1f %10.2f %10.2E %10.2E %10.2E %10.2E %10.2E\n",\
		it,dt,time_sum,czas,temperature_eff(),l_mfp_num,l_mfp_teo,cv_coeff,diffusion_coeff,press_x[0]);
		
		fprintf(fp," %15.7E   %15.7E  \n",par[0].x,par[0].y);
		fflush(fp);
		
	}
	

	fclose(fp);
	return;	
}





/*
 ************************   CELL/PARTICLES SORTING   ***********************************
 *  indx0[i+j*nx][0] - liczba wszystkich czastek w komorce (i,j)
 *  indx0[i+j*nx][1] - poczatek indeksow dla komorki (i,j) w tablicy indx[ntot]
 *  indx0[i+j*nx][2] - koniec indeksow dla komorki (i,j) w tablicy indx[ntot]
 *  indx0[i+j*nx][3] - aktualna liczba indeksow dla komorki (i,j) w tablicy indx[ntot]
 * 
 * **************************************************************************************
 */
inline void DSMC_2D::sort()
{
		
	int i,j,ix,iy,nsum,numer,m;
	double x,y;
	
	/*
	 * zerujemy liczbe czastek w komorkach
	 */
	 
	for(i=0;i<nx*ny;i++){
		for(j=0;j<4;j++){
			indx0[i][j]=0; 
		}
	}
	
	
	/*
	 * zabezpieczenie przed wyjsciem czastki poza bariery zewnetrzne
	 * gdyby ktoras wyszla poza zakres to przenosimy ja do srodka blisko brzegu
	 * 
	 */
	
	for(i=0;i<ntot;i++){
		x=par[i].x;
		y=par[i].y;
		if(x<xmin)par[i].x=xmin+delta_x*1.0E-3;
		if(x>xmax)par[i].x=xmax-delta_x*1.0E-3;
		if(y<ymin)par[i].y=ymin+delta_y*1.0E-3;
		if(y>ymax)par[i].y=ymax-delta_y*1.0E-3;
	}
	
	
	
	/*
	 * kasujemy stare indeksy czastek
	 */
	for(i=0;i<ntot;i++) indx[i]=-1;
	
	
	/*
	 * okreslamy przynaleznosc czastek do komorki i zliczamy czastki w komorce
	*/
	for(i=0;i<ntot;i++){
		//ix=(int) ((par[i].x-xmin)/delta_x);
		//iy=(int) ((par[i].y-ymin)/delta_y);
		
		cell_ix_iy(par[i].x,par[i].y,&ix,&iy);
		
		par[i].ix=ix;
		par[i].iy=iy;
		m=global_index_cell(ix,iy);
		indx0[m][0]++; //sumujemy 
	}
	
	
	/*
	 * okreslamy start i koniec indeksow w tablic indx[ntot] dla kazdej komorki
	 * jesli komorka jest pusta to "indx_start>indx_end" co zapobiega iteracji po tej komorce
	*/
	
	nsum=0;
	for(m=0;m<nx*ny;m++){
			indx0[m][1]=nsum;	   //start
			nsum+=indx0[m][0];
			indx0[m][2]=nsum-1; //end
	}
	

	
	/*
	 * wypelniamy tablice indeksow czastek dla kazdej komorki: indx[ntot]
	 */
	for(i=0;i<ntot;i++){
		ix=par[i].ix;
		iy=par[i].iy;
		if(ix>=0&&ix<nx && iy>=0&&iy<ny){
			m=global_index_cell(ix,iy);
			indx0[m][3]++; //okreslamy aktualna liczbe znalezionych czastek
			numer=indx0[m][1]+indx0[m][3]-1; //numer w tablicy indx[]
			indx[numer]=i; //numer czastki
		}
	}

	
	/*
	 * sprawdzamy czy tablica indeksow czastek jest pelna
	*/
	nsum=0;
	for(i=0;i<ntot;i++){
		if(indx[i]>=0)nsum++;
	}
		
	if(ntot!=nsum){
		cout<< "ERROR =>  ntot, nsum=" << ntot <<"\t"<<nsum<<endl;
	}
	return;
}


/*
  ************************   CELL: IX, IY*********************************************
  * okreslamy numer komork w ktorej jest czastka 
  * **********************************************************************************
*/
inline void DSMC_2D::cell_ix_iy(double x,double y, int *ix,int *iy)
{
	*ix=(int) ((x-xmin)/delta_x);
	*iy=(int) ((y-ymin)/delta_y);
	return;
}



/*
  ************************   GLOBAL CELL INDEX 2D  ***********************************
  * okreslamy globalny numer komorki
  * 
  * **********************************************************************************
*/
inline int DSMC_2D::global_index_cell(int ix,int iy)
{
	int m=ix+iy*nx;
	return m;
}



/*
  ************************   RUN 1 STEP **********************************************
  * icol=0,1:  1 - zderzenia czastka-czastka, 0-brak zderzen (czastki nie widza sie)
  * **********************************************************************************
*/
inline void DSMC_2D::step()
{
	
	int i,j,k,ix,iy,m,im,ichunk;

	/*****************************************************
	* szukamy dt
	*/	
		double vmax=0.;
		for(i=0;i<ntot;i++) vmax=max(vmax,par[i].v);
		dt=min(delta_x,delta_y)*0.999/vmax; 
	
	
	/*
	 * flaga: irun=0 - czastka nie wykonala ruchu w danym kroku
	 */	
		for(i=0;i<ntot;i++){
			par[i].irun=0;
		}
	
	/*
	 * zerujemy tablice tensora cisnienia - liczony na biezaco
	 * pxx-0, pxy-1, pyx-2, pyy-3
	 */
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				for(k=0;k<4;k++) press_tens[i][j][k]=0.;
			}
		}
	
		
	
	
	/*
	 * lecimy po komorkach -> czastkach w srodku komorki
	 * 
	 */
	
	
	omp_set_num_threads(nthreads);
    
	ichunk=3; //conajmniej 4 aby nie bylo kolizji watkow	
#pragma omp parallel default(shared) private(ix,iy,m,im,i)
 {
  #pragma omp for schedule(dynamic,ichunk) ordered
	
	 	for(iy=0;iy<ny;iy++){
			for(ix=0;ix<nx;ix++){

			m=global_index_cell(ix,iy);
			// po czastkach
			for(im=indx0[m][1];im<=indx0[m][2];im++){
				i=indx[im];
		
				if((icol==0 || icol==1) && par[i].irun==0 ){
					particle_collisions(i,icol);
				}
				else if( icol!=0 && icol!=1 ){
					cout<< "bledny parametr zderzenia icol= "<< icol<< endl;
				}	
			}
		}
	}
 }	
 
 /*************************************************************
  * sotujemy czastki
  * wyznaczamy: temperature, cisnienie i gestosc w komorkach
  * 
  *************************************************************
  */
	sort();
 	temperature_cell();
	compute_pressure();
	density_cell();
	velocity_cell();
 	time_sum+=dt;
	
	
	return;
}






/*
 ***************** COLLISIONS *********************************************************
 * 
 * obsluga zderzen
 * 
 * i - indeks 1 czastki
 * i2 - indeks 2 czastki (te wyszukujemy)
 * 
 * itype=0,1:   0-brak zderzen, 1-uwzgledniamy zderzenia
 * 
 * ***************************************************************************************
 */
inline int DSMC_2D::particle_collisions(int i,int itype_scat)
{

	int i2,ix,iy,ix2,iy2,im,im2,m,m2;
	double x1,y1,x2,y2,vx1,vy1,vx2,vy2,vx_rel,vy_rel;
	double x_rel,y_rel;
	double rc1,rc2;
	double a,b,c,delta,dt0,dt1,dt2,dtmin;
	int i2_min;
	double rr_eff,rr_rel,vv_rel,rv_rel;
	int icount_1,icount_2;
	
	PARTICLE p1,p2;
	

	/*
	 * 
	 * czastki nie widza sie -> brak zderzen -> zwobodny lot + odbicia od brzegow
	 * 
	 */	
	
	if(itype_scat==0 && par[i].irun==0){
		
		double dti=dt;
		int icount; //(0-brak zderzenia z bariera, 1-zderzenia)
		
		do{
			icount=barrier_reflection(par[i],&dti); 
		}while(icount==1);
		
		if(dti>0)	free_flight(par[i],&dti);
		par[i].irun=1;
		
		return 0;
	}
	
	
	
	
	/*
	 * 
	 * czastki zderzaja sie -> musimy wyszukac druga czastke z ktora najszybciej zderzy sie i-ta
	 * moze sie zdarzyc ze czastka byla juz uwzgledniona w zderzeniu z inna -> wtedy irun==1
	 * 
	 */
	
	if(itype_scat>0 && par[i].irun==0){
			
		dtmin=dt*1.0; //szukamy minimalnego kroku (<=dt)
		i2_min=-1; //blokada rozpraszania
		ix=par[i].ix;
		iy=par[i].iy;
		for(ix2=max(0,ix-1);ix2<=min(nx-1,ix+1);ix2++){
			for(iy2=max(0,iy-1);iy2<=min(ny-1,iy+1);iy2++){
				m2=ix2+iy2*nx;
				for(im2=indx0[m2][1];im2<=indx0[m2][2];im2++){
				i2=indx[im2];
				
				
				
			/*
			 *
			 * okreslamy parametr zderzenia czastek o indeksach (i,i2)
			 * 
			 */
		
			
				if(par[i2].irun==0 && i2!=i){
								
						
					
					
					x1=par[i].x;
					y1=par[i].y;
					vx1=par[i].vx;
					vy1=par[i].vy;
					rc1=par[i].rc;
							
					x2=par[i2].x;
					y2=par[i2].y;
					vx2=par[i2].vx;
					vy2=par[i2].vy;
					rc2=par[i2].rc;
      						
					vx_rel=vx1-vx2;
					vy_rel=vy1-vy2;
							
					x_rel=x1-x2;
					y_rel=y1-y2;

				
					
							
				/* 
				* sprawdzamy warunek zblizenia obu czastek -> rownanie kwadratowe na dt0
				* szukamy pierwiastkow rzeczywistych i dodatnich: min(dt1,dt2)>0
				* 
				* ponizsze kilka linijek najbardziej spowalnia program (okolo 40 razy)
				* 
				* 
				*/
				
					rv_rel=x_rel*vx_rel+y_rel*vy_rel;
					rr_rel=pow(x_rel,2)+pow(y_rel,2);
					vv_rel=pow(vx_rel,2)+pow(vy_rel,2);
					
			
					rr_eff=pow(rc1+rc2,2);
					
					a=1.0;	
					b=2.0*rv_rel/vv_rel;
					c=(rr_rel-rr_eff)/vv_rel;
					delta=b*b-4*a*c;
					
					if(delta>=0){
						dt1=(-b-sqrt(delta))/2/a;  //UWAGA: b moze byc ujemne 
						dt2=(-b+sqrt(delta))/2/a; 
						
						if(dt1>0 && dt2>0){
							dt0=min(dt1,dt2);
						}
						else if(dt1>0 && dt2<0){
							dt0=dt1;
						}
						else if(dt1<0 && dt2>0){
							dt0=dt2;
						}
						else{
							dt0=-1.0; //blokada
						}
							
						if(dt0>0 && dt0<dtmin){
							dtmin=dt0;
							i2_min=i2;					
						}
					}
				   }
	
				}
			}
		}
					
	
	
	
	
		/*
		 * 
		 * obslugujemy zderzenie dwoch czastek
		 * zderzaja sie czastki dla dtmin
		 * sprawdzamy czy wczesniej nie wpadna na bariere
		 * jak nie wychodza na brzeg to je rozpraszamy w ukladzie CM
		 *				 
		 */
		
		

		
		if(i2_min>=0 && i2_min<ntot){
						
			/*
			 * sprawdzamy czy wczesniej czastki nie wpadna na bariere (dtmin)
			 *
			 */
			
			p1=par[i];
			p2=par[i2_min];		
						
			dt1=dtmin;
			dt2=dtmin;
			icount_1=barrier_reflection(p1,&dt1); 
			icount_2=barrier_reflection(p2,&dt2); 
						
			
								
			/*
			 * rozpraszanie w ukladzie CM
			 * czastki nie wpadaja na bariere
			 * 
			 */
			 
			if(icount_1==0 && icount_2==0){
						    
			    dt1=dtmin;
			    dt2=dtmin;
			    free_flight(par[i],&dt1);
			    free_flight(par[i2_min],&dt2);		
			    
		 /* kopie starych predkosci obu czastek */   
			    p1=par[i];
			    p2=par[i2_min];
			    
			    particle_particle_scattering_balls(par[i],par[i2_min],itype_scat);			      
			    
		 /* liczymy wklad 2 wyrazu (oddzialywania) do tensora cisnienia */ 
			 
			    compute_pressure_tensor_interaction(p1,p2,par[i],par[i2_min]);
			    
			    dt1=dt-dtmin;   //tyle czasu pozostalo na ruch po zderzeniu        
			    dt2=dt-dtmin;      
			}
			else{
				dt1=dt;  //wczesniej zderzenie z bariera wiec restart kroku dt
				dt2=dt;				
			}
						
				
		/*
		 * sprawdzamy czy po rozproszeniu nie wpadna na bariere
		 * 
		 */                         
			if(dt1>0){    
					do{
						icount_1=barrier_reflection(par[i],&dt1); 
					}while(icount_1==1);
						
					do{
        				    icount_2=barrier_reflection(par[i2_min],&dt2);                         
		                  }while(icount_2==1);					                         
			}                         
							                         	
            /*
             * dla pozostalego czasu czastki poruszaja sie swobodnie
		 * 
             */
					    
			if(dt1>0)free_flight(par[i],&dt1);
			if(dt2>0)free_flight(par[i2_min],&dt2);		
						
		/*
		 * zmieniamy status czastek -> wykonaly ruch:  irun=1 
		 */		
			par[i].irun=1;
			par[i2_min].irun=1;
						
		}
		else if(i2_min<0)
		{
			dt1=dt;
			do{
				icount_1=barrier_reflection(par[i],&dt1); 
			}while(icount_1==1);
			
			if(dt1>0)free_flight(par[i],&dt1);
			par[i].irun=1;
		}
		
		
	}//itype
	
	return 0;
}




/*
 ***************** SCATTERING-HARD BALLS *************************************************
 * 
 * zderzenie dwoch czastek, pedy rzutowane sa na kierunek laczacy srodki kul
 * na tym kierunku odbicie kul jest problemem jednowymiarowym
 * przeksztalcenia wykonne w programie MAPLE
 * 
 * itype=0,1:   0-klasyczne (odbicie w CM), 1-odbicie ze stochastyczna zmiana kierunku
 * 
 * ***************************************************************************************
 */
inline void DSMC_2D::particle_particle_scattering_balls(PARTICLE & p1, PARTICLE & p2,int itype_scat)
{

	double x1,y1,x2,y2,vx1,vy1,vx2,vy2;
	double mc1,mc2;
	double ux,uy,u;
	double t5,t6,t7,t8,t9,t10,t11,t13,t14,t15,t17,t18,t22,t35,t29,t23,t25;
	
	x1=p1.x;
	y1=p1.y;
	vx1=p1.vx;
	vy1=p1.vy;
	mc1=p1.mc;
	
	x2=p2.x;
	y2=p2.y;
	vx2=p2.vx;
	vy2=p2.vy;
	mc2=p2.mc;
	
	
 /*
  * rozpraszanie klasyczne w ukladzie CM
  * 
  */   	
	
	
    if(itype_scat==2){
		
	//vx1
	t5 = mc1*vx1;
      t6 = -x1+x2;
      t7 = fabs(t6);
      t8 = t7*t7;
      t9 = y1-y2;
      t10 = fabs(t9);
      t11 = t10*t10;
      t13 = sqrt(t8+t11);
      t14 = 1/t13;
      t15 = -t6*t14;
      t18 = t9*t14;
      t22 = -t6*t14*(t15*t5+t18*mc1*vy1);
      t35 = (2.0*(t22-t6*t14*(t15*mc2*vx2+t18*mc2*vy2))/(mc1+mc2)*mc1-2.0*t22+t5)/mc1;
	p1.vx=t35;
	
	//vy1
	t6 = -x1+x2;
      t7 = fabs(t6);
      t8 = t7*t7;
      t9 = -y1+y2;
      t10 = fabs(t9);
      t11 = t10*t10;
      t13 = sqrt(t8+t11);
      t14 = 1/t13;
      t15 = -t6*t14;
      t17 = mc1*vy1;
      t18 = -t9*t14;
      t22 = -t9*t14*(t15*mc1*vx1+t18*t17);
      t35 = (2.0*(t22-t9*t14*(t15*mc2*vx2+t18*mc2*vy2))/(mc1+mc2)*mc1-2.0*t22+t17)/mc1;
	p1.vy=t35;
	
	//vx2
	t6 = -x1+x2;
      t7 = fabs(t6);
      t8 = t7*t7;
      t9 = -y1+y2;
      t10 = fabs(t9);
      t11 = t10*t10;
      t13 = sqrt(t8+t11);
      t14 = 1/t13;
      t15 = -t6*t14;
      t18 = -t9*t14;
      t23 = mc2*vx2;
      t29 = -t6*t14*(t15*t23+t18*mc2*vy2);
      t35 = (2.0*(-t6*t14*(t15*mc1*vx1+t18*mc1*vy1)+t29)/(mc1+mc2)*mc2-2.0*t29+t23)/mc2;
	p2.vx=t35;
	
	//vy2
	t6 = x1-x2;
      t7 = fabs(t6);
      t8 = t7*t7;
      t9 = -y1+y2;
      t10 = fabs(t9);
      t11 = t10*t10;
      t13 = sqrt(t8+t11);
      t14 = 1/t13;
      t15 = t6*t14;
      t18 = -t9*t14;
      t25 = mc2*vy2;
      t29 = -t9*t14*(t15*mc2*vx2+t18*t25);
      t35 = (2.0*(-t9*t14*(t15*mc1*vx1+t18*mc1*vy1)+t29)/(mc1+mc2)*mc2-2.0*t29+t25)/mc2;
	p2.vy=t35;
	
    }
 /*
  *
  * rozpraszanie z losowo obracanymi pedami w ukladzie CM
  * 
  */   
    else if(itype_scat==1)
    {
	
	double vcm_x,vcm_y,teta,vx10,vy10,vx11,vy11,vx21,vy21,ss,cc;    
	vcm_x=(mc1*vx1+mc2*vx2)/(mc1+mc2);
	vcm_y=(mc1*vy1+mc2*vy2)/(mc1+mc2);
	
	vx10=vx1-vcm_x;
	vy10=vy1-vcm_y;
	teta=rng_U(generator)*2*M_PI;
	cc=cos(teta);
	ss=sin(teta);
	vx11=cc*vx10-ss*vy10;
	vy11=ss*vx10+cc*vy10;
	vx21=-vx11*mc1/mc2;
	vy21=-vy11*mc1/mc2;
	
	p1.vx=vx11+vcm_x;
	p1.vy=vy11+vcm_y;
	p2.vx=vx21+vcm_x;
	p2.vy=vy21+vcm_y;
	    
    }
	
	
	
	
    /* liczymy moduly predkosci po rozproszeniu */
	
	p1.v=sqrt(pow(p1.vx,2)+pow(p1.vy,2));			      
	p2.v=sqrt(pow(p2.vx,2)+pow(p2.vy,2));
	
    /* zliczamy zderzenia z innymi czastkami */
	p1.ncol++;
	p2.ncol++;
	

	return;
}




/*
 ***************** PRESSURE  *************************************************************
 * 
 * liczymy cisnienie w kazdej komorce
 * 
 * 
 * ***************************************************************************************
 */
inline void DSMC_2D::compute_pressure()
{
	int ix,iy;
	
	compute_pressure_tensor_kinetic();
	
	for(ix=0;ix<nx;ix++){
		for(iy=0;iy<ny;iy++){
			press_cell[ix][iy]=(press_tens[ix][iy][0]+press_tens[ix][iy][3])/2/delta_x/delta_y;
		}
	}
	
	for(ix=0;ix<nx;ix++){
		press_x[ix]=0.;
		for(iy=0;iy<ny;iy++){
			press_x[ix]+=press_cell[ix][iy]/ny;
		}
	}
	
	for(iy=0;iy<ny;iy++){
		press_y[iy]=0.;
		for(ix=0;ix<nx;ix++){	
			press_y[iy]+=press_cell[ix][iy]/nx;
		}
	}
		
	
	return;
}




/*
 ***************** PRESSURE TENSOR: KINETIC **********************************************
 * 
 * liczymy wklad kinetyczny do tensora cisnienia - dla wszystkich komorek i czastek 
 * 
 * 
 * ***************************************************************************************
 */
inline void DSMC_2D::compute_pressure_tensor_kinetic()
{
	double vx_av,vy_av,dvx,dvy,mc;
	int i,ix,iy,im,m,k;
	
	for(ix=0;ix<nx;ix++){
		for(iy=0;iy<ny;iy++){
			m=global_index_cell(ix,iy);
			// predkosc srednia w komorce

			 vx_av=0.;
			 vy_av=0.;
			for(im=indx0[m][1];im<=indx0[m][2];im++){
				i=indx[im];
				vx_av+=par[i].vx/indx0[m][0];
				vy_av+=par[i].vy/indx0[m][0];
			}
			
			for(im=indx0[m][1];im<=indx0[m][2];im++){
				i=indx[im];
				dvx=par[i].vx-vx_av;
				dvy=par[i].vy-vy_av;
				mc=par[i].mc;
				press_tens[ix][iy][0]+=mc*dvx*dvx;
				press_tens[ix][iy][1]+=mc*dvx*dvy;
				press_tens[ix][iy][2]+=mc*dvy*dvx;
				press_tens[ix][iy][3]+=mc*dvy*dvy;
			}
		}
	}
	return;
}

	


/*
 ***************** PRESSURE TENSOR: INTERACTION *************************************************
 * 
 * liczymy wklad oddzialywania do tensora cisnienia
 * 
 * 
 * ***************************************************************************************
 */
inline void DSMC_2D::compute_pressure_tensor_interaction(PARTICLE & p1_old, PARTICLE & p2_old,\
							  PARTICLE & p1_new, PARTICLE & p2_new)
{
		int ix1,iy1,ix2,iy2;
		double x1,y1,x2,y2;
		double vx1,vy1,vx2,vy2;
		double mc1,mc2;
		double x12,y12,x21,y21;
		double px1,py1,px2,py2;
	
		x1=p1_old.x;
		y1=p1_old.y;
		vx1=p1_old.vx;
		vy1=p1_old.vy;
		mc1=p1_old.mc;
	
		x2=p2_old.x;
		y2=p2_old.y;
		vx2=p2_old.vx;
		vy2=p2_old.vy;
		mc2=p2_old.mc;
	
		cell_ix_iy(x1,y1,&ix1,&iy1);
		cell_ix_iy(x2,y2,&ix2,&iy2);
	
		x12=x1-x2;
		y12=y1-y2;
	
		x21=-x12;
		y21=-y12;
	
		px1=mc1*(p1_new.vx-vx1);
		py1=mc1*(p1_new.vy-vy1);
		px2=mc2*(p2_new.vx-vx2);
		py2=mc2*(p2_new.vy-vy2);
	
		press_tens[ix1][iy1][0]+=x12*px1/2./dt;
		press_tens[ix1][iy1][1]+=x12*py1/2./dt;
		press_tens[ix1][iy1][2]+=y12*px1/2./dt;
		press_tens[ix1][iy1][3]+=y12*py1/2./dt;
	
		press_tens[ix2][iy2][0]+=x21*px2/2./dt;
		press_tens[ix2][iy2][1]+=x21*py2/2./dt;
		press_tens[ix2][iy2][2]+=y21*px2/2./dt;
		press_tens[ix2][iy2][3]+=y21*py2/2./dt;
		
		return;

}






/*
 ***************** FREE-FLIGHT ***********************************************
 * 
 * obsluga swobodnego lotu
 * 
 * ***************************************************************************
 */
inline void DSMC_2D::free_flight(PARTICLE & par,double *dti)
{
		par.x=par.x+par.vx*(*dti);
		par.y=par.y+par.vy*(*dti);
		par.path=par.path+par.v*(*dti);
		par.vxt+=par.vx*(*dti);
		par.vyt+=par.vy*(*dti);
		par.cv+=(*dti)*(par.vx0*par.vx+par.vy0*par.vy);
		*dti=0.0;
}
	
	
			
/*
 ***************** BARRIER REFLECTION ***********************************************
 * 
 * obsluga swobodnego lotu z uwzglednieniem odbic od brzegow
 * 
 * UWAGA: brzeg wewnetrzny zorientowany zgodnie ze wskazowkami zegara
 * 
 * (u1x,u1y) ,(u2x,u2y) - poczatek i koniec krawedzi brzegu
 * (x1,y1) ,(x2,y2) - poczatek i koniec pierwotnej trajektorii czastki
 * (wx,wy) - punkt przeciecia trajektorii z brzegiem
 * (zx,zy) - nowy kierunek predkosci czastki (wersor) po odbiciu od brzegu
 * 
 * return=0/1: 0-brak przeciecia, 1-znaleziono przeciecie trajektoria-brzeg
 * **********************************************************************************
 */
inline int DSMC_2D::barrier_reflection(PARTICLE & par,double *dti)
{
		double x1,y1,x2,y2,vx,vy,v;
		int i,inter;
		double u1x,u1y,u2x,u2y;
		double wx,wy,zx,zy;
		double dti2,xx,yy,r;
		double cross_z;
		
		
		/* 
		 * sprawdzamy odbicia od brzegu zewnetrznego 
		 * jesli sie odbije (return 1) to sprawdzamy jeszcze raz: mozliwe drugie odbicie
		 * odbicia sprawdzamy tylko gdy czastka jest w komorce przy brzegu (ix=0,nx-1; iy=0,ny)
		 * jesli jest w dalszych komorkach to najwyzej dotzre do brzegu ale jeszcze sie nie odbije
		 * 
		 */
		
		if(par.ix==0 || par.ix==(nx-1) || par.iy==0 || par.iy==(ny-1)){
		
			for(i=0;i<nodes_out;i++){
			
				x1=par.x;
				y1=par.y;
				vx=par.vx;
				vy=par.vy;
		
				x2=x1+vx*(*dti);
				y2=y1+vy*(*dti);
			
				u1x=edge_out[i][0];
				u1y=edge_out[i][1];
				u2x=edge_out[i+1][0];
				u2y=edge_out[i+1][1];
			
				
			/*
			 * sprawdzamy czy czastka "pada na bariere" od srodka czy juz sie odbila
			 * liczymy z-owa skaladowa il. wektorowego 
			 * (r_cross_u)_z <0 - to odbicie od srodka - obliczamy parametry zderzenia
			 * (r_cross_u)_z >=0 to czastka juz jest na barierze - brak akcji
			 */	
				
				
				cross_z=vx*(u2y-u1y)-vy*(u2x-u1x);
			      if(cross_z>=0) continue;
				
			
			/*  inter==1:  jest odbicie: 
			 *  zmieniamy polozenie na punkt przeciecia (lokalizujemy czastke na brzegu), 
			 *  predkosc w kierunku odbicia,
			 *  czas redukujemy o czas przelotu do przeciecia
			 */
			
				intersection(u1x,u1y,u2x,u2y,x1,y1,x2,y2,&wx,&wy,&zx,&zy,&inter);	
				
				if(inter==1){
				
				/*
				 * wyznaczamy czas lotu do bariery i obslugujemy go funkcja free_flight()
				 * 
				 */	
					xx=wx-x1;
					yy=wy-y1;
					r=sqrt(xx*xx+yy*yy);
					par.path=par.path+r;
					v=par.v;
					dti2=r/v;       //czas do uderzenia w bariere					
					free_flight(par,&dti2); 
					*dti=*dti-dti2; //zmniejszamy czas swobodnego lotu
					
				/*
				 * tempi[i]<0: rozpraszanie adiabatyczne (tylko zmiana kierunku)
				 * tempi[i]>0: rozklad MB na brzegu  dla temp=tempi[i] (losowanie)
				 * 
				 */
					
					if(tempi[i]<=1.0E-10){
						v=sqrt(vx*vx+vy*vy);
						par.vx=zx*v;
						par.vy=zy*v;
					}
					else{
						
						double vnew,vnew2,ux,uy,uv,u_norm;
						double vx_parallel,vy_parallel;
						double vx_perpendicular,vy_perpendicular;
						
					
						
						
												
						ux=u2x-u1x;
						uy=u2y-u1y;
						u_norm=sqrt(ux*ux+uy*uy);
						ux=ux/u_norm;
						uy=uy/u_norm;
					
					
						maxwell_distribution(&vnew,&vnew2,tempi[i],par.mc);
						
						vnew=sqrt(vnew*vnew+vnew2*vnew2);			
						vx_perpendicular=fabs(vnew)*uy;
						vy_perpendicular=-fabs(vnew)*ux;
						
						maxwell_distribution(&vnew,&vnew2,tempi[i],par.mc);
						
						vx_parallel=vnew*ux;
						vy_parallel=vnew*uy;
						
						par.vx=vx_parallel+vx_perpendicular;
						par.vy=vy_parallel+vy_perpendicular;
						par.v=sqrt(pow(par.vx,2)+pow(par.vy,2));
						
					}
					
					par.nbound_col++; // zliczamy liczbe zderzen ze sciankami
					
					return 1;
				}
			}
		}
		
		
	/*
	 * brak odbicia: realizujemy swobodny lot
	*/
	
	return 0; 
}



/*
 ***************** INTERSECTION  *************************************
 * 
 * dane sa odcinki: (u1x,u1y)->(u2x,u2y),    (v1x,v1y)->(v2x,v2y)
 * wektory przesuwamy o wektor u1 i obracamy o kat (pi-fi) - ustawiamy wektor u2||y"
 * inter=0,1 :    0-brak przeciecia, 1-przeciecie
 * punkt przeciecia (wx,wy)
 * kierunek odbicia (obrot wektora v wokol wekotra u): (zx,zy) w ukl LAB
 * 
 */
inline int DSMC_2D::intersection(double u1x,double u1y,double u2x,double u2y,
					   double v1x,double v1y,double v2x,double v2y,
					   double *wx,double* wy,double *zx,double *zy,int *inter)
{
	
	double u,cc,ss;
	double ux,uy,px,py;
	double v1x_p,v1y_p,v2x_p,v2y_p;
	double a,b,a_num,a_denom,skalar;

	/* -------------jesli wektory nie sa rownolegle to sprawdzamy czy sie nie przecinaja -------- */
	ux=u2x-u1x;
	uy=u2y-u1y;
	u=sqrt(ux*ux+uy*uy);
	
	cc=ux/u; //kosinus
	ss=uy/u; //sinus
	
	//obrot o kat alfa=(pi/2-fi): Op(alfa)*(vec{v}-vec{r})
	v1x_p=ss*(v1x-u1x)-cc*(v1y-u1y);
	v1y_p=cc*(v1x-u1x)+ss*(v1y-u1y);
	v2x_p=ss*(v2x-u1x)-cc*(v2y-u1y);
	v2y_p=cc*(v2x-u1x)+ss*(v2y-u1y);
	
	a_num=(v2y_p-v1y_p);
	a_denom=(v2x_p-v1x_p);
	a=a_num/a_denom; //wsp nachylenia prostej
	b=v2y_p-a*v2x_p;//punkt przeciecia z osia 0y
	
	/*---- wektory sa rownolegle: EXIT -----*/
	if(fabs(a_num)>fabs(1.E+10*a_denom)){
		*inter=0;
		*wx=0.;
		*wy=0.;
		return 0;
	}
	else if(b>=(-u*1.E-12) && b<=u*(1.0+1.0E-12) && v1x_p*v2x_p<=0){
		
	   /* zabezpieczenie przed ominieciem wierzcholka */	
		if(b<0)b=0.0;
		if(b>u)b=u;
		
		*inter=1;
	   /*punkt przeciecia w ukl LAB*/
		*wx=cc*b+u1x;
		*wy=ss*b+u1y;
		
	   /*---------  kierunek odbicia wekt v wzgledem wekt u w ukl LAB -----------------*/
		ux=ss*(-1)*v1x_p+cc*v1y_p;  
		uy=-cc*(-1)*v1x_p+ss*v1y_p;
		px=ss*(-1)*v2x_p+cc*v2y_p;  
		py=-cc*(-1)*v2x_p+ss*v2y_p;
		
		ux=px-ux;
		uy=py-uy;
		
		u=sqrt(ux*ux+uy*uy);
		*zx=ux/u;
		*zy=uy/u;
		
	}
	else{
		*inter=0;		
		*wx=0.;
		*wy=0.;
	}
	
	
	
	return 0;
}






/*
 ***************** ZAPIS POLOZENIA i PREDKOSCI CZASTEK XY *************************************
 */
void DSMC_2D::write_position_velocity(const char* text)
{
	
	FILE *fp;
	fp=fopen(text,"w");
	for(int i=0;i<ntot;i++){
		fprintf(fp," %15.7E  %15.7E  %15.7E  %15.7E  %15.7E  %15.7E \n  ",\
		par[i].x,par[i].y,par[i].vx,par[i].vy,par[i].rc,par[i].mc);
	}
	fclose(fp);
	return;
}



/*
 ***************** ZAPIS brzegow komorek *************************************
 */
void DSMC_2D::write_cell_bounds(const char* text)
{
	FILE *fp=fopen(text,"w");
	for(int i=0;i<=nx;i++){
		fprintf(fp,"%12.3g  %12.3g\n",delta_x*i,ymin);
		fprintf(fp,"%12.3g  %12.3g\n\n\n",delta_x*i,ymax);
	}
	for(int i=0;i<=ny;i++){
		fprintf(fp,"%12.3g  %12.3g\n",xmin,delta_y*i);
		fprintf(fp,"%12.3g  %12.3g\n\n\n",xmax,delta_y*i);
	}
	fclose(fp);
	return;
}






/*
 ***************** HISTOGRAM *************************************
 * 
 * vmax=v_multi*sigma  (sigma=sqrt(kb*Temp/masa))
 * nhist - liczba binow w histogramie
 * 
 *****************************************************************
 */

void DSMC_2D::hist_velocity_all(const char* text,double v_multi,int nhist)
{
	double sigma,dist,vmax,dv,vx,vy,v;
	int i,j,ic;
	double **hist_v_num; //histogram numeryczny
	double **hist_v_teo;  //histogram teoretyczny		
	
	hist_v_num = new double*[n_mix];
	hist_v_teo = new double*[n_mix];
	for(int i = 0; i < n_mix; i++) {
		hist_v_num[i] = new double[nhist];
		hist_v_teo[i] = new double[nhist];
	}
	
	
	sigma=sqrt(kb*temp/mc[0]);
	vmax=v_multi*sigma;
	dv=vmax/nhist;
	
	//HISTOGRAM NUMERYCZNY DLA KAZDEGO TYPU CZASTEK
	for(ic=0;ic<n_mix;ic++){
		for(i=0;i<nhist;i++){
			hist_v_num[ic][i]=0.;
		}
	}
	
	
	for(i=0;i<ntot;i++){
		vx=par[i].vx;	
		vy=par[i].vy;
		v=sqrt(vx*vx+vy*vy);
		j=(int)(v/dv);
		ic=par[i].ic;
		if(j<nhist){
			hist_v_num[ic][j]+=1./nc[ic]/dv;
		}else{
			hist_v_num[ic][nhist-1]+=1./nc[ic]/dv; //tu ladujemy te ktore wychodza poza zakres
		}
	}
	
		
	//HISTOGRAM TEORETYCZNY DLA KAZDEGO TYPU CZASTEK
	for(ic=0;ic<n_mix;ic++){
		sigma=sqrt(kb*temp/mc[ic]);
		for(i=0;i<nhist;i++){
			v=dv*(i+0.5);
			dist=1/sigma/sigma*v*exp(-v*v/2/sigma/sigma);
			hist_v_teo[ic][i]=dist;
		}
	}
	
	
	FILE *fp=fopen(text,"w");
	
	for(i=0;i<nhist;i++){
		v=dv*(i+0.5);
		fprintf(fp,"%12.4E ",v);
		for(ic=0;ic<n_mix;ic++)	fprintf(fp,"%12.4E ",hist_v_num[ic][i]);
		for(ic=0;ic<n_mix;ic++)	fprintf(fp,"%12.4E ",hist_v_teo[ic][i]);
		fprintf(fp,"\n");
	}
	
	fclose(fp);
	
	
	
	/*---- zwalniamy pamiec ----------------------------------*/
	for(int i = 0; i < n_mix; i++)delete [] hist_v_num[i];
	delete [] hist_v_num;
	
	for(int i = 0; i < n_mix; i++)delete [] hist_v_teo[i];
	delete [] hist_v_teo;
	
	return;
}


/*
 ********************WRITE npTv x-DISTRIBUTIONS ********************************************
 * 
 * zapis do pliku co k-iteracji
 * 
 * p*V/T=nR - jak dla gazu doskonalego
 * 
 * 
 ******************************************************************************************* 
 */
inline void  DSMC_2D::write_nptv(const char * filename,int msr0)
{
	double temp_av,dens_av,p_av,v_av,vx_av,vy_av,jx_av,nR_av;
	int msr,ksr;
	FILE *fp;
	fp=fopen(filename,"w");
		
	msr=min(msr0,nx);
	
	for(int ii=0;ii<nx;ii+=msr){
		temp_av=0.;
		dens_av=0.;
		p_av=0.;
		vx_av=0.;
		vy_av=0.;
		jx_av=0.;
		nR_av=0.;
		for(int ix=ii;ix<min(ii+msr,nx);ix++){
			ksr=min(ii+msr,nx)-ii;
			for(int iy=0;iy<ny;iy++){
				temp_av+=temp_cell[ix][iy]/ksr/ny;
				dens_av+=dens_cell[ix][iy]/ksr/ny;
				p_av+=press_cell[ix][iy]/ksr/ny;
				vx_av+=vx_cell[ix][iy]/ksr/ny;
				vy_av+=vy_cell[ix][iy]/ksr/ny;
				jx_av+=vx_cell[ix][iy]*dens_cell[ix][iy]/ksr/ny;
				
			}
		}
		v_av=sqrt(pow(vx_av,2)+pow(vy_av,2)); 
		nR_av=p_av*delta_x*delta_y/temp_av;
		fprintf(fp,"  %12.5E",(ii+0.5*msr)*delta_x);
		fprintf(fp,"  %12.5E",dens_av);
  	      fprintf(fp,"  %12.5E",p_av);
  	      fprintf(fp,"  %12.5E",temp_av);
  	      fprintf(fp,"  %12.5E",v_av);
		fprintf(fp,"  %12.5E",jx_av);
		fprintf(fp,"  %12.5E",nR_av);
		fprintf(fp,"\n");		    
		}
		fclose(fp);
	return;
}




/*
 ******************** AUTOCORRELATION C(t)**************************************************
 * 
 * wspolczynnik autokorelacji usredniony po czasie i po czastkach ktore nie uderzyly w scianki zewnetrzne
 * (z czasem bedzie ich coraz mniej)
 ******************************************************************************************* 
 */
inline void DSMC_2D::autocorrelation()
{
	int i,k;
	cv_coeff=0.;
	
	k=0;
	for(i=0;i<ntot;i++){
		if(par[i].nbound_col==0){
			cv_coeff+=par[i].cv;
			k++;	
		}
	}
	if(k>0) cv_coeff=cv_coeff/k/time_sum; 
	return;
}



/*
 ******************** MEAN FREE PATH *******************************************************
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::mean_free_path()
{
	
	double sigma,ro2d,pcf,gamma,v_rel_2d,ro2d_reduced;
	int k=0;
	
	l_mfp_num=0.;
	l_mfp_teo=0.;
	l_mfp_enskog=0.;
	
	for(int i=0;i<ntot;i++){
		if(par[i].ncol>0){	
			l_mfp_num+=par[i].path/par[i].ncol;
			k++;
		}
	}
	
	
	
	l_mfp_num=l_mfp_num/k;

	sigma=4*par[0].rc;
	ro2d=ntot/(xmax-xmin)/(ymax-ymin);  // =liczba czastek/objetosc
	l_mfp_teo=1/sigma/ro2d/sqrt(2.0);   
	
	
/*
 *   we wzorze na czestosc rozpraszania Enksoga jest blad - nie zgadzaja sie jednostki
 *
 *	
	ro2d_reduced=ntot*M_PI*pow(par[0].rc,2)/(xmax-xmin)/(ymax-ymin);
	pcf=(1.-7./16.*ro2d_reduced)/pow(1.-ro2d_reduced,2);
	gamma=4.*ro2d_reduced*ntot*pcf/(2.*par[0].rc*sqrt(kb*temp/M_PI/par[0].mc));
	v_rel_2d=sqrt(M_PI*kb*temp/par[0].mc);
	l_mfp_enskog=v_rel_2d/gamma;
*/
	return;
}




/*
 ******************** MAXWELL-BOLTZMANN 2D DISTRIBUTION **********************************
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::maxwell_distribution(double *vx, double *vy, double temp_i,double mci)
{
	double sigma=sqrt(kb*temp_i/mci);
	*vx=rng_N(generator)*sigma;
	*vy=rng_N(generator)*sigma;
	return;
}





/*
 ******************** EFFECTIVE TEMPERATURE **********************************
 * 
 ******************************************************************************************* 
 */
inline double DSMC_2D::temperature_eff()
{
	double temp1,ekin;
	ekin=0.;
	for(int i=0;i<ntot;i++){
		ekin+=par[i].mc*pow(par[i].v,2)/2.;
	}
	temp1=ekin/ntot/kb;
	return temp1;
}


/*
 ******************** EFFECTIVE TEMPERATURE IN CELLS  **************************************
 * 
 * okreslamy temperature efektywna na brzegach zewnetrznych tj. w komorkach przyleglych
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::temperature_cell()
{
	for(int ix=0;ix<nx;ix++){
		for(int iy=0;iy<ny;iy++){
			temp_cell[ix][iy]=0.;
			int m=global_index_cell(ix,iy);
			for(int im=indx0[m][1];im<=indx0[m][2];im++){
				int i=indx[im];
				temp_cell[ix][iy]+=par[i].mc*pow(par[i].v,2)/2;
			}
			if(indx0[m][0]>0)temp_cell[ix][iy]=temp_cell[ix][iy]/kb/indx0[m][0];
			
		}
	}
	return;
}


/*
 ******************** DENSITY IN CELLS  ****************************************************
 * 
 * okreslamy gestosc czastek w komorkach
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::density_cell()
{
	for(int ix=0;ix<nx;ix++){
		for(int iy=0;iy<ny;iy++){
			int m=global_index_cell(ix,iy);
			dens_cell[ix][iy]=indx0[m][0]/delta_x/delta_y;
		}
	}
	return;
}


/*
 ******************** VELOCITIES IN CELLS  ****************************************************
 * 
 * okreslamy srednie predkosci w komorkach
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::velocity_cell()
{
	for(int ix=0;ix<nx;ix++){
		for(int iy=0;iy<ny;iy++){
			vx_cell[ix][iy]=0.;
			vy_cell[ix][iy]=0.;
			int m=global_index_cell(ix,iy);
			int ile=indx0[m][0];
			for(int im=indx0[m][1];im<=indx0[m][2];im++){
				int i=indx[im];
				vx_cell[ix][iy]+=par[i].vx/ile;
				vy_cell[ix][iy]+=par[i].vy/ile;
			}
		}
	}
	return;
}



/*
 ******************** TEST  ****************************************************
 * 
 * procedura testowa do sprawdzania roznych rzeczy
 * 
 ******************************************************************************************* 
 */
inline void DSMC_2D::test()
{
	int i,j,k,n;
	double mci,temp_i,vx,vy,etot;
	n=10000;
	mci=par[0].mc;
	temp_i=tempi[0];
	etot=0.;
	for(i=0;i<n;i++){
		maxwell_distribution(&vx,&vy,temp_i,mci);	
		etot+=mci*(vx*vx+vy*vy)/2;
	}
	etot=etot/n/kb;
	//printf("%10.4E\n",etot)	;
	return;
}



/*
 ******************** INITIAL SPACE/VELOCITY DISTRIBUTION **********************************
 * 
 * init_dist=0   - wczytujemy z pliku
 * init_dist=1   - rozklad jednorodny predkosci [v=sqrt(2*kB*Temp/masa) ]
 * init_dist=2   - rozklad Maxwella-Boltzmanna  [z rozkladu normalnego]
 * 
 ******************************************************************************************* 
 */
void DSMC_2D::init_distribution()
{
	
	if(init_dist==0){
		double x,y,vx,vy,rci,mci;
		int k,ieof;
		
		printf("===  wczytujemy polozenia z pliku   ====\n\n");
		
		FILE *fp_in;
		fp_in=fopen("pos_vel_start.dat","r");
		
		if (fp_in == NULL) 
            {   
              printf("Brak pliku z polozeniami czastek -> KONIEC\n"); 
              exit(-1); 
            } 
            else {
		  k=0;
			for(int ic=0;ic<n_mix;ic++){
				for(int i=0;i<nc[ic];i++){
				
					ieof=fscanf(fp_in, "%lf %lf %lf  %lf %lf  %lf", &x, &y, &vx, &vy, &rci, &mci); 	
					if(ieof!=EOF){
						par[k].x=x;	
						par[k].y=y;	
						par[k].vx=vx;	
						par[k].vy=vy;
						par[k].x0=x;	
						par[k].y0=y;	
						par[k].vx0=vx;	
						par[k].vy0=vy;
						par[k].v=sqrt(vx*vx+vy*vy);	
						par[k].rc=rci;	
						par[k].mc=mci;	
						par[k].ic=ic;
						par[k].irun=0;
						par[k].ix=-1;
						par[k].iy=-1;
						k++;
					}
				}
			}
			
			if(k!=ntot){
				printf("za malo danych: ntot= %d, wczytano k=%d\n",ntot,k);
				printf("EXIT\n");
				exit(0);
			}
		}
	}	
	else{ 
		// (1) v=const, kierunek izotropowy;   (2) Maxwell-Boltzmann 2D
		
		int k=0;
		double x,y,eki,v,vx,vy,vx0,vy0,fi,sigma;
		
		
		k=0;
		for(int ic=0;ic<n_mix;ic++){
			for(int i=0;i<nc[ic];i++){
				if(init_dist==1||init_dist==2){
					x=rng_U(generator)*(xmax-xmin)+xmin;
					y=rng_U(generator)*(ymax-ymin)+ymin;
				}else if(init_dist==3||init_dist==4){
					x=rng_U(generator)*delta_x+xmin;
					y=rng_U(generator)*delta_y+ymin;
				}	
				   //MAXWELL-BOLTZMANN 
				sigma=sqrt(kb*temp/mc[ic]);
				vx=rng_N(generator)*sigma;
				vy=rng_N(generator)*sigma;
				  //IZOTROPOWY-SFERYCZNIE KONTUROWANY (v=const)	
				if(init_dist==1 || init_dist==3){
					v=sqrt(2.0)*sigma;
					vx0=v*vx/sqrt(vx*vx+vy*vy);
					vy0=v*vy/sqrt(vx*vx+vy*vy);
					vx=vx0;
					vy=vy0;
				}
				
				
			par[k].x=x;	
			par[k].y=y;	
			par[k].vx=vx;	
			par[k].vy=vy;
			par[k].x0=x;	
			par[k].y0=y;	
			par[k].vx0=vx;	
			par[k].vy0=vy;
			par[k].v=sqrt(vx*vx+vy*vy);	
			par[k].rc=rc[ic];	
			par[k].mc=mc[ic];	
			par[k].ic=ic;
			par[k].irun=0;
			par[k].ix=-1;
			par[k].iy=-1;
			par[k].indx=k; //zachowujemy indeks globalny - aby kopia obiektu wiedziala skad jest 
			par[k].ncol=0;
			k++;
			}
		}
		
		/*
		 * niewielka RENORMALIZACJA predkosci dla rozkladu MB do wlasciwej temperatury 
		 * (zazwyczaj jest zawyzona ze wzgledu na losowy charakter)
		 * dla ekin=const temperatura jest dobrze okreslona
		 * 
		 */
			double temp1=temperature_eff();
			for(int i=0;i<ntot;i++){
				vx=par[i].vx*sqrt(temp/temp1);
				vy=par[i].vy*sqrt(temp/temp1);
				v=sqrt(vx*vx+vy*vy);
				par[i].vx=vx;
				par[i].vy=vy;
				par[i].v=v;
			}
			
	}	
	return;
}









/*
 * ************************* WCZYTUJEMY DANE Z PLIKU *************************************************
 * 
 * ***************************************************************************************************
 */

void DSMC_2D::read(const char * text)
{
	FILE *fp_in;
	fp_in=fopen(text,"r");
	
          if (fp_in == NULL) 
            {   
              printf("Brak pliku -> KONIEC\n"); 
              exit(-1); 
            } else {
		  fscanf(fp_in, "%lf %lf %lf  %lf", &xmin, &xmax, &ymin, &ymax); 	
		  fscanf(fp_in, "%d %d", &nx, &ny); 		
		  fscanf(fp_in, "%lf ", &kb);
		  fscanf(fp_in, "%lf %lf %lf %lf %lf ", &temp, &tempi[0], &tempi[1], &tempi[2], &tempi[3]); 			
		  fscanf(fp_in, "%d ", &init_dist); 
		  fscanf(fp_in, "%d ", &n_mix); 
		  for(int i=0;i<n_mix;i++) 
			  fscanf(fp_in, "%d  %lf %lf ", &nc[i],&mc[i], &rc[i]);

		  fscanf(fp_in, "%d ", &nodes); 
		  for(int i=0;i<nodes;i++) 
			  fscanf(fp_in, "%lf %lf ", &edge[i][0],&edge[i][1]);
			// kopiujemy pierwszy wezel jako pierwszy - zeby utworzyc pare z dwoch kolejnych komorek
				edge[nodes][0]=edge[0][0];
				edge[nodes][1]=edge[0][1];
		  		
		  
		  printf("xmin, xmax, ymin, ymax=  %10.2g  %10.2g %10.2g %10.2g \n ",xmin,xmax,ymin,ymax);	
		  printf("nx,ny= %d %d\n",nx,ny);
		  printf("kb= %10.3E \n",kb);
		  printf("Temp, tempi[0-3]= %10.2f    %10.2f %10.2f %10.2f %10.2f\n",temp,tempi[0],tempi[1],tempi[2],tempi[3]);
		  printf("init_dist= %d\n",init_dist);
		  printf("n_mix= %d\n",n_mix);
		  for(int i=0;i<n_mix;i++) printf("nc, mc, rc = %d  %10.2g %10.2g \n", nc[i],mc[i], rc[i]); 			
		  
		  printf("nodes= %d\n",nodes);
		  if(nodes>0) for(int i=0;i<nodes;i++) printf("(x%d,y%d)=  %10.2g %10.2g \n",i,i, edge[i][0],edge[i][1]); 			
		  
		}
	fclose(fp_in);
	return;
}
