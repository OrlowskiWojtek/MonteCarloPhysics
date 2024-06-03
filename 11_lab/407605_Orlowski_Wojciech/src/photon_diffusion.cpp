#include<cmath>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<fstream>

using namespace std;

#include "photon_difussion.hpp"

/*
 * konstruktor z typowymi wartosciami
 * 
 */
PHOTON_DIFFUSION_2D::PHOTON_DIFFUSION_2D(){
		write_all_paths=0;
		write_source_detection_paths=0;
		nlayers=0;
		nx=0;
		ny=0;
		xmax=0.;
		ymax=0.;
		dx=0;
		dy=0;
		p_min=0.1;
		w_min=1.0E-4;
		abs_specular=0.; //akumulator dla odbic na WEJSCIU
		
		rx0=0.;
		ry0=1.;
		
	//layers=[mu_abs, mu_scat, width, g_anizo, n_refraction,y1,y2]
		int M=20;
		layers_data.resize(M,vector<double>(M,0.));		
			for(int i=0;i<M;i++){
				layers_data[i][3]=0.0; //g_anizo
				layers_data[i][4]=1.0; //n_refraction
			}
}		



/***********************************************************************************
 * 
 *	 			INITIALIZATION
 * 
 ***********************************************************************************/

void PHOTON_DIFFUSION_2D::init(){
	
	absorption.resize(nx+1,vector<double>(ny+1,0));
	reflectance.resize(nx+1,0.);
	transmittance.resize(nx+1,0.);
	
	dx=xmax/nx;
	ymax=0.;
	for(int i=1;i<=nlayers;i++)ymax+=layers_data[i][2];
	dy=ymax/ny;
	
	// okreslamy granice warstw;	
		for(int i=1;i<=nlayers+1;i++){
				layers_data[i][5]=layers_data[i-1][6];  //poczatek warstwy
				layers_data[i][6]+=layers_data[i][5]+layers_data[i][2];   //koniec warstwy
		}
	
	
	
	/*
	 * warunki
	 */
	if( (x_source-dx_source/2) <0 || (x_source+dx_source/2)>xmax){
		printf("source beyond region - STOP\n");
		printf("%15.5E   %15.5E   \n",x_source-dx_source/2,x_source+dx_source);
		exit(0);
	}
	printf("*******************************************************************************\n");
	printf("xmax = %15.7f \n",xmax);
	printf("ymax = %15.7f \n",ymax);
	printf("dx = %15.7f \n",dx);
	printf("dy = %15.7f \n",dy);
	printf("nx = %10d \n",nx);
	printf("ny = %10d \n",ny);
	printf("x_source = %15.7f \n",x_source);
	printf("dx_source = %15.7f \n",dx_source);	
	printf("x_detect = %15.7f \n",x_detect);
	printf("dx_detect = %15.7f \n",dx_detect);	
	
	
	
	printf("layer   wsp_absorption    wsp_scattering     thickness\t     g_aniz\t  n_refraction\t     y_bottom\t     y_top  \n");
	for(int i=1;i<=nlayers;i++){
		printf("%5d  ",i);
		for(int j=0;j<7;j++)printf("%15.7f ",layers_data[i][j]);
		printf("\n");
	}
	
	printf("*******************************************************************************\n");
	
	
	// kasujemy zawartosc plikow ze sciezkami
	FILE *fp;
	fp=fopen("source_detection_paths.dat","w");
	fclose(fp);
	fp=fopen("all_paths.dat","w");
	fclose(fp);
	
	
	
	return;
}


/**********************************
 * 	UNIFORM
 **********************************/
double PHOTON_DIFFUSION_2D::uniform(){
	return (double)rand()/RAND_MAX;
}


/**********************************************************************
 * 
 * 					SINGLE_PATH
 * 
 * 
 *  symulacja pojedynczej sciezki
 * 
 *  w: waga paczki fotonow
 * 
 **********************************************************************/
void PHOTON_DIFFUSION_2D::single_path(){
	
	
	// start paczki - inicjalizacja polozenia
	if( dx_source>=0.  && x_source >=0 && x_source<=xmax){
		beam.x=x_source+dx_source/2*(2*uniform()-1);
		beam.y=1.0E-10;
		beam.rx=rx0;
		beam.ry=ry0;
		beam.w=1.0;
		beam.alive=true;
		beam.layer=1;
		
		beam.path.clear();
		beam.path.push_back(beam.x);
		beam.path.push_back(beam.y);
		beam.path.push_back(beam.w);
		
		//pierwsze odbicie na interfejsie proznia/layer_1
		//double n0=layers_data[0][4];
		//double n1=layers_data[1][4];
		//double R_specular=pow(n0-n1,2)/pow(n0+n1,2);
		//beam.w=beam.w-R_specular;
		//abs_specular=abs_specular+R_specular;
	}

	int l=0;
	while(beam.alive==true){
		calculate_new_position(); //x_new,y_new
		scatter_in_layer(); 
		scatter_up_down_boundary();
		roulette();
		// zachowujemy sciezke 
		beam.path.push_back(beam.x);
		beam.path.push_back(beam.y);
		beam.path.push_back(beam.w);
		l++;
		//printf("%10d   %15.5E  %15.5E \n",l,beam.x,beam.y);
		
	}
	
		write_paths_to_file();
		
	
	return;
}

/********************************************************************************
 * 
 *  			WRITE_SOURCE_DETECTION_PATHS_TO_FILE
 * 
 ********************************************************************************/

void PHOTON_DIFFUSION_2D:: write_paths_to_file(){
	
	
	
	
	if(write_source_detection_paths==1 && fabs(beam.x-x_detect)<= dx_detect/2 ){
		FILE *fp;
		fp=fopen("source_detection_paths.dat","a");
		fprintf(fp,"\n");
		
		for(decltype(beam.path.size()) i=0;i<beam.path.size();i+=3){
			double x=beam.path[i];
			double y=beam.path[i+1];
			double w=beam.path[i+2];
			fprintf(fp,"%15.5E  %15.5E   %15.5E \n",x,y,w);
		}
		fclose(fp);
	}
	
	
	
	if(write_all_paths==1){
		FILE *fp;
		fp=fopen("all_paths.dat","a");
		fprintf(fp,"\n");
		
		for(decltype(beam.path.size()) i=0;i<beam.path.size();i+=3){
			double x=beam.path[i];
			double y=beam.path[i+1];
			double w=beam.path[i+2];
			fprintf(fp,"%15.5E  %15.5E   %15.5E \n",x,y,w);
		}
		fclose(fp);
	}
	
	
	
	return;
}






/********************************************************************************
 * 
 */

double PHOTON_DIFFUSION_2D:: sign(double x){
	if(x<0)return -1.0;
	else return 1.0;
}




/********************************************************************************
 *
 * 		  SEGMENT INTERSECTION
 * 
 * icross=0,1:   0-brak przeciecia, 1-przeciecie (0<alfa<1  oraz 0<beta<1)
 * 
 ********************************************************************************/

void PHOTON_DIFFUSION_2D:: segment_intersection(double x1,double y1, double x2,double y2,double x3,double y3, double x4,double y4,
								 double & x_cross,double & y_cross, int & icross){

	double a,b,c,d,b1,b2,alfa,beta,det;
	a=x2-x1;
	b=-(x4-x3);
	c=y2-y1;
	d=-(y4-y3);
	b1=x3-x1;
	b2=y3-y1;
	
	det=a*d-b*c;
	icross=0;
	if(fabs(det)>1.0E-50){
		alfa=(d*b1-b*b2)/det;
		beta=(-c*b1+a*b2)/det;
		if(alfa>=0 && alfa <=1.0 && beta>=0 && beta <=1.0){
			x_cross=x1+alfa*(x2-x1);
			y_cross=y1+alfa*(y2-y1);
			icross=1;
		}
	}
	return;
}





/********************************************************************************
 * 
 * 				SCATTER_UP_DOWN_BOUNDARY
 * 
 *  sprawdzamy czy mozemy wykonac przesuniecie: length=1
 *  jesli wiazka dochodzi (x_new,y_new) do gornego/dolnego brzegu warstwy 
 *  to przesuwamy ja na intgerfejs i zmieniamy kierunek
 * 
 *			 length=0
 *  
 ********************************************************************************/
void PHOTON_DIFFUSION_2D:: scatter_up_down_boundary(){
	
	if(beam.length==1){
		
		double y_bottom=layers_data[beam.layer][5];
		double y_top=layers_data[beam.layer][6];
		double x_cross,y_cross;
		int i,j;
		
		double x1,y1,x2,y2,x3,y3,x4,y4;
		int icross,ktory;
		
	/*************************************************
	 * sprawdzamy ktory brzeg przecina wiazka
	 * ktory=0 brak przeciecia
	 * ktory=1 LEFT
	 * ktory=2 RIGHT
	 * ktory=3 TOP
	 * ktory=4 BOTTOM
	 *************************************************/ 
			x1=beam.x;
			y1=beam.y;
			x2=beam.x_new;
			y2=beam.y_new;		
			ktory=0; //brak przeciecia+blokada dla naroznikow
		
		//lewy brzeg 
		if(ktory==0){
			x3=0;
			y3=y_bottom;
			x4=0;
			y4=y_top;
			segment_intersection(x1,y1,x2,y2,x3,y3,x4,y4,x_cross,y_cross,icross);	
			if(icross==1)ktory=1;
		}
		//prawy brzeg 
		if(ktory==0){
			x3=xmax;
			y3=y_bottom;
			x4=xmax;
			y4=y_top;
			segment_intersection(x1,y1,x2,y2,x3,y3,x4,y4,x_cross,y_cross,icross);				
			if(icross==1)ktory=2;
		}
		//gorny brzeg 
		if(ktory==0){
			x3=0;
			y3=y_top;
			x4=xmax;
			y4=y_top;
			segment_intersection(x1,y1,x2,y2,x3,y3,x4,y4,x_cross,y_cross,icross);				
			if(icross==1)ktory=3;
		}
		//dolny brzeg 
		if(ktory==0){
			x3=0;
			y3=y_bottom;
			x4=xmax;
			y4=y_bottom;
			segment_intersection(x1,y1,x2,y2,x3,y3,x4,y4,x_cross,y_cross,icross);				
			if(icross==1)ktory=4;
		}
		
		i=round(x_cross/dx);
		j=round(y_cross/dy);
		
		
	//************** LEFT-RIGHT BOUNDARY - ABSORPTION   *******************************
		if( ktory==1 || ktory==2){
				absorption[i][j]+=beam.w;	
				beam.x=x_cross;
				beam.y=y_cross;
				beam.length=0;
				beam.alive=false;
				return;
		}
				
			
	//************* TOP-BOTTOM BOUNDARY - REFLECTION/TRANSMISSION   **************************
		if( ktory==3 || ktory==4 ){
			
			//check total reflection
					
				double no,nn,sin_alfa_crit,sin_alfa;
				
				no=layers_data[beam.layer][4];
				if(ktory==3)nn=layers_data[beam.layer+1][4];
				if(ktory==4)nn=layers_data[beam.layer-1][4];
				
				sin_alfa_crit=nn/no;
				sin_alfa=fabs(beam.rx);
				
				if(sin_alfa>sin_alfa_crit){ //total reflection
					beam.ry=-beam.ry;
					beam.x=x_cross;
					beam.y=y_cross+beam.ry*1.0E-10; //lekko przesuwamy w dol aby nie bylo przeciecia
					beam.length=0;
					return;
				}else{  // reflection or transmission to upper layer	
					double alfa_o=fabs(asin(beam.rx));
					double alfa_n=fabs(asin(no/nn*beam.rx));
					double reflection=(pow(sin(alfa_o-alfa_n),2)/pow(sin(alfa_o+alfa_n),2) + pow(tan(alfa_o-alfa_n),2)/pow(tan(alfa_o+alfa_n),2)  )/2.;
					double u1=uniform();
					if(u1<reflection){ //reflection
						beam.ry=-beam.ry;
						beam.x=x_cross;
						beam.y=y_cross+beam.ry*1.0E-10; //lekko przesuwamy aby nie bylo przeciecia w kolejym kroku
						beam.length=0;
						return;
					}else{ //transmission 
						double rx=beam.rx; //stara wartosc: sin_alfa_old
						double ry=beam.ry; //stara wartosc: sin_alfa_old
						beam.rx=no/nn*rx; //kierunek uwzgledniony; sin_alfa_new
						beam.ry=sign(ry)*sqrt(1-pow(no/nn*rx,2)); //kierunek 
						beam.x=x_cross;
						beam.y=y_cross+beam.ry*1.0E-10; //lekko przesuwamy gora/dol aby nie bylo przeciecia w kolejnym kroku
						beam.length=0.;
						
						if(ktory==3)beam.layer++;
						else if(ktory==4) beam.layer--;
						
						if(beam.layer==(nlayers+1)){
							transmittance[i]+=beam.w;
							beam.w=0.;
							beam.alive=false;
						}
						else if(beam.layer==0){
							reflectance[i]+=beam.w;
							beam.w=0.;
							beam.alive=false;
						}
						return;
					}
				}
		}			
	}
	return;
}








/********************************************************************************
 * 
 * 				SCATTER_IN_LAYER
 * 
 *  sprawdzamy czy wiazka rozprasza sie w danej warstwie 
 *  length=1 i rozpraszanie zachodzi w warstwie to ja rozpraszamy 
 *   i ustawiamy:
 *			 length=0
 *  
 ********************************************************************************/
void PHOTON_DIFFUSION_2D:: scatter_in_layer(){
	double y1=layers_data[beam.layer][5];
	double y2=layers_data[beam.layer][6];
	
	if(beam.x_new>0 && beam.x_new<xmax && beam.y_new>y1 && beam.y_new<y2 && beam.length==1){
		
		double g,u1,u2,cos_teta,sign,sin_teta,rx,ry;
		u1=uniform();
		g=layers_data[beam.layer][3];
		if(g>1.0E-3){
			cos_teta=(1+g*g-pow((1-g*g)/(1-g+2.*g*u1),2))/2./g;
		}else{
			cos_teta=2.*u1-1.;
		}
		u2=uniform();
		if(u2<=0.5){
			sign=-1.0;
		}else{
			sign=1.0;
		}
		
	//rozpraszanie zachodzi zawsze - obrot wersora o kat (sign*teta)
		sin_teta=sign*sqrt(1-cos_teta*cos_teta);
		rx=beam.rx;
		ry=beam.ry;
		beam.rx=cos_teta*rx-sin_teta*ry;
		beam.ry=sin_teta*rx+cos_teta*ry;
	//absorpcja zachodzi zawsze - kasujemy czesc wiazki
		double dw=0.;
		double wsp_abs=layers_data[beam.layer][0];
		double wsp_scat=layers_data[beam.layer][1];
		dw=wsp_abs/(wsp_abs+wsp_scat)*beam.w;
		beam.w=beam.w-dw;
		int i,j;
		i=round(beam.x_new/dx);
		j=round(beam.y_new/dy);
		absorption[i][j]+=dw;
		beam.x=beam.x_new;
		beam.y=beam.y_new;
	//kasujemy mozliwosc przesuniecia w danej iteracji	
		beam.length=0;
	}
	return;
}











/*******************************************************************************
 * 
 *		 		CALCULATE_NEW_POSITION
 *  
 *  szukamy nowego proponowanego polozenia
 *  ustawiamy: beam.length=1   - wskazuje ze nalezy przesunac wiazke
 *******************************************************************************/
void PHOTON_DIFFUSION_2D:: calculate_new_position(){
	
	double wsp_abs=layers_data[beam.layer][0];
	double wsp_scat=layers_data[beam.layer][1];
	double wsp=wsp_abs+wsp_scat;
	double s=-log(uniform())/wsp;
	beam.x_new=beam.x+beam.rx*s;	
	beam.y_new=beam.y+beam.ry*s;
	beam.length=1;
	return;
}


/********************************************************************************
 * 
 * 						ROULETTE
 * 
 ********************************************************************************/
void PHOTON_DIFFUSION_2D:: roulette(){
	if(beam.w<w_min && beam.alive==true){
		double p=uniform();
		if(p<=p_min){
			beam.w=beam.w/p;
		}else{
			beam.w=0;
			beam.alive=false;
		}
	}
	return;
}

void PHOTON_DIFFUSION_2D::fill_m_layer(int m, double a, double s, double d, double g, double n)
{
	if(m > nlayers || m < 0){
		printf("Wrong layer\n");
		return;
	}

	layers_data[m][0] = a;
	layers_data[m][1] = s; 
	layers_data[m][2] = d;
	layers_data[m][3] = g;
	layers_data[m][4] = n;
}

void PHOTON_DIFFUSION_2D::export_to_files(string name)
{
	ofstream file;
	file.open("data/abs_"+name);

	// export absorbtion data
	file << "#\t" << nx << "\t" << ny << "\n";
	for(long unsigned int k = 0; k < absorption.size(); k++){
        for(long unsigned int l = 0; l < absorption[k].size(); l++){
            file << absorption[k][l] << "\t";
        }
        file << "\n";
    }

	file.close();

	// export transmitance data
	file.open("data/trans"+name);

	for(long unsigned int k = 0; k < transmittance.size(); k++){
		file << transmittance[k] << "\n";
	}
	
	file.close();

	// export reflectance data
	file.open("data/refl"+name);

	for(long unsigned int k = 0; k < reflectance.size(); k++){
		file << reflectance[k] << "\n";
	}
	
	file.close();
}