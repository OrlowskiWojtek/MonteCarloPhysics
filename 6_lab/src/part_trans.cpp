#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<cmath>

#include "part_trans.hpp"

using namespace std;

/*
 * obrot wektora r1=[x,y,z] wokol osi t=[tx,ty,tz] o kat teta przy uzyciu kwaternionow
 * 
 * r2=p*r1*p^{*}
 * 
 */
void vector_rotation(double tx, double ty, double tz, double & x, double & y, double & z, double teta){
	double p0,p1,p2,p3;
	double q0,q1,q2,q3;
	double r0,r1,r2,r3;
	
	double t=sqrt(tx*tx+ty*ty+tz*tz);
	double cs=cos(teta/2);
	double sn=sin(teta/2);
	// kwaternion obrotu
	p0=cs;
	p1=tx/t*sn;
	p2=ty/t*sn;
	p3=tz/t*sn;
	//kwaternion obracanego wektora
	r0=0.;
	r1=x;
	r2=y;
	r3=z;
	// q=r_1*p^{*}
	q0=r0*p0-r1*(-p1)-r2*(-p2)-r3*(-p3);
	q1=r1*p0+r0*(-p1)-r3*(-p2)+r2*(-p3);
	q2=r2*p0+r3*(-p1)+r0*(-p2)-r1*(-p3);
	q3=r3*p0-r2*(-p1)+r1*(-p2)+r0*(-p3);
	
	// r_2=p*q
	r0=p0*q0-p1*q1-p2*q2-p3*q3;
	r1=p1*q0+p0*q1-p3*q2+p2*q3;
	r2=p2*q0+p3*q1+p0*q2-p1*q3;
	r3=p3*q0-p2*q1+p1*q2+p0*q3;
	
	// zwracamy wspolrzedne wektora obroconego
	x=r1;
	y=r2;
	z=r3;
	
	
	return;
}


/*
 * przemieszczamy czastke z punktu P1 do P2, jesli trajektoria przecina brzeg (okrag) w punkcie P3
 * to czastke umieszczamy na brzegu w P1=P3 a nowym punktem koncowym jest P2=P4 ktory powstaje 
 * z obrocenia P1-P3 wokol wektora (0,0)-P3  i przeskalowania go do dlugosci |P2-P3|
 * 
 * jesli odcinek P1-P2 przecina obszar absorbera (okrag o promieniu Ra) to ustawiamy flage: exist=0    
 *
 * zmienna length okresla dlugosc nowego odcinka P1-P2, 
 * jesli length=0 to konczymy przesuwanie czastki
 *
 * WEJSCIE:
 * (x1,y1) - poczatek odcinka
 * (x2,y2) - koniec odcinka
 * (xr,yr), Rr - polozenie srodka i promien obszaru/okregu
 * (xa,ya), Ra - polozenie srodka i promien absorbera
 * 
 * WYJSCIE:
 * (x1,y1) - poczatek sciezki po odbiciu
 * (x2,y2) - koniec  sciezki po odbiciu
 * exist=0:  0-przejscie przez absorber 
 * length - dlugosc odbitej lub koncowej sciezki (jesli nie ma odbicia)
 * 
 * 
 */


void particle_translation(double & x1,double & y1, double & x2,double & y2, 
				  double xr,double yr, double Rr, double xa,double ya, double Ra, 
				  int & exist ,double & length)
{
	double a,b,c,delta,beta1,beta2,alfa1,alfa2,alfa;
	/*
	 * liczymy bete: jesli pierwiastki sa rzeczywiste, dodatnie i beta=[0,1] to czastka przechodzi przez obszar absorpcji 
	 * 
	 */	
		a=pow(x2-x1,2)+pow(y2-y1,2);
		b=2*( (x2-x1)*(x1-xa)+(y2-y1)*(y1-ya) );
		c=pow(x1-xa,2)+pow(y1-ya,2)-pow(Ra,2);
		delta=b*b-4*a*c;
		if(delta>=0){
			beta1=(-b-sqrt(delta))/2/a;
			beta2=(-b+sqrt(delta))/2/a;
			if(beta1>=0 && beta1<=1 || beta2>=0 && beta2<=1) exist=0; //wskazuje ze czastka znika
		}
	
	
	/*
	 * liczymy alfe: jesli jest jeden pierwiastek rzeczywisty, dodatni i alfa=[0,1] to czastka przechodzi przez brzeg obszaru/okregu 
	 * 
	 */   
	      alfa=-1.0; //blokada
	      
		a=pow(x2-x1,2)+pow(y2-y1,2);
		b=2*( (x2-x1)*(x1-xr)+(y2-y1)*(y1-yr) );
		c=pow(x1-xr,2)+pow(y1-yr,2)-pow(Rr,2);
		delta=b*b-4*a*c;
		if(delta>=0){
			alfa1=(-b-sqrt(delta))/2/a;
			alfa2=(-b+sqrt(delta))/2/a;
			/*
			 * sprawdzamy czy ktorys pierwiastek pasuje 
			 */
			if(alfa1>=0 && alfa1 <=1){
				alfa=alfa1;
			}else if(alfa2>=0 && alfa2 <=1){
				alfa=alfa2;
			}
		}
	
		
	/*
	 * alfa<0:  brak przeciecia z brzegiem - koniec odcinka jest w obszarze
	 * alfa=[0,1]: wyznaczamy punkt przeciecia z okregiem P3->(x1,y1) i nowy punkt koncowy P4->(x2,y2)
	 * 
	 */	
	if(alfa<0){
		x1=x2;
		y1=y2;
		length=0;
		return;
	} else{
		double norm;
		double x3,y3,z3,x4,y4,z4;
		x3=x1+alfa*(x2-x1); //punkt na brzegu
		y3=y1+alfa*(y2-y1);
		z3=0.;
		
		double x13=x1-x3; //wektor ktory obracamy: P1-P3 - to bedzie nowy kierunek
		double y13=y1-y3;
		double z13=0;
		norm=sqrt(pow(x13,2)+pow(y13,2)+pow(z13,2));
		x13=x13/norm;
		y13=y13/norm;
		z13=z13/norm;
		
		double tx=0-x3; //wspolrzedne osi obrotu - wektor: (0,0)-P3
		double ty=0-y3;
		double tz=0-z3;
		
		// obrot P13 wokol osi (tx,ty,tz) o kat Pi
		vector_rotation(tx,ty,tz,x13,y13,z13,M_PI);
		
		
		length=sqrt(pow(x2-x3,2)+pow(y2-y3,2)); //dlugosc wektora po odbiciu od brzegu
		// nowy punkt P1=P3
		x1=x3+x13*1.0E-6; //przesuwamy MINIMALNIE punkt z brzegu do srodka - aby nie bylo 2 przeciec w nastepnej iteracji
		y1=y3+y13*1.0E-6;
		// nowy punkt P2=P4=P3+length*[x13,y13]
		x2=x3+length*x13;
		y2=y3+length*y13;
		
		return;
	}	
	return;
}


