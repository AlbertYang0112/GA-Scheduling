#include <iostream>
#include <math.h>
using namespace std;

double mod(double x,double y){
	double z=fmod(x,y);
	if(z<0){
		z=z+y;
	}
	return z;
}

double dubins(double x1,double y1,double alpha,double x2,double y2,double beta,double r){
	//convert
	double pi=M_PI;
	double theta=atan((y2-y1)/(x2-x1));
	double d=(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)))/r;
	double alpha_0=alpha-theta;
	double beta_0=beta-theta;

	//LSL
	double T_lsl=mod(atan((cos(beta_0)-cos(alpha_0))/(d+sin(alpha_0)-sin(beta_0))),(2*pi)) - alpha_0;
	double P_lsl=sqrt(2+d*d-2*cos(alpha_0-beta_0)+2*d*(sin(alpha_0)-sin(beta_0)));
	double Q_lsl=beta_0 - mod(atan((cos(beta_0)-cos(alpha_0))/(d+sin(alpha_0)-sin(beta_0))),(2*pi));
	double LSL=T_lsl+P_lsl+Q_lsl;
	double length=LSL;
//	cout<<T_lsl<<" "<<P_lsl<<" "<<Q_lsl<<endl;
//	cout<<LSL*r<<endl;

	//RSR
	double T_rsr=alpha_0-mod(atan((cos(alpha_0)-cos(beta_0))/(d+sin(beta_0)-sin(alpha_0))),(2*pi));
	double P_rsr=sqrt(2+d*d-2*cos(alpha_0-beta_0)+2*d*(sin(beta_0)-sin(alpha_0)));
	double Q_rsr=mod((-beta_0),(2*pi)) + mod(atan((cos(alpha_0)-cos(beta_0))/(d-sin(alpha_0)+sin(beta_0))),(2*pi));
	double RSR=T_rsr+P_rsr+Q_rsr;
	if(RSR>0)
		length=min(length,RSR);
//	cout<<T_rsr<<" "<<P_rsr<<" "<<Q_rsr<<endl;
//	cout<<RSR*r<<endl;

	//LSR(???)
	double P_lsr=sqrt(-2+d*d+2*cos(alpha_0-beta_0)+2*d*(sin(alpha_0)+sin(beta_0)));
	double T_lsr=mod((atan((-cos(alpha_0)-cos(beta_0))/(d+sin(alpha_0)+sin(beta_0)))-alpha_0-atan(-2/P_lsr)),(2*pi));
	double Q_lsr=-mod(beta_0,(2*pi))+atan((-cos(alpha_0)-cos(beta_0))/(d+sin(alpha_0)+sin(beta_0)))-mod(atan(-2/P_lsr),(2*pi));
	double LSR=P_lsr+T_lsr+Q_lsr;
//	if(LSR>0)
//		length=min(LSR,length);
//	cout<<T_lsr<<" "<<P_lsr<<" "<<Q_lsr<<endl;
//	cout<<LSR*r<<endl;

	//RSL(???)
	double P_rsl=sqrt(d*d-2+2*cos(alpha_0-beta_0)-2*d*(sin(alpha_0)+sin(beta_0)));
	double T_rsl=alpha_0-atan((cos(alpha_0)+cos(beta_0))/(d-sin(alpha_0)-sin(beta_0)))+mod((atan(2/P_rsl)),(2*pi));
	double Q_rsl=mod(beta_0,(2*pi))-atan((cos(alpha_0)+cos(beta_0))/(d-sin(alpha_0)-sin(beta_0)))+mod(atan(2/P_rsl),(2*pi));
	double RSL=P_rsl+T_rsl+Q_rsl;
//	if(RSL>0)
//		length=min(length,RSL);
//	cout<<T_rsl<<" "<<P_rsl<<" "<<Q_rsl<<endl;
//	cout<<RSL*r<<endl;

	//RLR(???)
	double P_rlr=acos((6-d*d+2*cos(alpha_0-beta_0)+2*d*(sin(alpha_0)-sin(beta_0)))/8);
	double T_rlr=alpha_0-atan((cos(alpha_0)-cos(beta_0))/d-sin(alpha_0)+sin(beta_0))+mod((P_rlr/2),(2*pi));
	double Q_rlr=alpha_0-beta_0-T_rlr+mod(P_rlr,(2*pi));
	double RLR=P_rlr+T_rlr+Q_rlr;
	if(RLR>0)
		length=min(length,RLR);
//	cout<<T_rlr<<" "<<P_rlr<<" "<<Q_rlr<<endl;
//	cout<<RLR*r<<endl;

	//LRL(???)
	double P_lrl=mod(acos((6-d*d+2*cos(alpha_0-beta_0)+2*d*(sin(alpha_0)-sin(beta_0)))/8),(2*pi));
	double T_lrl=mod((-alpha_0+atan((-cos(alpha_0)+cos(beta_0))/d+sin(alpha_0)-sin(beta_0))+(P_lrl/2)),(2*pi));
	double Q_lrl=-alpha_0+mod(beta_0,(2*pi))-T_lrl+mod((2*P_lrl),(2*pi));
	double LRL=P_lrl+T_lrl+Q_lrl;
	if(LRL>0)
		length=min(length,LRL);
//	cout<<T_lrl<<" "<<P_lrl<<" "<<Q_lrl<<endl;
//	cout<<LRL*r<<endl;

	return (length*r);
}

//int main()
//{
//	int a1,a2,b1,b2,the1,the2,r;
//	cout<<"please input the position of UAV"<<endl;
//	cin>>a1>>b1>>the1;
//	cout<<"please input the position of target"<<endl;
//	cin>>a2>>b2>>the2;
//	cout<<"please input the radius of UAV"<<endl;
//	cin>>r;
//	cout<<"the length is:"<<endl;
//	cout<<dubins(a1,b1,a2,b2,the1,the2,r);
//	return 0;
//}
