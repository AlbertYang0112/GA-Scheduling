//============================================================================
// Name        : raw_kd.cpp
// Author      : Z
// Version     : 0.1
// Copyright   : Your copyright notice
// Description : KD_tree in C++, Ansi-style
//============================================================================

#include "kd_tree.h"
#include <algorithm>
using namespace std;

KD_tree::KD_tree(int num, Point* points)
{
	this->KD=points;
	this->length=num;
	build(0, num-1, 0);
	cout<<endl;
}

KD_tree::~KD_tree()
{
}

Point KD_tree::find_near(Point a, int layer, int start, int end, bool log){
	int num=end-start+1;
	Point son;
	double best_now;
	if(num==1){
		if(this->KD[start+(num/2)].affiliate!=a.affiliate)
			return this->KD[start+(num/2)];
		else
			return far_point;

	}
	if(num==0){
		return far_point;
	}
	Point current=this->KD[start+(num/2)];
	//if affiliate to the same point, it is invalid
	if(log)
		cout<<"track:"<<current.x<<" "<<current.y<<" "<<current.affiliate<<endl;


	if(layer%2==0){
		if(a.x<current.x)
		{
			if(log)
				cout<<"left"<<layer<<endl;
			son=find_near(a, layer+1, start, (end+start-1)/2, log);
		}
		else
		{
			if(log)
				cout<<"right"<<layer<<endl;
			son=find_near(a, layer+1, (start+end+3)/2, end, log);
		}
		}
	else{
		if(a.y<current.y)
		{
			if(log)
				cout<<"left"<<layer<<endl;
			son=find_near(a, layer+1, start, (end+start-1)/2, log);
		}
		else{
			if(log)
				cout<<"right"<<layer<<endl;
			son=find_near(a, layer+1, (start+end+3)/2, end, log);
		}
	}
	//if the distance of the hyper plan and the point is smaller than current best, check another half
	best_now=distance(son,a);
	if(distance(a,current,layer)<best_now){
		Point son2;
		if(layer%2==0){
			if(a.x>=current.x){
				if(log)
					cout<<"check left"<<layer<<endl;
				son2=find_near(a, layer+1, start, (end+start-1)/2, log);
			}
			else{
				if(log)
					cout<<"check right"<<layer<<endl;
				son2=find_near(a, layer+1, (start+end+3)/2, end, log);
			}
		}
		else{
			if(a.y>=current.y){
				if(log)
					cout<<"check left"<<layer<<endl;
				son2=find_near(a, layer+1, start, (end+start-1)/2, log);
			}
			else{
				if(log)
					cout<<"check right"<<layer<<endl;
				son2=find_near(a, layer+1, (start+end+3)/2, end, log);
			}
		}
		if(distance(son2,a)<best_now){
			son=son2;
			best_now=distance(son2,a);
		}
	}
	if(current.affiliate==a.affiliate)
		current=far_point;
	if(distance(current,a)<best_now)
		return current;
	else
		return son;
}
void KD_tree::build(int start, int end, int layer){
	int num=end-start+1;
	if(num<=1)
		return;
	this->point_sort(start, end , layer, this->KD);
	build(start, (end-1+start)/2, layer+1);
	build((start+end+3)/2, end, layer+1);
}

bool KD_tree::cmpx(Point a, Point b){
	return (a.x<b.x);
}

bool KD_tree::cmpy(Point a, Point b){
	return (a.y<b.y);
}

void KD_tree::point_sort(int start, int end, int dim, Point* points){
	if(dim%2==0)
	    sort(points+start, points+end+1, cmpx);
	else
	    sort(points+start, points+end+1, cmpy);
}

double KD_tree::distance(Point a, Point b, int layer){
	if(layer==-1){
		double dis=0;
		dis=(a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);
		dis=sqrt(dis);
		return dis;
	}
	else{
		double dis=0;
		if(layer%2==0)
			dis=abs(a.x-b.x);
		else
			dis=abs(a.y-b.y);
		return dis;
	}
	return -1;
}

Point KD_tree::find(Point a, bool log){
	return this->find_near(a,0,0,this->length-1,log);
}
