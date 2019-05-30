/*
 * collision.cpp
 *
 *  Created on: 2019年5月25日
 *      Author: 15619
 */

#include "collision.h"

using namespace std;

void Collision::Add_point(Point a){
	points[point_num] = a;
	point_num += 1;
	if(point_num>=max_point)
		throw "Too many points!";
}

int Collision::recall_for_dubins(double q[3], double x, void* user_data){
	Point a;
	Collision* p = static_cast<Collision *>(user_data);
	a.x = q[0];
	a.y = q[1];
	a.t = x;
	a.affiliate = p->who;
	p->Add_point(a);
	//cout<<"points:"<<p->who<<" "<<a.x<<" "<<a.y<<endl;
	return 0;
}

double Collision::distance(Point a, Point b){
	double dis=0;
	dis=(a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);
	dis=sqrt(dis);
	if(abs(a.t-b.t)>sample_rate)
		dis=LONG_MAX;
	if(a.affiliate==b.affiliate)
		dis=LONG_MAX;
	return dis;
}

Collision::Collision(DubinsPath* path, int num_path, double sample, double min_distance){
	sample_rate=sample;
	min_dis=min_distance;
	points=new Point[max_point];
	co_points=new co_Point[max_collision];
	point_num=0;
	collision_num=0;
	for(int i=0;i<num_path;i++){
		who=i;
		dubins_path_sample_many( &path[i], sample_rate, recall_for_dubins, this);
	}
	KD=new KD_tree(point_num, sample, points);
	Point temp;
	for(int i=0;i<point_num;i++){
		temp=KD->find(points[i]);
		if(distance(temp,points[i])<min_dis){
			co_Point t{points[i],temp};
			co_points[collision_num]=t;
			collision_num++;
			if(collision_num>max_collision)
				throw "Too many collisions!";
		}
	}
}

Collision::~Collision(){
	delete[] points;
	delete[] co_points;
	delete KD;
}

int Collision::get_collision_num(){
	return collision_num;
}

co_Point* Collision::get_collisions(){
	return co_points;
}

void Collision::debug(){
	for(int i=0;i<point_num;i++){
		cout<<"points: affiliate:"<<points[i].affiliate<<" x:"<<points[i].x<<" y:"<<points[i].y<<" t:"<<points[i].t<<endl;
	}
}

void Collision::debug2(Point a){
	Point temp=KD->find(a);
	cout<<"points: affiliate:"<<temp.affiliate<<" x:"<<temp.x<<" y:"<<temp.y<<" t:"<<temp.t<<endl;
}
