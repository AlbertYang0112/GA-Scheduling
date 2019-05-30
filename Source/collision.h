/*
 * collision.h
 *
 *  Created on: 2019年5月25日
 *      Author: 15619
 */

#ifndef COLLISION_H_
#define COLLISION_H_

#include "kd_tree.h"
extern "C" {
#include "dubins.h"
};

const int max_collision=10000;
const int max_point=100000;

typedef struct co_Point{
	Point a;
	Point b;
}co_Point;

class Collision{
private:
	//Variable
	Point* points;
	co_Point* co_points;
	int collision_num;
	int point_num;
	//A temporary variable for affiliation
	int who;
	double sample_rate;
	double min_dis;
	KD_tree* KD;

	//Function
	double distance(Point a, Point b);
	void Add_point(Point a);
	void Add_collision(Point a, Point b);
	static int recall_for_dubins(double q[3], double x, void* user_data);
public:
	Collision(DubinsPath* paths, int num_path, double sample, double min_distance);
	~Collision();
	co_Point* get_collisions();
	int get_collision_num();
	void debug();
	void debug2(Point a);
};



#endif /* COLLISION_H_ */
