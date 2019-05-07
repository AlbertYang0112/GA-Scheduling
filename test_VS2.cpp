// test_VS2.cpp : 定义控制台应用程序的入口点。
//

//============================================================================
// Name        : test.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : test for calculate dubins length
//============================================================================

#include <iostream>
#include <math.h>
#include <algorithm>
#include "dubins.h"
#include <limits.h>


using namespace std;


typedef struct collision_POINT {
	double x1;
	double y1;
	double x2;
	double y2;
	double theta1;
	double theta2;
	double t;
	collision_POINT *next;
};



class COLLISION {
private:
	typedef struct POINT {
		double x;
		double y;
		double theta;
		double t;
	}POINT;
	POINT *points;
	int point_num;
	double min_distance;
	double sample_rates;

	bool cmpx(POINT a, POINT b) {
		return (a.x>b.x);
	}
	bool cmpy(POINT a, POINT b) {
		return (a.y>b.y);
	}
	double distance(POINT *a, POINT *b) {
		if ((b == NULL) || (a == NULL))
			return LONG_MAX;
		return((a->x - b->x)*(a->x - b->x)) + ((a->y - b->y)*(a->y - b->y));
	}
	void collision_add(POINT a, POINT b) {
		if ((a.t - b.t) <= sample_rates) {
			collision_POINT *x = new collision_POINT;
			x->x1 = a.x;
			x->x2 = b.x;
			x->y1 = a.y;
			x->y2 = b.y;
			x->t = a.t;
			x->next = NULL;
			p->next = x;
			p = x;
		}
	}

	int collision_detect(int level, int num, POINT *father, POINT *a) {
		//kd-tree build
		POINT current_point;
		if (num == 1) {
			current_point = a[0];
			if (distance(&current_point, father)<min_distance) {
				collision_add(current_point, *father);
			}
			return 0;
		}
		if (num == 0) {
			return 0;
		}

		if (level % 2) {
			sort(a, a + num, this->cmpy);
		}
		else {
			sort(a, a + num, this->cmpx);
		}
		current_point = a[num / 2];
		if (distance(&current_point, father)<min_distance) {
			collision_add(current_point, *father);
		}
		//generate sons
		POINT *left_child = new POINT[num / 2];
		POINT *right_child = new POINT[(num - 1) / 2];
		//这里好像可以更高效？
		for (int i = 0; i<num / 2; i++)
			left_child[i] = a[i];
		for (int i = (1 + (num / 2)); i<num; i++)
			right_child[i - (1 + (num / 2))] = a[i];
		collision_detect(level + 1, num / 2, &current_point, left_child);
		collision_detect(level + 1, (num - 1) / 2, &current_point, right_child);
		return 0;
	}


public:
	collision_POINT *collision_points, *p;


	int RecallForCollision(double q[3], double x, void* user_data) {
		POINT a;
		a.x = q[0];
		a.y = q[1];
		a.theta = q[2];
		a.t = x;
		points[point_num] = a;
		point_num++;
		return 0;
	}
	COLLISION(DubinsPath *path, int path_num, double sample_rate, double min_radius) {
		//create all the points
		double length = 0;
		for (int i = 0; i<path_num; i++) {
			length += dubins_path_length(&path[i]) / sample_rate;
		}
		point_num = int(length);
		points = new POINT[point_num];
		point_num = 0;
		sample_rates = sample_rate;

		for (int i = 0; i<path_num; i++) {
			dubins_path_sample_many(&path[i], sample_rate, this->RecallForCollision, NULL);
		}
		min_distance = min_radius*min_radius;
		collision_points = NULL;
		p = collision_points;
	}

	collision_POINT *get_collision() {
		collision_detect(0, point_num, NULL, points);
		return collision_points;
	}
};


void show_coll(collision_POINT *a) {
	if (a != NULL) {
		cout << a->x1 << a->x2 << endl;
	}
}

int main() {
	double q0[] = { 0,0,0 };
	double q1[] = { 4,4,3.142 };
	//double q2[] = { 0,4,1.56 };
	//double q3[] = { 4,0,1.56 };
	double turning_radius = 1.0;
	double sample_rate = 0.1;
	DubinsPath *path = new DubinsPath[2];
	dubins_shortest_path(&path[0], q0, q1, turning_radius);
	dubins_shortest_path(&path[1], q1, q0, turning_radius);
	//dubins_shortest_path( &path[2], q2, q3, turning_radius);
	//dubins_shortest_path( &path[3], q3, q1, turning_radius);
	COLLISION *A = new COLLISION(path, 2, sample_rate, turning_radius);
	collision_POINT *a;
	a = A->get_collision();
	show_coll(a);
	return 0;
}

