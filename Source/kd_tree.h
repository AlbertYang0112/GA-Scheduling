/*
 * kd_tree.h
 *
 *  Created on: 2019年5月24日
 *      Author: 15619
 */

#ifndef KD_TREE_H_
#define KD_TREE_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <ctime>

typedef struct Point
{
    double x;
    double y;
    double t;
    int affiliate;
} Point;

class KD_tree{
private:
	Point* KD;
	int length;
	static double distance(Point a, Point b, int layer=-1);
	static bool cmpx(Point a, Point b);
	static bool cmpy(Point a, Point b);
	static void point_sort(int start, int end, int dim, Point* points);
	void build(int start, int end, int layer);
	Point find_near(Point a, int layer, int start, int end, bool log=false);
public:
	KD_tree(int num, Point* points);
	Point find(Point a, bool log=false);
	~KD_tree();
};

const Point far_point{LONG_MAX, LONG_MAX, LONG_MAX, -1};

#endif /* KD_TREE_H_ */
