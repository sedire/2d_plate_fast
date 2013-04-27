#ifndef _PLATE_2D_RUNGEKUTTA_
#define _PLATE_2D_RUNGEKUTTA_ 1

#include <vector>
#include <iostream>
#include <time.h>
#include "plate_var_types.h"

using std::cout;
using std::endl;
using std::vector;

class RungeKutta
{
public:
	RungeKutta( int _varNum );
	~RungeKutta();
	void calc( /*const vector<vector<PL_NUM>>& A*/PL_NUM A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES], PL_NUM *f, PL_NUM dx, int thrNum, int hom, vector<PL_NUM>* x );			//method for solving a system of ODE like dy/dx = Ax + f
	void calc( /*const vector<vector<PL_NUM>>& A*/PL_NUM A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES], PL_NUM *f, PL_NUM dx, int thrNum, int hom, const vector<PL_NUM>& x, vector<PL_NUM>* x1 );

private:
	RungeKutta();
	int varNum;
	PL_NUM rgk_u;
	PL_NUM rgk_v;
	PL_NUM rgk_C1, rgk_C2, rgk_C3, rgk_C4;
	PL_NUM rgk_d21, rgk_d32, rgk_d31, rgk_d43, rgk_d42, rgk_d41;

	PL_NUM f1[NUM_OF_THREADS][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM f2[NUM_OF_THREADS][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM f3[NUM_OF_THREADS][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM f4[NUM_OF_THREADS][EQ_NUM * NUMBER_OF_LINES];

	int lb[EQ_NUM * NUMBER_OF_LINES];
	int rb[EQ_NUM * NUMBER_OF_LINES];
	void buildBorders();
};

#endif