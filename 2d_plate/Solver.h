#ifndef _PLATE_2D_SOLVER_
#define _PLATE_2D_SOLVER_ 1

#include <time.h>
#include "plate_var_types.h"
#include "RungeKutta.h"
#include "OrthoBuilder.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <omp.h>

using std::cout;
using std::vector;
using std::ofstream;
using std::string;
using std::stringstream;
using std::setprecision;

class Solver
{
public:
	Solver();
	~Solver();

	void setTask();

	void calc_nonlin_system_run_test( long  _x, long _t );

	void pre_step();
	void do_step();
	void dump_sol();
	void dump_whole_sol( int var );
	void dump_check_sol();
	void dump_check_sol2D();
	void dump_border_vals();


	void testLUsolve();		//DELETEME

	PL_NUM cur_t;
	PL_NUM dt;			//time step
	int curTimeStep;
private:
	PL_NUM E1;				//Young's modulus
	PL_NUM E2;				//Young's modulus
	PL_NUM nu21;			//Poisson's ratio	
	PL_NUM nu23;			//Poisson's ratio	
	PL_NUM rho;				//composite's density
	PL_NUM G23;				//shear modulus

	PL_NUM mu;

	PL_NUM sigma_x;			//electric conductivity
	PL_NUM sigma_x_mu;
	PL_NUM sigma_y;
	PL_NUM sigma_y_mu;
	PL_NUM sigma_z;

	PL_NUM hp;				//thickness of the plate
	PL_NUM ap;				//width of the plate
	PL_NUM bp;				//length of the plate

	int eq_num;				//number of equations in main system

	PL_NUM J0;
	PL_NUM omega;
	PL_NUM  p0;				//constant mechanical load

	PL_NUM B11;
	PL_NUM B22;
	PL_NUM B12;
	PL_NUM B66;
	PL_NUM By0;
	PL_NUM By1;                                      // in considered boundary-value problem
	PL_NUM By2;   

	PL_NUM eps_0;
	PL_NUM eps_x;
	PL_NUM eps_x_0;

	int Km;				//number of steps by space
	int nx;				//number of lines in Method of Lines

	int varNum;
	int newtonIt;
	int maxNewtonIt;

	PL_NUM dx;
	PL_NUM dy;

	PL_NUM al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density
	PL_NUM betta;		//parameter at Newmark's time integration scheme

	vector<VarVect> mesh;		//2d mesh for solution.
	//vector<vector<PL_NUM>> matr_A;		//matrix A for the nonlinear system at certain t and x
	//vector<PL_NUM> vect_f;		//vector f on right part of nonlinear system at certain t and x
	//vector<PL_NUM> newmark_A;
	//vector<PL_NUM> newmark_B;
	//PL_NUM matr_A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM matr_A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM vect_f[EQ_NUM * NUMBER_OF_LINES];		//vector f on right part of nonlinear system at certain t and x
	PL_NUM newmark_A[EQ_NUM * NUMBER_OF_LINES];
	PL_NUM newmark_B[EQ_NUM * NUMBER_OF_LINES];

	RungeKutta* rungeKutta;
	OrthoBuilder* orthoBuilder;

	void calcConsts();
	void calc_Newmark_AB( int _x, int mode );		//don't really know why we need this stuff with mode need to FIX
	void calc_system( int _x );
	void walkthrough( int mode );
	void updateDerivs();
	void dumpMatrA( int _x );

	//vector<SolInfo> solInfoMap;
	//void orthonorm( int baseV, int n );		//baseV - number of the basis vector to orthonormalize 
											//n - spatial node (i.e, x coordinate). this orthonormalizes N1, ... , N5 and builds Omega
	//void buildSolution();					//builds solution for the current time step
	int checkConv();
};

#endif