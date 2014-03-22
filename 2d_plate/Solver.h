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
#include <Eigen/Eigen>

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

	void dump_Amir_sol();

	PL_NUM increaseTime();
	PL_NUM getCurTime();

private:
	const PL_NUM E1;				//Young's modulus
	const PL_NUM E2;				//Young's modulus
	const PL_NUM nu21;			//Poisson's ratio	
	const PL_NUM nu23;			//Poisson's ratio	
	const PL_NUM rho;				//composite's density
	const PL_NUM G23;				//shear modulus

	const PL_NUM B11;
	const PL_NUM B22;
	const PL_NUM B12;
	const PL_NUM B66;
	const PL_NUM By0;
	const PL_NUM By1;				// in considered boundary-value problem
	const PL_NUM By2;   
	const PL_NUM betta;		//parameter at Newmark's time integration scheme

	const PL_NUM mu;
	const PL_NUM sigma_x;			//electric conductivity
	const PL_NUM sigma_x_mu;
	const PL_NUM sigma_y;
	const PL_NUM sigma_y_mu;
	const PL_NUM sigma_z;

	const PL_NUM J0;
	const PL_NUM tauC;
	const PL_NUM tauP;
	PL_NUM omega;
	const PL_NUM  p0;				//constant mechanical load
	const PL_NUM impRadSq;

	const PL_NUM eps_0;
	const PL_NUM eps_x;
	const PL_NUM eps_x_0;

	const PL_NUM hp;				//thickness of the plate
	const PL_NUM ap;				//width of the plate
	const PL_NUM bp;				//length of the plate

	const int Km;				//number of steps by space
	const int nx;				//number of lines in Method of Lines
	const int eq_num;				//number of equations in main system
	const int varNum;		//number of variables on the whole line. equals eq_num * number_of_lines

	const int maxNewtonIt;
	int newtonIt;

	const PL_NUM dx;
	const PL_NUM dy;

	const PL_NUM dt;			//time step
	PL_NUM cur_t;
	int curTimeStep;

	const PL_NUM al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density

	vector<VarVect> mesh;		//2d mesh for solution.
	//vector<vector<PL_NUM>> matr_A;		//matrix A for the nonlinear system at certain t and x
	//vector<PL_NUM> vect_f;		//vector f on right part of nonlinear system at certain t and x
	//vector<PL_NUM> newmark_A;
	//vector<PL_NUM> newmark_B;
	PL_NUM matr_A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM vect_f[EQ_NUM * NUMBER_OF_LINES];		//vector f on right part of nonlinear system at certain t and x
	PL_NUM newmark_A[EQ_NUM * NUMBER_OF_LINES];
	PL_NUM newmark_B[EQ_NUM * NUMBER_OF_LINES];

	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES, RowMajor> Ma;

	RungeKutta* rungeKutta;
	OrthoBuilder* orthoBuilder;

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