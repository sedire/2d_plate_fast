#ifndef _PLATE_2D_ORTHOBUILDER_
#define _PLATE_2D_ORTHOBUILDER_ 1

#include "plate_var_types.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::bad_alloc;

class SolInfo
{
public:
	SolInfo();
	~SolInfo();

	vector<PL_NUM> o;		//omega matrix to restore the solution
	vector<vector<PL_NUM>> zi;		//basis vectors of the solution, orthonormalized
	//vector<PL_NUM> z5;

	vector<PL_NUM> C;

	void setup( int _varNum );
	void flushO();
};

class OrthoBuilder
{
public:
	vector<SolInfo> solInfoMap;
	OrthoBuilder();
	virtual ~OrthoBuilder();
	virtual void setParams( int _Km );
	virtual void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt ) {};
	virtual void buildSolution( vector<VarVect>* _mesh ) {};
	virtual void flushO( int x );
	virtual void setInitVects( const vector<PL_NUM>& N1, const vector<PL_NUM>& N2, const vector<PL_NUM>& N3, const vector<PL_NUM>& N4, const vector<PL_NUM>& N5 );
	void LUsolve( vector<vector<PL_NUM>>& AA, vector<PL_NUM>& ff, vector<PL_NUM>* xx );
	PL_NUM z5[NODES_ON_Y][EQ_NUM * NUMBER_OF_LINES];
protected:
	int eq_num;
	int Km;
	vector<vector<PL_NUM>> LL;
	vector<vector<PL_NUM>> UU;
};

class OrthoBuilderGodunov : public OrthoBuilder			//not working yet!
{
public:
	OrthoBuilderGodunov( int _varNum );
	~OrthoBuilderGodunov() {};
	void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );
};

class OrthoBuilderGSh : public OrthoBuilder			//use this
{
public:
	OrthoBuilderGSh( int _varNum);
	~OrthoBuilderGSh() {};
	void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );
};

#endif