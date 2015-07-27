#ifndef _PLATE_2D_ORTHOBUILDER_
#define _PLATE_2D_ORTHOBUILDER_ 1

#include "plate_var_types.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>
#include "omp.h"
#include <Eigen/Eigen>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::bad_alloc;
using std::complex;
using namespace Eigen;

class SolInfo
{
public:
	SolInfo();
	~SolInfo();

	vector<PL_NUM> o;		//omega matrix to restore the solution
	//vector<vector<PL_NUM>> zi;		//basis vectors of the solution, orthonormalized
	//vector<PL_NUM> z5;

	//vector<PL_NUM> C;

	void setup( int _varNum );
	void flushO();
};

class OrthoBuilder
{
public:
	vector<SolInfo> solInfoMap;
	OrthoBuilder( int _varNum, int Km );
	virtual ~OrthoBuilder();
	virtual void setParams();
	virtual void orthonorm( int baseV, int n, PL_NUM* NtoOrt ) {};
	virtual void buildSolution( vector<VarVect>* _mesh ) {};
	virtual void flushO( int x );
	virtual void setInitVects( const vector<PL_NUM>& N1, const vector<PL_NUM>& N2, const vector<PL_NUM>& N3, const vector<PL_NUM>& N4, const vector<PL_NUM>& N5 );
	virtual void setOrthoDoneInfo( int y );
	virtual void resetOrthoDoneInfo();
	inline virtual void setNextSolVects( int n, const PL_NUM (*decompVect)[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] ) {};
	inline virtual int checkOrtho( int n, 
									PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
									PL_NUM vectSetOrig[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] ) { return 1; };

	inline PL_NUM getInfNorm( PL_NUM* vect, int vectSize );
	void LUsolve( vector<vector<PL_NUM>>& AA, vector<PL_NUM>& ff, vector<PL_NUM>* xx );
	//PL_NUM zi[NODES_ON_Y][EQ_NUM * NUMBER_OF_LINES / 2][EQ_NUM * NUMBER_OF_LINES];
	vector<vector<vector<PL_NUM> > > zi;
	vector<vector<PL_NUM> > z5;//[NODES_ON_Y][EQ_NUM * NUMBER_OF_LINES];
protected:
	const int varNum;
	const int Km;
	vector<vector<PL_NUM>> LL;
	vector<vector<PL_NUM>> UU;
	vector<bool> orthoDone;
	PL_NUM omega2[EQ_NUM * NUMBER_OF_LINES];
	PL_NUM omegaPar[NUM_OF_THREADS];
	PL_NUM Cx[EQ_NUM * NUMBER_OF_LINES / 2];
	PL_NUM Cx1[EQ_NUM * NUMBER_OF_LINES / 2];
private:
	OrthoBuilder();
};

class OrthoBuilderGodunov : public OrthoBuilder			//not working yet!
{
public:
	OrthoBuilderGodunov( int _varNum, int _Km );
	~OrthoBuilderGodunov() {};
	void orthonorm( int baseV, int n, PL_NUM* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );
};

class OrthoBuilderGSh : public OrthoBuilder			//use this
{
public:
	OrthoBuilderGSh( int _varNum, int Km );
	~OrthoBuilderGSh() {};
	void orthonorm( int baseV, int n, PL_NUM* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );

	void calcScalarProdsPar( int baseV, int n, vector<PL_NUM>* NtoOrt );
	void calcScalarProdsPar2( int baseV, int n, vector<PL_NUM>* NtoOrt );

	inline void setNextSolVects( int n, const PL_NUM (*decompVect)[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] );
	inline int checkOrtho( int n, 
							PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
							PL_NUM vectSetOrig[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] );
private:
	OrthoBuilderGSh();
};

#endif