#ifndef  _PLATE_2D_SOLVER_THERMO_
#define _PLATE_2D_SOLVER_THERMO_ 1

#include "plate_var_types.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;

class SolverThermo
{
public:
	//const string name;
	SolverThermo();
	virtual ~SolverThermo();

//	virtual int setParams() = 0;
	virtual int doStep() = 0;
	virtual void dumpSol() = 0;
	virtual void dumpSolBin() = 0;
	virtual PL_NUM* getSolDx() = 0;
protected:
	PL_NUM aa;		//lenght of the plate in x-dir
	PL_NUM bb;		//width of the plate in y-dir
	PL_NUM hh;		//thickness of the plate

	PL_NUM J0;		//current density
	PL_NUM tauExp;
	PL_NUM tauSin;
	PL_NUM Jx;		//current at the actual time
	PL_NUM Tamb;	//ambient temperature

	PL_NUM Rc;		//contact resistance

	//plate material params
	PL_NUM kx;
	PL_NUM hInf;
	PL_NUM sigmaX;
	PL_NUM rho;
	PL_NUM cc;

	PL_NUM AA;
	PL_NUM BB;
	PL_NUM CC;
	PL_NUM CC1;

	int NN;		//number of nodes in the plate (-1). This should be a multiple of the number of lines +2 in the main problem
	PL_NUM dt;		//time step, should be the same as in the main code
	PL_NUM dx;		//spatial step

	PL_NUM* plateNodesN;
	PL_NUM* plateNodesN1;
	PL_NUM* solDx;

	//a matrix and the rhs of the system of equations to solve
	PL_NUM* ai;
	PL_NUM* bi;
	PL_NUM* ci;
	PL_NUM* di;

	PL_NUM curTime;
	int curTimeStep;

	virtual int fillABSystem() = 0;
	virtual int updateSystem() = 0;
};

class SolverThermoWElectrodes : public SolverThermo
{
public:
	SolverThermoWElectrodes( int _NN,
				PL_NUM _aa, PL_NUM _bb, PL_NUM _aaEl, PL_NUM _hh, PL_NUM _hhEl,
				PL_NUM _J0, PL_NUM _tauExp, PL_NUM _tauSin, PL_NUM _Tamb,
				PL_NUM _Rc,
				PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
				PL_NUM _kEl, PL_NUM _hInfEl, PL_NUM _sigmaEl, PL_NUM _rhoEl, PL_NUM _ccEl,
				PL_NUM _dt );
	//int setParams( int _NN,
	//			PL_NUM _aa, PL_NUM _bb, PL_NUM _aaEl, PL_NUM _hh, PL_NUM _hhEl,
	//			PL_NUM _J0, PL_NUM _Tamb,
	//			PL_NUM _Rc,
	//			PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
	//			PL_NUM _kEl, PL_NUM _hInfEl, PL_NUM _sigmaEl, PL_NUM _rhoEl, PL_NUM _ccEl,
	//			PL_NUM _dt );
	~SolverThermoWElectrodes();
	int doStep();
	void dumpSol();
	void dumpSolBin();
	PL_NUM* getSolDx();
protected:
	PL_NUM aaEl;	//length of the electrode in x-dir
	PL_NUM hhEl;	//thickness of the electrode
	//electrode material params
	PL_NUM kEl;
	PL_NUM hInfEl;
	PL_NUM sigmaEl;
	PL_NUM rhoEl;
	PL_NUM ccEl;
	//coeffs in the PDE
	PL_NUM AAel;
	PL_NUM BBel;
	PL_NUM CCel;
	PL_NUM CCel1;

	int MM;		//number of nodes in the electrode (-1). We want the same dx everywhere, so it will depend on NN

	PL_NUM TelRightG;	//right ghost node for the electrode at time n
	PL_NUM TleftG;		//left ghost node for the plate at time n
	PL_NUM TelRightG1;	//right ghost node for the electrode at time n + 1
	PL_NUM TleftG1;		//left ghost node for the plate at time n + 1

	PL_NUM* electrNodesN;
	PL_NUM* electrNodesN1;

	int fillABSystem();
	int updateSystem();
};

class SolverThermoConstFlow : public SolverThermo
{
public:
	SolverThermoConstFlow( int _NN,
				PL_NUM _aa, PL_NUM _bb, PL_NUM _hh,
				PL_NUM _J0, PL_NUM _tauExp, PL_NUM _tauSin, PL_NUM _Tamb,
				PL_NUM _Rc,
				PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
				PL_NUM _dt );
	~SolverThermoConstFlow();
	//int setParams( int _NN,
	//			PL_NUM _aa, PL_NUM _bb, PL_NUM _hh,
	//			PL_NUM _J0, PL_NUM _Tamb,
	//			PL_NUM _Rc,
	//			PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
	//			PL_NUM _dt ) {};
	int doStep();
	void dumpSol();
	void dumpSolBin();
	PL_NUM* getSolDx() { cout << "WARNING: I am void\n"; return 0; };
protected:
	int fillABSystem();
	int updateSystem();
};

#endif