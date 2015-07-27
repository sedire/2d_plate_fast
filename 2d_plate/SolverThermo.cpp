#include "SolverThermo.h"

SolverThermo::SolverThermo() :
	plateNodesN( 0 ),
	plateNodesN1( 0 ),
	ai( 0 ),
	bi( 0 ),
	ci( 0 ),
	di( 0 )
{

}

SolverThermo::~SolverThermo()
{
	if( plateNodesN != 0 )
	{
		delete[] plateNodesN;
	}
	if( plateNodesN1 != 0 )
	{
		delete[] plateNodesN1;
	}
	if( ai != 0 )
	{
		delete[] ai;
	}
	if( bi != 0 )
	{
		delete[] bi;
	}
	if( ci != 0 )
	{
		delete[] ci;
	}
	if( di != 0 )
	{
		delete[] di;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

SolverThermoWElectrodes::SolverThermoWElectrodes( int _NN,
				PL_NUM _aa, PL_NUM _bb, PL_NUM _aaEl, PL_NUM _hh, PL_NUM _hhEl,
				PL_NUM _J0, PL_NUM _tauExp, PL_NUM _tauSin, PL_NUM _Tamb,
				PL_NUM _Rc,
				PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
				PL_NUM _kEl, PL_NUM _hInfEl, PL_NUM _sigmaEl, PL_NUM _rhoEl, PL_NUM _ccEl,
				PL_NUM _dt ) :
	//name( "withEl" ),
	electrNodesN( 0 ),
	electrNodesN1( 0 )
{
	aa = _aa;
	bb = _bb;
	hh = _hh;
	J0 = _J0;
	tauExp = _tauExp;
	tauSin = _tauSin;
	Jx = 0.0;
	Tamb = _Tamb;
	Rc = _Rc;
	kx = _kx;
	hInf = _hInf;
	sigmaX = _sigmaX;
	rho = _rho;
	cc = _cc;
	NN = _NN;
	dt = _dt;

	aaEl = _aaEl;
	hhEl = _hhEl;
	kEl = _kEl;
	hInfEl = _hInfEl;
	sigmaEl = _sigmaEl;
	rhoEl = _rhoEl;
	ccEl = _ccEl;

	//now, calculating the others
	dx = aa / 2.0 / ( NN - 1 );		//since we are solving for a half of the plate only
	MM = ceil( aaEl / dx );
	//cout << "----------------\n";
	//cout << aaEl / dx << " " << aaEl / ( aa / 2.0 / ( NN - 1 ) ) << " " << aa / aaEl << endl;
	//cout << aaEl << endl;
	if( aaEl < MM * dx )
	{
		aaEl = MM * dx;				//we are going to change (increase) the length of the electrode in x-dir a little
									//so that we can use the same spatial step dx
	}
	//cout << aaEl << endl;
	//cout << NN << " " << MM << endl;
	//cout << "----------------\n";

	electrNodesN = new PL_NUM[MM];
	electrNodesN1 = new PL_NUM[MM];
	plateNodesN = new PL_NUM[NN];
	plateNodesN1 = new PL_NUM[NN];

	for( int i = 0; i < MM; ++i )
	{
		electrNodesN[i] = 0.0;	//because we solve for dT
	}
	for( int i = 0; i < NN; ++i )
	{
		plateNodesN[i] = 0.0;	//because we solve for dT
	}
	TelRightG = 0.0;
	TleftG = 0.0;

	ai = new PL_NUM[MM + NN - 1];
	bi = new PL_NUM[MM + NN - 1];
	ci = new PL_NUM[MM + NN - 1];
	di = new PL_NUM[MM + NN - 1];

	//PL_NUM JxInit = J0 * exp( -dt / tauExp ) * sin( M_PI * dt / tauSin );

	AA = kx / rho / cc;
	BB = 2.0 * hInf / rho / cc / hh;
	CC = 0.0;
	CC1 = 0.0;//JxInit * JxInit / sigmaX / rho / cc;

	AAel = kEl / rhoEl / ccEl;
	BBel = 2.0 * hInfEl / rhoEl / ccEl / hhEl;
	CCel = 0.0;
	CCel1 = 0.0;//JxInit * JxInit / sigmaEl / rhoEl / ccEl; 

	curTime = 0.0;
	curTimeStep = 0;

	fillABSystem();
}

SolverThermoWElectrodes::~SolverThermoWElectrodes()
{
	if( electrNodesN != 0 )
	{
		delete[] electrNodesN;
	}
	if( electrNodesN1 != 0 )
	{
		delete[] electrNodesN1;
	}
}

int SolverThermoWElectrodes::fillABSystem()
{
	bi[0] = 1.0 / dt + AAel / dx / dx + BBel / 2.0 + AAel / dx * hInfEl / kEl;
	for( int i = 1; i < MM - 1; ++i )
	{
		ai[i] = -AAel / 2.0 / dx / dx;
		bi[i] = 1.0 / dt + AAel / dx / dx + BBel / 2.0;
	}
	ai[MM - 1] = -AAel / dx / dx;
	bi[MM - 1] = 1.0 / dt + AAel / dx / dx + BBel / 2.0 + AAel / AA * kx / kEl * ( 1.0 / dt + AA / dx / dx + BB / 2.0 );
	for( int i = MM; i < MM + NN - 2; ++i )
	{
		ai[i] = -AA / 2.0 / dx / dx;
		bi[i] = 1.0 / dt + AA / dx / dx + BB / 2.0;
	}
	ai[MM + NN - 2] = -AA / dx / dx;
	bi[MM + NN - 2] = 1.0 / dt + AA / dx / dx + BB / 2.0;
	
	return 0;
}

int SolverThermoWElectrodes::updateSystem()
{
	CC = CC1;
	CCel = CCel1;
	Jx = J0 * exp( -curTime / tauExp ) * sin( M_PI * curTime / tauSin );
	CC1 = Jx * Jx / sigmaX / rho / cc;
	CCel1 = Jx * Jx / sigmaEl / rhoEl / ccEl; 

	ci[0] = -AAel / dx / dx;
	di[0] = AAel / dx / dx * electrNodesN[1] + ( 1.0 / dt - AAel / dx / dx - BBel / 2.0 - AAel / dx * hInfEl / kEl ) * electrNodesN[0]
			+ 0.5 * ( CCel + CCel1 );
	for( int i = 1; i < MM - 1; ++i )
	{
		ci[i] = -AAel / 2.0 / dx / dx;
		di[i] = AAel / 2.0 / dx / dx * electrNodesN[i + 1] + ( 1.0 / dt - AAel / dx / dx - BBel / 2.0 ) * electrNodesN[i]
				+ AAel / 2.0 / dx / dx * electrNodesN[i - 1]
				+ 0.5 * ( CCel + CCel1 );
	}
	ci[MM - 1] = -AAel / dx / dx * kx / kEl;

	//PL_NUM Jx = J0 * exp( -curTime / tauExp ) * sin( M_PI * curTime / tauSin );
	di[MM - 1] = AAel / 2.0 / dx / dx * TelRightG + ( 1.0 / dt - AAel / dx / dx - BBel / 2.0 ) * electrNodesN[MM - 1]
					+ AAel / 2.0 / dx / dx * electrNodesN[MM - 2]
					+ 0.5 * ( CCel + CCel1 )
					+ AAel / dx * Jx * Jx * Rc / kEl * bb * hh + AAel / AA * kx / kEl
					* ( AA / 2.0 / dx / dx * plateNodesN[1] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * electrNodesN[MM - 1]
					+ AA / 2.0 / dx / dx * TleftG
					+ 0.5 * ( CC + CC1 ) );
	for( int i = MM; i < MM + NN - 2; ++i )
	{
		ci[i] = -AA / 2.0 / dx / dx;
		int j = i - MM + 1;
		di[i] = AA / 2.0 / dx / dx * plateNodesN[j + 1] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[j]
				+ AA / 2.0 / dx / dx * plateNodesN[j - 1]
				+ 0.5 * ( CC + CC1 );
	}
	di[MM + NN - 2] = AA / dx / dx * plateNodesN[NN - 2] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[NN - 1] + 0.5 * ( CC + CC1 ); 

	return 0;
}

int SolverThermoWElectrodes::doStep()
{
	curTime += dt;
	curTimeStep += 1;
	updateSystem();		//return system to the unmodified form

	ci[0] = ci[0] / bi[0];
	di[0] = di[0] / bi[0];
	for( int i = 1; i < MM + NN - 2; ++i )
	{
		ci[i] = ci[i] / ( bi[i] - ai[i] * ci[i - 1] );
		di[i] = ( di[i] - ai[i] * di[i - 1] ) / ( bi[i] - ai[i] * ci[i - 1] );
	}
	di[MM + NN - 2] = ( di[MM + NN - 2] - ai[MM + NN - 2] * di[MM + NN - 3] ) / ( bi[MM + NN - 2] - ai[MM + NN - 2] * ci[MM + NN - 3] );

	//now fill the solution
	plateNodesN1[NN - 1] = di[MM + NN - 2];
	for( int i = NN - 2; i > 0; --i )
	{
		int j = i + MM - 1;
		plateNodesN1[i] = di[j] - ci[j] * plateNodesN1[i + 1];
	}
	for( int i = MM - 1; i >= 0; --i )
	{
		electrNodesN1[i] = di[i] - ci[i] * electrNodesN1[i + 1];
	}
	plateNodesN1[0] = electrNodesN1[MM - 1];	//one of the boundary conditions

	//PL_NUM Jx = J0 * exp( -curTime / tauExp ) * sin( M_PI * curTime / tauSin );

	TelRightG1 = kx / kEl * 2.0 * plateNodesN1[1] + 2.0 * dx / kEl * Jx * Jx * Rc * bb * hh + electrNodesN1[MM - 2]
				- 2.0 * dx * dx / AA * kx / kEl * ( 1.0 / dt + AA / dx / dx + BB / 2.0 ) * plateNodesN1[0]
				+ 2.0 * dx * dx / AA * kx / kEl * ( AA / 2.0 / dx / dx * plateNodesN[1]
				+ ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[0] + AA / 2.0 / dx / dx * TleftG + CC );
	TleftG1 = plateNodesN1[1] + 2.0 * dx / kx * Jx * Jx * Rc * bb * hh - kEl / kx * TelRightG1 + kEl / kx * electrNodesN1[MM - 2];

	TelRightG = TelRightG1;
	TleftG = TleftG1;

	/*PL_NUM JxNext = J0 * exp( -( curTime + dt ) / tauExp ) * sin( M_PI * ( curTime + dt ) / tauSin );
	CC = CC1;
	CC1 = JxNext * JxNext / sigmaX / rho / cc;
	CCel = CCel1;
	CCel1 = JxNext * JxNext / sigmaEl / rhoEl / ccEl;*/

	for( int i = 0; i < MM; ++i )
	{
		electrNodesN[i] = electrNodesN1[i];
	}
	for( int i = 0; i < NN; ++i )
	{
		plateNodesN[i] = plateNodesN1[i];
	}

	return 0;
}

void SolverThermoWElectrodes::dumpSol()	//dump numerical soln + the soln obtained analytically for 1d case
{
	stringstream ss;
	ss << "./res/testSolThermo_" << curTimeStep << ".txt";
	ofstream of1( ss.str(), ofstream::app );
	for( int i = 0; i < MM; ++i )
	{
		of1 << electrNodesN1[i] + Tamb << endl;
	}
	//of1 << "----\n";
	for( int i = 0; i < NN; ++i )
	{
		of1 << plateNodesN1[i] + Tamb<< endl;
	}
	of1.close();
}

void SolverThermoWElectrodes::dumpSolBin()
{
	stringstream ss1;
	ss1 << "./res/testSolThermo_" << curTimeStep << ".bin";
	ofstream of1( ss1.str(), ofstream::out | ofstream::binary );

	stringstream ss2;
	ss2 << "./res/testSolThermoEl_" << curTimeStep << ".bin";
	ofstream of2( ss2.str(), ofstream::out | ofstream::binary );

	for( int i = 0; i < NN; ++i )
	{
		PL_NUM val = plateNodesN1[i] + Tamb;
		of1.write( reinterpret_cast<char*>( &( val ) ), sizeof( PL_NUM ) );
	}
	for( int i = 0; i < MM; ++i )
	{
		PL_NUM val = electrNodesN1[i] + Tamb;
		of2.write( reinterpret_cast<char*>( &( val ) ), sizeof( PL_NUM ) );
	}

	of1.close();
	of2.close();
}

//void SolverThermoWElectrodes::dumpSolElectrBin()
//{
//	stringstream ss;
//	ss << "./res/testSolThermoEl_" << curTimeStep << ".bin";
//	ofstream of1( ss.str(), ofstream::out | ofstream::binary );
//	for( int i = 0; i < MM; ++i )
//	{
//		PL_NUM val = electrNodesN1[i] + Tamb;
//		of1.write( reinterpret_cast<char*>( &( val ) ), sizeof( PL_NUM ) );
//	}
//}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

SolverThermoConstFlow::SolverThermoConstFlow( int _NN,
				PL_NUM _aa, PL_NUM _bb, PL_NUM _hh,
				PL_NUM _J0, PL_NUM _tauExp, PL_NUM _tauSin, PL_NUM _Tamb,
				PL_NUM _Rc,
				PL_NUM _kx, PL_NUM _hInf, PL_NUM _sigmaX, PL_NUM _rho, PL_NUM _cc,
				PL_NUM _dt )
{
	aa = _aa;
	bb = _bb;
	hh = _hh;
	J0 = _J0;
	tauExp = _tauExp;
	tauSin = _tauSin;
	Jx = 0.0;
	Tamb = _Tamb;
	Rc = _Rc;
	kx = _kx;
	hInf = _hInf;
	sigmaX = _sigmaX;
	rho = _rho;
	cc = _cc;
	NN = _NN;
	dt = _dt;

	//now, calculating the others
	dx = aa / 2.0 / ( NN - 1 );		//since we are solving for a half of the plate only

	plateNodesN = new PL_NUM[NN];
	plateNodesN1 = new PL_NUM[NN];

	for( int i = 0; i < NN; ++i )
	{
		plateNodesN[i] = 0.0;	//because we solve for dT
	}

	ai = new PL_NUM[NN];
	bi = new PL_NUM[NN];
	ci = new PL_NUM[NN];
	di = new PL_NUM[NN];

	AA = kx / rho / cc;
	BB = 2.0 * hInf / rho / cc / hh;
	CC = J0 * J0 / sigmaX / rho / cc; 

	curTime = 0.0;
	curTimeStep = 0;

	fillABSystem();	//TODO: re-check this when current is not const
}

SolverThermoConstFlow::~SolverThermoConstFlow()
{

}

int SolverThermoConstFlow::fillABSystem()
{
	bi[0] = ( 1.0 / dt + AA / dx / dx + BB / 2.0 );
	for( int i = 1; i < NN - 1; ++i )
	{
		ai[i] = -AA / 2.0 / dx / dx;
		bi[i] = 1.0 / dt + AA / dx / dx + BB / 2.0;
	}
	ai[NN - 1] = -AA / dx / dx;
	bi[NN - 1] = 1.0 / dt + AA / dx / dx + BB / 2.0;
	
	return 0;
}

int SolverThermoConstFlow::updateSystem()
{
	ci[0] = -AA / dx / dx;
	di[0] = AA / dx / dx * plateNodesN[1] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[0]
			+ CC + AA / 2.0 / dx / dx * Rc * bb * hh * dx / kx * 2.0 * J0 * J0;		//CAUTION: this is for const current only!
	for( int i = 1; i < NN - 1; ++i )
	{
		ci[i] = -AA / 2.0 / dx / dx;
		di[i] = AA / 2.0 / dx / dx * plateNodesN[i + 1] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[i]
				+ AA / 2.0 / dx / dx * plateNodesN[i - 1]
				+ CC;
	}
	di[NN - 1] = AA / dx / dx * plateNodesN[NN - 2] + ( 1.0 / dt - AA / dx / dx - BB / 2.0 ) * plateNodesN[NN - 1] + CC; 

	return 0;
}

int SolverThermoConstFlow::doStep()
{
	curTime += dt;
	curTimeStep += 1;
	updateSystem();		//return system to the unmodified form

	ci[0] = ci[0] / bi[0];
	di[0] = di[0] / bi[0];
	for( int i = 1; i < NN - 1; ++i )
	{
		ci[i] = ci[i] / ( bi[i] - ai[i] * ci[i - 1] );
		di[i] = ( di[i] - ai[i] * di[i - 1] ) / ( bi[i] - ai[i] * ci[i - 1] );
	}
	di[NN - 1] = ( di[NN - 1] - ai[NN - 1] * di[NN - 2] ) / ( bi[NN - 1] - ai[NN - 1] * ci[NN - 2] );

	//now fill the solution
	plateNodesN1[NN - 1] = di[NN - 1];
	for( int i = NN - 2; i >= 0; --i )
	{
		plateNodesN1[i] = di[i] - ci[i] * plateNodesN1[i + 1];
	}

	for( int i = 0; i < NN; ++i )
	{
		plateNodesN[i] = plateNodesN1[i];
	}

	return 0;
}

void SolverThermoConstFlow::dumpSol()	//dump numerical soln + the soln obtained analytically for 1d case
{
	stringstream ss;
	ss << "./res/testSolThermo_" << curTimeStep << ".txt";
	ofstream of1( ss.str(), ofstream::app );
	for( int i = 0; i < NN; ++i )
	{
		of1 << plateNodesN1[i] + Tamb<< endl;
	}
	of1.close();
}

void SolverThermoConstFlow::dumpSolBin()
{
	stringstream ss1;
	ss1 << "./res/testSolThermo_" << curTimeStep << ".bin";
	ofstream of1( ss1.str(), ofstream::out | ofstream::binary );

	for( int i = 0; i < NN; ++i )
	{
		PL_NUM val = plateNodesN1[i] + Tamb;
		of1.write( reinterpret_cast<char*>( &( val ) ), sizeof( PL_NUM ) );
	}

	of1.close();
}