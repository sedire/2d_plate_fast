#include "Solver.h"

Solver::Solver():
	E1( 102970000000 ),
	//E2( 7550000000 ),
	E2( 102970000000 ),
	nu21( 0.3 ),
	nu23( 0.3 ),
	rho( 1594 ),
	G23( E2 / 2.0 / ( 1 + nu23 ) ),

	B11( E1 * E1 / ( E1 - nu21 * nu21 * E2 ) ),
	B22( E2 / ( 1 - nu21 * nu21 * E2 / E1 ) ),
	B12( nu21 * E2 * E1 / ( E1 - nu21 * nu21 * E2 ) ),
	B66( G23 ),

	By0( 0.0 ),
	By1( 2.0 * By0 ),
	By2( 0.0 ),

	betta( 0.25 ),

	mu( 0.00000125664 ),
	sigma_x( 39000 ),
	sigma_x_mu( sigma_x * mu ),
	sigma_y( sigma_x * 0.0001 ),
	//sigma_y( sigma_x ),
	sigma_y_mu( sigma_y * mu ),
	sigma_z( sigma_y ),

	//J0( 100000.0 ),
	J0( 0.0 ),
	omega( 0.0 ),
	tauC( 0.01 ),
	tauP( 0.01 ),
	p0( 1000000.0 ),
	impRadSq( 64.0 ),

	eps_0( 0.000000000008854 ),
	eps_x( 0.0000000002501502912 ),
	eps_x_0( eps_x - eps_0 ),

	hp( 0.0021 ),
	ap( 0.1524 ),		//len in x-dir
	bp( 0.1524 ),		//width in y-dir

	Km( NODES_ON_Y ),
	nx( NUMBER_OF_LINES ),
	eq_num( EQ_NUM ),
	varNum( nx * eq_num ),

	maxNewtonIt( MAX_NEWTON_IT ),
	newtonIt( 0 ),

	dx( ap / ( nx + 1 ) ),
	dy( bp / ( Km - 1 ) ),

	dt( 0.0001 ),
	cur_t( 0.0 ),
	curTimeStep( 0 ),

	al( 1.0 ),

	rungeKutta( 0 ),
	orthoBuilder( 0 )
{
	setTask();
}

Solver::~Solver()
{
	if( rungeKutta != 0 )
	{
		delete( rungeKutta );
	}
	if( orthoBuilder != 0 )
	{
		delete( orthoBuilder );
	}
}

PL_NUM Solver::increaseTime()
{
	cur_t += dt;
	++curTimeStep;
	return cur_t;
}

PL_NUM Solver::getCurTime()
{
	return cur_t;
}

void Solver::setTask()
{
	ofstream of1( "test_sol.txt" );
	of1.close();

	omega = _MMM_PI / tauC;
	//al = 1.0;
	//E1 = 102970000000;
	//E2 = E1;
	//E2 = 7550000000;
	//nu21 = 0.3;
	//nu23 = 0.3;
	//rho = 1594;
	//hp = 0.0021;
	//ap = 0.1524 * 100.0;	//len
	//bp = 0.1524;	//width
	//mu = 4 * _MMM_PI / 10000000;
	//mu = 0.00000125664;
	//sigma_x = 39000;
	//sigma_x_mu = sigma_x * mu;
	//sigma_y = sigma_x * 0.0001;
	//sigma_y_mu = sigma_y * mu;
	//sigma_z = sigma_y;

	//eq_num = EQ_NUM;
	//J0 = 0.0;
	//omega = 314.16;
	//p0 = 100.0;
	//By0 = 0.0;
	//eps_0 = 0.000000000008854;			// electric permittivity in the vacuum
	//eps_x = 0.0000000002501502912;		// electric permittivity in the fiber direction

	//Km = NODES_ON_Y;
	//nx = NUMBER_OF_LINES;
	//dt = 0.00005;
	//dx = ap / ( nx + 1 );
	//dy = bp / ( Km - 1 );
	//betta = 0.25;
	//varNum = nx * eq_num;

	rungeKutta = new RungeKutta( varNum );
	orthoBuilder = new OrthoBuilderGSh( varNum, Km );
	orthoBuilder->setParams();			//NOTE: it takes a lot of time to initialize so much memory

	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( varNum );
	}

	//matr_A.resize( varNum, vector<PL_NUM>( varNum, 0.0) );
	//cout << "A resized\n";
	/*vect_f.resize( varNum, 0.0 );
	cout << "f resized\n";*/

	for( int i = 0; i < nx * eq_num; ++i )
	{
		for( int j = 0; j < nx * eq_num; ++j )
		{
			matr_A[i][j] = 0.0;
		}
		vect_f[i] = 0.0;
	}

	//newmark_A.resize( varNum, 0.0 );
	//newmark_B.resize( varNum, 0.0 );
	//cout << "newmark resized\n";

	//B11 = E1 * E1 / ( E1 - nu21 * nu21 * E2 );
	//B22 = E2 / ( 1 - nu21 * nu21 * E2 / E1 );
	//B12 = nu21 * E2 * E1 / ( E1 - nu21 * nu21 * E2 );
	//G23 = E2 / 2.0 / ( 1 + nu23 );
	//B66 = G23;
	//By1 = 2.0 * By0;                                      // in considered boundary-value problem
	//By2 = 0.0;
	//eps_x_0 = eps_x - eps_0;
}

void Solver::calc_Newmark_AB( int _x, int mode )
{
	if( mode == 0 )
	{
		for( int i = 0; i < varNum; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk0[i] / betta / dt / dt - mesh[_x].d1N0[i] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N0[i];
			newmark_B[i] = -0.5 * mesh[_x].Nk0[i] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N0[i]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N0[i];
		}
	}
	else
	{
		for( int i = 0; i < varNum; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk1[i] / betta / dt / dt - mesh[_x].d1N[i] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N[i];
			newmark_B[i] = -0.5 * mesh[_x].Nk1[i] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N[i]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N[i];
		}
	}
}

void Solver::calc_system( int _x )
{
	PL_NUM h = hp;
	PL_NUM Btdt = 2 * dt * betta;
	PL_NUM Jx = J0 * exp( -( cur_t ) / tauC ) * sin( omega * ( cur_t ) );  
	PL_NUM Pimp = p0 * sin( 100.0 * _MMM_PI * ( cur_t ) );

	//strip load
	//PL_NUM cur_X = _x * dy - bp / 2.0;
	//if( fabs( cur_X ) <= h / 10.0 && cur_t + dt <= tauP )
	//{
	//	Pimp = p0 * sqrt( 1.0 - cur_X / h * 10.0 * cur_X / h * 10.0 ) * sin( _MMM_PI * ( cur_t + dt ) / tauP );
	//}

	PL_NUM Rad2 = ap * ap / impRadSq;

	int i = 0;
	int r = i + 1;
	int rr = i + 2;
	int j = i - 1;
	int jj = i - 2;
	
	//for the left line:
	matr_A[0 + i * eq_num][1 + r * eq_num] = -1.0 / ( 2.0 * dx ) / al;
	matr_A[0 + i * eq_num][2 + i * eq_num] = 1.0 / ( h  *B66 ) / al;

	matr_A[1 + i * eq_num][0 + r * eq_num] = - B12 / ( B22 * 2.0 * dx ) / al;
	matr_A[1 + i * eq_num][3 + i * eq_num] = 1.0 / ( h * B22 ) / al;

	matr_A[2 + i * eq_num][0 + r * eq_num] = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx ) / al;
	matr_A[2 + i * eq_num][0 + i * eq_num] = rho * h / ( betta * dt * dt ) - ( B12 * B12 / B22 - B11 ) * 2.0 * h / ( dx * dx ) + h * sigma_z * By1 * By1 / ( 4.0 * Btdt )
				+ ( B12 / B22 / 2.0 / dx * h * B12 / dx );
	matr_A[2 + i * eq_num][1 + r * eq_num] = eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
	matr_A[2 + i * eq_num][3 + r * eq_num] = -B12 / ( B22 * 2.0 * dx );
	matr_A[2 + i * eq_num][4 + r * eq_num] = -eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
	matr_A[2 + i * eq_num][8 + i * eq_num] = ( eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[1 + r * eq_num] ) 
				+newmark_B[1 + r * eq_num] ) - eps_x_0 * h * By1 / ( 4.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] ) 
				+newmark_B[4 + r * eq_num] ) );							//SIMPLIFY!
	matr_A[2 + i * eq_num][9 + i * eq_num] = ( h / ( mu * 2.0 * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - 1.0 * mesh[_x].Nk[9 + i * eq_num] ) + eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[8 + i * eq_num] 
				* ( 1.0 / Btdt * mesh[_x].Nk[1 + r * eq_num]
				+ newmark_B[1 + r * eq_num] ) )
				+ ( - h / mu / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] ) / al;
	matr_A[2 + i * eq_num][9 + r * eq_num] = ( h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );

	matr_A[3 + i * eq_num][0 + r * eq_num] = ( -B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][1 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / Btdt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) / al
				+ h * B66 / dx / dx / 2.0 / al;
	matr_A[3 + i * eq_num][2 + r * eq_num] = -1.0 / ( 2.0 * dx );
	matr_A[3 + i * eq_num][3 + i * eq_num] = ( eps_x_0 / ( B22 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][4 + i * eq_num] = ( -sigma_x * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] );
	matr_A[3 + i * eq_num][8 + i * eq_num] = ( mesh[_x].Nk[9 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * mesh[_x].Nk[0 + r * eq_num] + newmark_B[0 + r * eq_num] ) )
				- eps_x_0 * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) );
	matr_A[3 + i * eq_num][9 + i * eq_num] = ( mesh[_x].Nk[8 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] ) + newmark_B[0 + r * eq_num] ) )
				+ 2.0 * sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				- sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[4 + i * eq_num] + newmark_B[4 + i * eq_num] ) + h * Jx );

	matr_A[4 + i * eq_num][5 + i * eq_num] = 1.0 / al;

	matr_A[5 + i * eq_num][4 + r * eq_num] = -B12 / ( B22 * dx * dx ) / al;
	matr_A[5 + i * eq_num][4 + i * eq_num] = 2.0 * B12 / ( B22 * dx * dx ) / al;
	matr_A[5 + i * eq_num][6 + i * eq_num] = -12.0 / ( B22 * h * h * h ) / al;

	matr_A[6 + i * eq_num][4 + i * eq_num] = ( -h * h * h * B12 / B22 * eps_x_0 / ( 6.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] / Btdt );
	matr_A[6 + i * eq_num][4 + r * eq_num] = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] / Btdt );
	matr_A[6 + i * eq_num][5 + i * eq_num] = ( -rho * h * h * h / ( 6.0 * Btdt * dt ) - h * h * h * B66 / ( 3.0 * dx * dx )
				 - h * h * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( 12.0 * Btdt ) );
	matr_A[6 + i * eq_num][5 + r * eq_num] = h * h * h * B66 / ( 6.0 * dx * dx );
	matr_A[6 + i * eq_num][6 + i * eq_num] = ( eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( Btdt * B22 ) ) / al;
	matr_A[6 + i * eq_num][7 + i * eq_num] = 1.0 / al;
	matr_A[6 + i * eq_num][8 + i * eq_num] = ( eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * B12 * eps_x_0 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) 
				 + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] ) );
	matr_A[6 + i * eq_num][9 + i * eq_num] = ( -h * h * h / 6.0 * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 + eps_x_0 / B22 * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * eps_x_0 * B12 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) 
				 + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] ) );
	
	matr_A[7 + i * eq_num][1 + i * eq_num] = ( -sigma_x * h * By1 / Btdt / 2 * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[7 + i * eq_num][4 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / ( 4.0 * Btdt ) * ( By1 * By1 + 1.0 / 3.0 * By2 * By2 ) + rho * h * h * h / ( 3.0 * dx * dx * Btdt * dt )
				 + h * h * h / 6.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) )

				 +(h*h*h/3.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al 
					+(5.0*h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;

				 //+ ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al
				 //- ( h * h * h / 3.0 * ( B11 - B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + r * eq_num] = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 - h * h * h / ( 24.0 * dx * dx ) * sigma_y / Btdt * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] )		//TODO * N10i??
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) )

				  -(h*h*h/2.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al
				+(-4.0*h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;

				 //- ( h * h * h / 2.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al
				 //+ ( h * h * h / 2.0 * ( B11 - B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + rr * eq_num] = ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al

				+(h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;		
		
				//+ ( -h * h * h / 3.0 * ( B11 - B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + ( rr + 1 ) * eq_num] = 
		
				-(h*h*h/12.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al;

				//-( h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h / Btdt * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 + h * h * h / 6.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt ) );
	matr_A[7 + i * eq_num][5 + r * eq_num] = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 - h * h * h / 12.0 * eps_x_0 * ( ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num]) 
				 * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
	matr_A[7 + i * eq_num][6 + i * eq_num] = ( 2.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
	matr_A[7 + i * eq_num][6 + r * eq_num] = ( -( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
	matr_A[7 + i * eq_num][8 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 - eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - 1.0 * mesh[_x].Nk[9 + i * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] ) 
				 + newmark_B[5 + r * eq_num] )
				 - h * h * h * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] ) 
				 + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] ) )
				 + ( h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * mesh[_x].Nk[5 + r * eq_num]
				 + newmark_B[5 + r * eq_num] ) ) / al;			
	matr_A[7 + i * eq_num][8 + r * eq_num] = ( -h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] ) 
				 + newmark_B[5 + r * eq_num] ) );
	matr_A[7 + i * eq_num][9 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				 - eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[8 + r * eq_num] - 1.0 * mesh[_x].Nk[8 + i * eq_num] ) * ( 1.0 / Btdt * mesh[_x].Nk[5 + r * eq_num]
				 + newmark_B[5 + r * eq_num] )
				 - h * h * h * sigma_y / ( 24.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - 1.0 * mesh[_x].Nk[9 + i * eq_num] ) * ( 1.0 / Btdt * mesh[_x].Nk[4 + r * eq_num]
				 + newmark_B[4 + r * eq_num] )
				 - h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] / ( 6.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
				 + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] )
				 - h * h * h * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] / ( 12.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] )
				 + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] ) )
				 + ( h * h * h / 24.0 / dx * sigma_y / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * mesh[_x].Nk[4 + r * eq_num]
				 + newmark_B[4 + r * eq_num] ) + h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * mesh[_x].Nk[5 + r * eq_num]
				 + newmark_B[5 + r * eq_num] ) ) / al; 
	matr_A[7 + i * eq_num][9 + r * eq_num] = ( -h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] )
				 + newmark_B[4 + r * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * (mesh[_x].Nk[5 + r * eq_num] )
				 + newmark_B[5 + r * eq_num] ) );

	matr_A[8 + i * eq_num][0 + i * eq_num] = ( 1.0 / ( 2.0 * dx * Btdt ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] ) ) / al;
	matr_A[8 + i * eq_num][0 + r * eq_num] = ( 1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[8 + i * eq_num][9 + i * eq_num] = ( 1.0 / ( sigma_y_mu * dx * dx ) + 1.0 / Btdt + 1.0 / ( 2.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] )
				 + newmark_B[0 + r * eq_num] ) - 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );
	matr_A[8 + i * eq_num][9 + r * eq_num] = ( -1.0 / ( sigma_y_mu * dx * dx ) + 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) ) / al;

	matr_A[9 + i * eq_num][1 + i * eq_num] = sigma_x_mu / Btdt * mesh[_x].Nk[9 + i * eq_num];
	matr_A[9 + i * eq_num][4 + i * eq_num] = -sigma_x_mu * By1 / 4.0 / betta / dt;
	matr_A[9 + i * eq_num][8 + i * eq_num] = sigma_x_mu;
	matr_A[9 + i * eq_num][9 + i * eq_num] = sigma_x_mu * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] );

	vect_f[2 + i * eq_num] = h * sigma_z * By1 * By1 / 4.0 * newmark_B[0 + i * eq_num] 
			+ ( rho * h * newmark_A[0 + i * eq_num] - eps_x_0 * h / 4.0 / betta / dt / dx * ( mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + r * eq_num] 
			- By1 / 2.0 * mesh[_x].Nk[4 + r * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] - eps_x_0 * h / 2.0 / dx * mesh[_x].Nk[8 + i * eq_num]
			* ( 1.0 / Btdt * mesh[_x].Nk[1 + r * eq_num] + newmark_B[1 + r * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] )
			- h / dx / mu / 2.0 * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] );

	vect_f[3 + i * eq_num] = rho * h * newmark_A[1 + i * eq_num] - sigma_x * h * ( 1.0 / betta / dt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] ) 
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[4 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / Btdt / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			* mesh[_x].Nk[3 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] 
			* ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
			+ B12 / B22 * eps_x_0 * h / 4.0 / betta / dt / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[0 + r * eq_num]
			+ B12 / B22 * eps_x_0 * h / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[0 + r * eq_num]
			+ newmark_B[0 + r * eq_num] ) + eps_x_0 * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			- sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num];

	vect_f[6 + i * eq_num] = - rho * h * h * h / 12.0 * newmark_A[5 + i * eq_num] + sigma_x * h * h * h / 12.0 
			* ( 1.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] / Btdt  
			- eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[6 + i * eq_num] / Btdt + newmark_B[6 + i * eq_num] ) 
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx * ( 1.0 / Btdt
			* ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] );

	vect_f[7 + i * eq_num] = rho * h * newmark_A[4 + i * eq_num] + Pimp + sigma_x * h / 4.0 * By1 * By1 * newmark_B[4 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h / Btdt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) - h / 2.0 * By1 * Jx
			- rho * h * h * h / 12.0 / dx / dx * ( newmark_A[4 + r * eq_num] - 2.0 * newmark_A[4 + i * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / Btdt * mesh[_x].Nk[4 + r * eq_num]
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[4 + r * eq_num]
			+ newmark_B[4 + r * eq_num] ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] )
			+ h * h * h / 12.0 * sigma_y / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
			+ newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] / Btdt * mesh[_x].Nk[5 + r * eq_num]
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / Btdt * mesh[_x].Nk[5 + r * eq_num]
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + i * eq_num] ) * ( 1.0 / Btdt * mesh[_x].Nk[5 + r * eq_num]
			+ newmark_B[5 + r * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] ) * ( 1.0 / Btdt * mesh[_x].Nk[5 + r * eq_num]
			+ newmark_B[5 + r * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] ) )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] )
			+ newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] )
			- ( By1 * By1 * h * h * h * sigma_z / 48.0 * ( newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] ) ) / dx / dx / al;

	vect_f[8 + i * eq_num] = newmark_B[9 + i * eq_num] - mesh[_x].Nk[0 + i * eq_num] / 2.0 / dx / Btdt * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + i * eq_num] )
			- mesh[_x].Nk[9 + i * eq_num] / 2.0 / dx / Btdt * ( mesh[_x].Nk[0 + r * eq_num] );

	vect_f[9 + i * eq_num] = By2 / h - sigma_x_mu / Btdt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - 0.5 * sigma_x_mu * By1 * newmark_B[4 + i * eq_num];

	i = nx - 1;
	r = i + 1;
	rr = i + 2;
	j = i - 1;
	jj = i - 2;
	
	//for the right line:
	matr_A[0 + i * eq_num][1 + j * eq_num] = 1.0 / ( 2.0 * dx ) / al;
	matr_A[0 + i * eq_num][2 + i * eq_num] = 1.0 / ( h * B66 ) / al;

	matr_A[1 + i * eq_num][0 + j * eq_num] = B12 / ( B22 * 2.0 * dx ) / al;
	matr_A[1 + i * eq_num][3 + i * eq_num] = 1.0 / ( h * B22 ) / al;

	matr_A[2 + i * eq_num][0 + j * eq_num] = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx );
	matr_A[2 + i * eq_num][0 + i * eq_num] = rho * h * 2.0 / ( Btdt * dt ) - ( B12 * B12 / B22 - B11 ) * 2.0 * h / ( dx * dx ) + h * sigma_z * By1 * By1 / ( 4.0 * Btdt )
				+ ( B12 / B22 / 2.0 / dx * h * B12 / dx );		//SHOULD BE -?? - NO, LOOKS LIKE A MISTAKE IN AMIR'S CODE
	matr_A[2 + i * eq_num][1 + j * eq_num] = -eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
	matr_A[2 + i * eq_num][3 + j * eq_num] = B12 / ( B22 * 2.0 * dx ) / al;
	matr_A[2 + i * eq_num][4 + j * eq_num] = eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
	matr_A[2 + i * eq_num][8 + i * eq_num] = ( eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( -mesh[_x].Nk[1 + j * eq_num] ) 
				-newmark_B[1 + j * eq_num] ) - eps_x_0 * h * By1 / ( 4.0 * dx ) * ( 1.0 / Btdt * ( -mesh[_x].Nk[4 + j * eq_num] ) 
				-newmark_B[4 + j * eq_num] ) );
	matr_A[2 + i * eq_num][9 + i * eq_num] = ( h / ( mu * 2.0 * dx ) * ( 1.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) + eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[8 + i * eq_num] 
				* ( -1.0 / Btdt * mesh[_x].Nk[1 + j * eq_num]
				- newmark_B[1 + j * eq_num] ) ) + ( h / mu / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] ) / al;
	matr_A[2 + i * eq_num][9 + j * eq_num] = ( -h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );

	matr_A[3 + i * eq_num][0 + j * eq_num] = ( B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][1 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / Btdt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) / al
				+ h * B66 / 2.0 / dx / dx / al;			//SHOULD BE -?? - NO, LOOKS LIKE A MISTAKE IN AMIR'S CODE
	matr_A[3 + i * eq_num][2 + j * eq_num] = 1.0 / ( 2.0 * dx );
	matr_A[3 + i * eq_num][3 + i * eq_num] = ( eps_x_0 / ( B22 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][4 + i * eq_num] = ( -sigma_x * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[3 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] );
	matr_A[3 + i * eq_num][8 + i * eq_num] = ( mesh[_x].Nk[9 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( -1.0 / Btdt * mesh[_x].Nk[0 + j * eq_num] - newmark_B[0 + j * eq_num] ) )
				- eps_x_0 * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) );
	matr_A[3 + i * eq_num][9 + i * eq_num] = ( mesh[_x].Nk[8 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( -mesh[_x].Nk[0 + j * eq_num] ) - newmark_B[0 + j * eq_num] ) )
				+ 2.0 * sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				- sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[4 + i * eq_num] + newmark_B[4 + i * eq_num] ) + h * Jx );

	matr_A[4 + i * eq_num][5 + i * eq_num] = 1.0 / al;

	matr_A[5 + i * eq_num][4 + j * eq_num] = -B12 / ( B22 * dx * dx ) / al;
	matr_A[5 + i * eq_num][4 + i * eq_num] = 2.0 * B12 / ( B22 * dx * dx ) / al;
	matr_A[5 + i * eq_num][6 + i * eq_num] = -12.0 / ( B22 * h * h * h ) / al;

	matr_A[6 + i * eq_num][4 + j * eq_num] = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
	matr_A[6 + i * eq_num][4 + i * eq_num] = ( -h * h * h * B12 / B22 * eps_x_0 / ( 6.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
	matr_A[6 + i * eq_num][5 + i * eq_num] = ( -rho * h * h * h / ( 6.0 * Btdt * dt ) - 2.0 * h * h * h * B66 / ( 6.0 * dx * dx )
				 - h * h * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( 12.0 * Btdt ) );
	matr_A[6 + i * eq_num][5 + j * eq_num] = h * h * h * B66 / ( 6.0 * dx * dx );
	matr_A[6 + i * eq_num][6 + i * eq_num] = ( eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( Btdt * B22 ) );
	matr_A[6 + i * eq_num][7 + i * eq_num] = 1.0 / al;
	matr_A[6 + i * eq_num][8 + i * eq_num] = ( eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * B12 * eps_x_0 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) 
				 + newmark_B[4 + j * eq_num] - 2.0 * newmark_B[4 + i * eq_num] ) );
	matr_A[6 + i * eq_num][9 + i * eq_num] = ( -h * h * h / 6.0 * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 + eps_x_0 / B22 * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * eps_x_0 * B12 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) 
				 + newmark_B[4 + j * eq_num] - 2.0 * newmark_B[4 + i * eq_num] ) );
	
	matr_A[7 + i * eq_num][1 + i * eq_num] = ( -sigma_x * h * By1 / Btdt / 2.0 * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[7 + i * eq_num][4 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / ( 4.0 * Btdt ) * ( By1 * By1 + 1.0 / 3.0 * By2 * By2 ) + rho * h * h * h / ( 3.0 * dx * dx * Btdt * dt )
				 + h * h * h / 6.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) )
				 
				 +(h*h*h/3.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al
             +(5.0*h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;
				 
				 //+ ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al
				 //- ( h * h * h / 3.0 * ( B11 - B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + j * eq_num] = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 + h * h * h / ( 24.0 * dx * dx ) * sigma_y / Btdt * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) )
				 
				 -(h*h*h/2.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al
             +(-4.0*h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;

				 //- ( h * h * h / 2.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al
				 //+ ( h * h * h / 2.0 * ( B11 - B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + jj * eq_num] = 
		
		+(h*h*h/3.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al
             +(h*h*h/12.0*(B11-B12*B12/B22)/dx/dx/dx/dx)/al;
		
		//( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 ) / dx / dx / dx / dx ) / al
			//	 + ( h * h * h / 3.0 * ( -B11 + B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;			//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][4 + ( jj - 1 ) * eq_num] = 
				-(h*h*h/12.0*(2.0*B66*B12/B22)/dx/dx/dx/dx)/al ;
		
		//- ( h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / dx / dx / dx / dx ) / al;		//CAUTION: ANOTHER THING IN AMIR'S VERSION
	matr_A[7 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h / Btdt * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 + h * h * h / 6.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt ) );
	matr_A[7 + i * eq_num][5 + j * eq_num] = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 + h * h * h / 12.0 * eps_x_0 * ( ( mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num]) 
				 * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
	matr_A[7 + i * eq_num][6 + i * eq_num] = ( 2.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
	matr_A[7 + i * eq_num][6 + j * eq_num] = ( -1.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
	matr_A[7 + i * eq_num][8 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 - eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( 1.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] ) 
				 - newmark_B[5 + j * eq_num] )
				 - h * h * h * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + j * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] ) 
				 + newmark_B[5 + j * eq_num] - 2.0 * newmark_B[5 + i * eq_num] ) )
				 + ( -h * h * h / 48.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] ) - newmark_B[5 + j * eq_num] ) );
	matr_A[7 + i * eq_num][8 + j * eq_num] = ( h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] ) 
				 - newmark_B[5 + j * eq_num] ) );
	matr_A[7 + i * eq_num][9 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				 - eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( 1.0 * mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( -1.0 / Btdt * mesh[_x].Nk[5 + j * eq_num]
				 - newmark_B[5 + j * eq_num] )
				 - h * h * h * sigma_y / ( 24.0 * dx * dx ) * ( 1.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( -1.0 / Btdt * mesh[_x].Nk[4 + j * eq_num]
				 - newmark_B[4 + j * eq_num] )
				 - h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] / ( 6.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
				 + newmark_B[4 + j * eq_num] - 2.0 * newmark_B[4 + i * eq_num] )
				 - h * h * h * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] / ( 12.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + j * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] )
				 + newmark_B[5 + j * eq_num] - 2.0 * newmark_B[5 + i * eq_num] ) )
				 + ( -h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / betta / dt 
				 * ( -mesh[_x].Nk[4 + j * eq_num] ) - newmark_B[4 + j * eq_num] )
				 - h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / betta / dt 
				 * ( -mesh[_x].Nk[5 + j * eq_num] ) - newmark_B[5 + j * eq_num] ) ) / al;
	matr_A[7 + i * eq_num][9 + j * eq_num] = ( h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( -mesh[_x].Nk[4 + j * eq_num] )
				 - newmark_B[4 + j * eq_num] )
				 + h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] )
				 - newmark_B[5 + j * eq_num] ) );

	matr_A[8 + i * eq_num][0 + i * eq_num] = ( 1.0 / ( 2.0 * dx * Btdt ) * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) );
	matr_A[8 + i * eq_num][0 + j * eq_num] = ( -1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
	matr_A[8 + i * eq_num][9 + i * eq_num] = ( 1.0 / ( sigma_y_mu * dx * dx ) + 1 / Btdt + 1.0 / ( 2.0 * dx ) * ( 1.0 / Btdt * ( -mesh[_x].Nk[0 + j * eq_num] )
				 - newmark_B[0 + j * eq_num] ) 
				 + 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );
	matr_A[8 + i * eq_num][9 + j * eq_num] = ( -1.0 / ( sigma_y_mu * dx * dx ) - 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );

	matr_A[9 + i * eq_num][1 + i * eq_num] = sigma_x_mu / Btdt * mesh[_x].Nk[9 + i * eq_num];
	matr_A[9 + i * eq_num][4 + i * eq_num] = -sigma_x_mu * By1 / 4.0 / betta / dt;
	matr_A[9 + i * eq_num][8 + i * eq_num] = sigma_x_mu;
	matr_A[9 + i * eq_num][9 + i * eq_num] = sigma_x_mu * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] );

	vect_f[2 + i * eq_num] = h * sigma_z * By1 * By1 / 4.0 * newmark_B[0 + i * eq_num] 
			+ ( rho * h * newmark_A[0 + i * eq_num] - eps_x_0 * h / 4.0 / betta / dt / dx * ( -mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + j * eq_num] 
			+ By1 / 2.0 * mesh[_x].Nk[4 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] - eps_x_0 * h / 2.0 / dx * mesh[_x].Nk[8 + i * eq_num]
			* ( -1.0 / Btdt * mesh[_x].Nk[1 + j * eq_num] - newmark_B[1 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] )
			- h / dx / mu / 2.0 * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] );

	vect_f[3 + i * eq_num] = rho * h * newmark_A[1 + i * eq_num] - sigma_x * h * ( 1.0 / betta / dt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] ) 
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[4 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / Btdt / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			* mesh[_x].Nk[3 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] 
			* ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
			+ B12 / B22 * eps_x_0 * h / 4.0 / betta / dt / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( -mesh[_x].Nk[0 + j * eq_num] )
			+ B12 / B22 * eps_x_0 * h / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( -1.0 / Btdt * mesh[_x].Nk[0 + j * eq_num]
			- newmark_B[0 + j * eq_num] ) + eps_x_0 * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			- sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num];

	vect_f[6 + i * eq_num] = - rho * h * h * h / 12.0 * newmark_A[5 + i * eq_num] + sigma_x * h * h * h / 12.0 
			* ( 1.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] / Btdt  
			- eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[6 + i * eq_num] / Btdt + newmark_B[6 + i * eq_num] ) 
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx / Btdt * ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx * ( 1.0 / Btdt
			* ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] ) + newmark_B[4 + j * eq_num] - 2.0 * newmark_B[4 + i * eq_num] );

	vect_f[7 + i * eq_num] = rho * h * newmark_A[4 + i * eq_num] + Pimp + sigma_x * h / 4.0 * By1 * By1 * newmark_B[4 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h / Btdt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) - h / 2.0 * By1 * Jx
			- rho * h * h * h / 12.0 / dx / dx * ( newmark_A[4 + j * eq_num] - 2.0 * newmark_A[4 + i * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / Btdt * ( -mesh[_x].Nk[4 + j * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * ( -1.0 / Btdt * mesh[_x].Nk[4 + j * eq_num]
			- newmark_B[4 + j * eq_num] ) * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
			+ h * h * h / 12.0 * sigma_y / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 / Btdt * ( mesh[_x].Nk[4 + j * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] )
			+ newmark_B[4 + j * eq_num] - 2.0 * newmark_B[4 + i * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / Btdt * ( -mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( -1.0 / Btdt * mesh[_x].Nk[5 + j * eq_num]
			- newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( -1.0 / Btdt * mesh[_x].Nk[5 + j * eq_num]
			- newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + j * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] ) )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + j * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] )
			+ newmark_B[5 + j * eq_num] - 2.0 * newmark_B[5 + i * eq_num] )
			- ( By1 * By1 * h * h * h / 48.0 * sigma_z * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) ) / dx / dx / al;

	vect_f[8 + i * eq_num] = newmark_B[9 + i * eq_num] - mesh[_x].Nk[0 + i * eq_num] / 2.0 / dx / Btdt * ( mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
			- mesh[_x].Nk[9 + i * eq_num] / 2.0 / dx / Btdt * ( -mesh[_x].Nk[0 + j * eq_num] );

	vect_f[9 + i * eq_num] = By2 / h - sigma_x_mu / Btdt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - 0.5 * sigma_x_mu * By1 * newmark_B[4 + i * eq_num];

	for( int i = 1; i < nx - 1; ++i )
	{
		//Pimp = p0 * sin( 100.0 * _MMM_PI * ( cur_t ) );
		//PL_NUM rad2 = ( ( Km - 1 ) / 2 - _x ) * dy * ( ( Km - 1 ) / 2 - _x ) * dy + ( ( nx - 1 ) / 2 - i ) * dx * ( ( nx - 1 ) / 2 - i ) * dx;
		//if( rad2 < Rad2 )
		//{
		//	Pimp = p0 * sqrt( 1 - rad2 / Rad2 ) * sin( 100.0 * _MMM_PI * ( cur_t ) );
		//}
		//else if( _x == ( Km - 1 ) / 2 )
		//{
		//	cout << " == line " << i << " is out\n";
		//}

		r = i + 1;
		rr = i + 2;
		j = i - 1;
		jj = i - 2;

		matr_A[0 + i * eq_num][1 + r * eq_num] = -1.0 / ( 2.0 * dx );
		matr_A[0 + i * eq_num][1 + j * eq_num] = 1.0 / ( 2.0 * dx );
		matr_A[0 + i * eq_num][2 + i * eq_num] = 1.0 / ( h * B66 );

		matr_A[1 + i * eq_num][0 + r * eq_num] = - B12 / ( B22 * 2.0 * dx );
		matr_A[1 + i * eq_num][0 + j * eq_num] = B12 / ( B22 * 2.0 * dx );
		matr_A[1 + i * eq_num][3 + i * eq_num] = 1.0 / ( h * B22 );

		matr_A[2 + i * eq_num][0 + r * eq_num] = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx );
		matr_A[2 + i * eq_num][0 + j * eq_num] = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx );
		matr_A[2 + i * eq_num][0 + i * eq_num] = rho * h * 2.0 / ( Btdt * dt ) - ( B12 * B12 / B22 - B11 ) * 2.0 * h / ( dx * dx ) + h * sigma_z * By1 * By1 / ( 4.0 * Btdt );
		matr_A[2 + i * eq_num][1 + r * eq_num] = eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
		matr_A[2 + i * eq_num][1 + j * eq_num] = -eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
		matr_A[2 + i * eq_num][3 + r * eq_num] = -B12 / ( B22 * 2.0 * dx );
		matr_A[2 + i * eq_num][3 + j * eq_num] = B12 / ( B22 * 2.0 * dx );
		matr_A[2 + i * eq_num][4 + r * eq_num] = -eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
		matr_A[2 + i * eq_num][4 + j * eq_num] = eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
		matr_A[2 + i * eq_num][8 + i * eq_num] = ( eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
				+newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] ) - eps_x_0 * h * By1 / ( 4.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) 
				+newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) );
		matr_A[2 + i * eq_num][9 + i * eq_num] = ( h / ( mu * 2.0 * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) + eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[8 + i * eq_num] 
				* ( 1.0 / Btdt * ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
				+ newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] ) );
		matr_A[2 + i * eq_num][9 + r * eq_num] = ( h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[2 + i * eq_num][9 + j * eq_num] = ( -h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );

		matr_A[3 + i * eq_num][0 + r * eq_num] = ( -B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[3 + i * eq_num][0 + j * eq_num] = ( B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[3 + i * eq_num][1 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / Btdt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[3 + i * eq_num][2 + r * eq_num] = -1.0 / ( 2.0 * dx );
		matr_A[3 + i * eq_num][2 + j * eq_num] = 1.0 / ( 2.0 * dx );
		matr_A[3 + i * eq_num][3 + i * eq_num] = ( eps_x_0 / ( B22 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[3 + i * eq_num][4 + i * eq_num] = ( -sigma_x * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[3 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] );
		matr_A[3 + i * eq_num][8 + i * eq_num] = ( mesh[_x].Nk[9 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) )
				- eps_x_0 * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) );
		matr_A[3 + i * eq_num][9 + i * eq_num] = ( mesh[_x].Nk[8 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) )
				+ 2.0 * sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				- sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[4 + i * eq_num] + newmark_B[4 + i * eq_num] ) + h * Jx );

		matr_A[4 + i * eq_num][5 + i * eq_num] = 1.0;

		matr_A[5 + i * eq_num][4 + r * eq_num] = -B12 / ( B22 * dx * dx );
		matr_A[5 + i * eq_num][4 + i * eq_num] = 2.0 * B12 / ( B22 * dx * dx );
		matr_A[5 + i * eq_num][4 + j * eq_num] = -B12 / ( B22 * dx * dx );
		matr_A[5 + i * eq_num][6 + i * eq_num] = -12.0 / ( B22 * h * h * h );

		matr_A[6 + i * eq_num][4 + i * eq_num] = ( -h * h * h * B12 / B22 * eps_x_0 / ( 6.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matr_A[6 + i * eq_num][4 + r * eq_num] = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matr_A[6 + i * eq_num][4 + j * eq_num] = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matr_A[6 + i * eq_num][5 + i * eq_num] = ( -rho * h * h * h / ( 6.0 * Btdt * dt ) - 2.0 * h * h * h * B66 / ( 6.0 * dx * dx )
				 - h * h * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( 12.0 * Btdt ) );
		matr_A[6 + i * eq_num][5 + r * eq_num] = h * h * h * B66 / ( 6.0 * dx * dx );
		matr_A[6 + i * eq_num][5 + j * eq_num] = h * h * h * B66 / ( 6.0 * dx * dx );
		matr_A[6 + i * eq_num][6 + i * eq_num] = ( eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( Btdt * B22 ) );
		matr_A[6 + i * eq_num][7 + i * eq_num] = 1.0;
		matr_A[6 + i * eq_num][8 + i * eq_num] = ( eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * B12 * eps_x_0 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) );
		matr_A[6 + i * eq_num][9 + i * eq_num] = ( -h * h * h / 6.0 * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 + eps_x_0 / B22 * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * eps_x_0 * B12 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) );

		matr_A[7 + i * eq_num][1 + i * eq_num] = ( -sigma_x * h * By1 / Btdt / 2.0 * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[7 + i * eq_num][4 + i * eq_num] = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / ( 4.0 * Btdt ) * ( By1 * By1 + 1.0 / 3.0 * By2 * By2 ) + rho * h * h * h / ( 3.0 * dx * dx * Btdt * dt )
				 + h * h * h / 6.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );
		matr_A[7 + i * eq_num][4 + r * eq_num] = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 - h * h * h / ( 24.0 * dx * dx ) * sigma_y / Btdt * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );
		matr_A[7 + i * eq_num][4 + j * eq_num] = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 + h * h * h / ( 24.0 * dx * dx * Btdt ) * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );              
		matr_A[7 + i * eq_num][4 + i * eq_num] = matr_A[7 + i * eq_num][4 + i * eq_num] - ( h * h * h / 2.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		matr_A[7 + i * eq_num][4 + r * eq_num] = matr_A[7 + i * eq_num][4 + r * eq_num] + ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		matr_A[7 + i * eq_num][4 + j * eq_num] = matr_A[7 + i * eq_num][4 + j * eq_num] + ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		if( i != nx - 2 )
		{
			matr_A[7 + i * eq_num][4 + rr * eq_num] = ( -h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		}
		if( i != 1 )
		{
			matr_A[7 + i * eq_num][4 + jj * eq_num] = ( -h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );        
		}
		matr_A[7 + i * eq_num][5 + i * eq_num] = ( -eps_x_0 * h / Btdt * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 + h * h * h / 6.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt ) );
		matr_A[7 + i * eq_num][5 + r * eq_num] = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 - h * h * h / 12.0 * eps_x_0 * ( (mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + r * eq_num]
				 - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
		matr_A[7 + i * eq_num][5 + j * eq_num] = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 + h * h * h / 12.0 * eps_x_0 * ( ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + r * eq_num]
				 - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
		matr_A[7 + i * eq_num][6 + i * eq_num] = ( 2.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matr_A[7 + i * eq_num][6 + r * eq_num] = ( -1.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matr_A[7 + i * eq_num][6 + j * eq_num] = ( -1.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matr_A[7 + i * eq_num][8 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 - eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
				 - h * h * h * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num]
				 + mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) );
		matr_A[7 + i * eq_num][8 + r * eq_num] = ( -h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );
		matr_A[7 + i * eq_num][8 + j * eq_num] = ( h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );           
		matr_A[7 + i * eq_num][9 + i * eq_num] = ( -sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				 - eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
				 - h * h * h * sigma_y / ( 24.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] )
				 - h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] / ( 6.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] )
				 - h * h * h * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] / ( 12.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num]
				 + mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) );
		matr_A[7 + i * eq_num][9 + r * eq_num] = ( -h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * (mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );
		matr_A[7 + i * eq_num][9 + j * eq_num] = ( h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * (mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) 
				 + h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) ); 

		matr_A[8 + i * eq_num][0 + i * eq_num] = ( 1.0 / ( 2.0 * dx * Btdt ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) );
		matr_A[8 + i * eq_num][0 + r * eq_num] = ( 1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[8 + i * eq_num][0 + j * eq_num] = ( -1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matr_A[8 + i * eq_num][9 + i * eq_num] = ( 2.0 / ( sigma_y_mu * dx * dx ) + 1 / Btdt + 1.0 / ( 2.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] )
				 + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) );
		matr_A[8 + i * eq_num][9 + r * eq_num] = ( -1.0 / ( sigma_y_mu * dx * dx ) + 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );
		matr_A[8 + i * eq_num][9 + j * eq_num] = ( -1.0 / ( sigma_y_mu * dx * dx ) - 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );

		matr_A[9 + i * eq_num][1 + i * eq_num] = sigma_x_mu / Btdt * mesh[_x].Nk[9 + i * eq_num];
		matr_A[9 + i * eq_num][4 + i * eq_num] = -sigma_x_mu * By1 / 4.0 / betta / dt;
		matr_A[9 + i * eq_num][8 + i * eq_num] = sigma_x_mu;
		matr_A[9 + i * eq_num][9 + i * eq_num] = sigma_x_mu * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] );

		vect_f[2 + i * eq_num] = ( - ( ( 2 * (B11 - B12 * B12 / B22) * h) / dx / dx + ( h * rho ) / ( betta * dt * dt ) + ( By1 * By1 * h * sigma_z ) / ( 8 * betta * dt ) ) ) 
			* mesh[_x].Nk[0 + i * eq_num] +	h * rho * ( newmark_A[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( betta * dt * dt) ) + ( 1 / 4 ) * By1 * By1 * h * sigma_z 
			* ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2 * betta * dt ) ) + ( ( B11 - B12 * B12 / B22 ) * h * mesh[_x].Nk[0 + r * eq_num] ) / dx / dx 
			+ ( ( B11 - B12 * B12 / B22 ) * h * mesh[_x].Nk[0 + j * eq_num] ) / dx / dx - ( ( B11 - B12 * B12 / B22 ) * h * ( -2 * mesh[_x].Nk[0 + i * eq_num] + mesh[_x].Nk[0 + r * eq_num] 
			+ mesh[_x].Nk[0 + j * eq_num] ) ) / dx / dx - ( h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] ) / ( 2 * dx * mu ) 
			+ ( h * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2 * dx * mu) 
			+ ( h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] ) / ( 2 * dx * mu ) + ( B12 * mesh[_x].Nk[3 + r * eq_num] ) / ( 2 * B22 * dx )
			- ( B12 * ( mesh[_x].Nk[3 + r * eq_num] - mesh[_x].Nk[3 + j * eq_num] ) ) / ( 2 * B22 * dx ) - ( B12 * mesh[_x].Nk[3 + j * eq_num] ) / ( 2 * B22 * dx )
			- ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * betta * dt * dx ) + ( eps_x_0 * h 
			* mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num]) / ( 2 * dx ) + ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * betta * dt * dx) 
			+ ( By1 * eps_x_0 * h * mesh[_x].Nk[4 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 8 * betta * dt * dx ) - ( ( eps_x_0 
			* h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
			/ ( 2 * betta * dt ) ) ) / ( 2 * dx ) - ( By1 * eps_x_0 * h * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] 
			+ ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 4 * dx ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- ( By1 * eps_x_0 * h * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] + ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * dx ) - ( By1 * eps_x_0 * h * mesh[_x].Nk[4 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 8 * betta * dt * dx )
			- mesh[_x].Nk[9 + i * eq_num] * ( ( h * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2 * dx * mu ) + ( eps_x_0 * h * ( newmark_B[1 + r * eq_num] 
			- newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * dx ) );

		vect_f[3 + i * eq_num] = h * Jx * mesh[_x].Nk[9 + i * eq_num] - ( ( h * rho ) / ( betta * dt * dt ) + ( h *sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] )
			/ ( 2 * betta * dt ) ) * mesh[_x].Nk[1 + i * eq_num] + h * rho * ( newmark_A[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( betta * dt * dt ) ) 
			+ h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2 * betta * dt ) ) 
			+ mesh[_x].Nk[2 + r * eq_num] / ( 2 * dx ) - ( mesh[_x].Nk[2 + r * eq_num] - mesh[_x].Nk[2 + j * eq_num] ) / ( 2 * dx ) - mesh[_x].Nk[2 + j * eq_num]
			/ ( 2 * dx ) + ( By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num] ) / ( 4 * betta * dt ) - ( 1 / 2 ) * By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num]
			* ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2 * betta * dt ) ) + h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 
			( B12 * eps_x_0 * h * mesh[_x].Nk[0 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * betta * B22 * dt * dx ) 
			- ( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[9 + i * eq_num] 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * B22 * dx ) - ( B12 * eps_x_0 * h * mesh[_x].Nk[0 + j * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) 
			/ ( 4 * betta * B22 * dt * dx ) - ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[3 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * betta * B22 * dt ) 
			+ ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * (newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / B22 
			+ ( By1 * eps_x_0 * h * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * betta * dt ) - ( 1 / 2 ) * By1 * eps_x_0 * h 
			* ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] - ( h * sigma_x 
			* mesh[_x].Nk[9 + i * eq_num] - ( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[9 + i * eq_num] ) / ( 2 * B22 * dx ) + ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2 * betta * dt ) ) )
			/ B22 - ( 1 / 2 ) * By1 * eps_x_0 * h * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- mesh[_x].Nk[9 + i * eq_num] * ( h * Jx + 2 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2 * betta * dt ) )
			- ( 1 / 2 ) * By1 * h * sigma_x * ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2 * betta * dt ) ) + h * sigma_x * mesh[_x].Nk[8 + i * eq_num] -
			( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * B22 * dx ) + ( eps_x_0 * ( newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / B22 );

		vect_f[6 + i * eq_num] = /*( - ( - ( ( B66 * h * h * h ) / ( 3 * dx * dx ) ) - ( h * h * h * rho ) / ( 12 * betta * dt * dt ) 
			- ( h * h * h* sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) / ( 24 * betta * dt ) ) ) * mesh[_x].Nk[5 + i * eq_num] - ( 1 / 12 ) 
			* h * h * h * rho * ( newmark_A[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( betta * dt * dt ) ) - ( 1 / 12 ) * h * h * h *
			sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) 
			- ( B66 * h * h * h * mesh[_x].Nk[5 + r * eq_num] ) / ( 6 * dx * dx ) - ( B66 * h * h * h * mesh[_x].Nk[5 + j * eq_num] ) / ( 6 * dx * dx ) 
			+ ( B66 * h * h * h * ( -2 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + r * eq_num] + mesh[_x].Nk[5 + j * eq_num] ) ) / ( 6 * dx * dx ) + ( B12 * eps_x_0 *
			h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 12 * betta * B22 * dt * dx * dx ) 
			- ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num]) / ( 24 * betta * B22 * dt * dx * dx ) 
			- ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num]) / ( 24 * betta * B22 * dt 
			* dx * dx ) + ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] + newmark_B[4 + j * eq_num] 
			+ ( -2 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 12 * B22 * dx * dx ) 
			- ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * betta * B22 * dt ) 
			+ ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) 
			/ B22 - ( ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] + newmark_B[4 + j * eq_num] 
			+ ( -2 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) 
			/ ( 12 * B22 * dx * dx ) + ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] / ( 2 * betta * dt ) ) ) / B22 ) 
			* mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[9 + i * eq_num] * ( ( -( 1 / 6 ) ) * h * h * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] 
			+ mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) + ( B12 * eps_x_0 * h * h * h * ( -2 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] 
			+ newmark_B[4 + j * eq_num] + ( -2 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) 
			/ ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 12 * B22 * dx * dx ) + ( eps_x_0 * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] 
			/ ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / B22 );*/
			( -rho * h * h * h / 12.0 * newmark_A[5 + i * eq_num] + sigma_x * h * h * h / 12.0 * ( 1.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] / 2.0 / betta / dt
			- eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[6 + i * eq_num] / 2.0 / betta / dt + newmark_B[6 + i * eq_num] )
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx / 2.0 / betta / dt
			* ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] )
			- eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx * ( 1.0 / 2.0 / betta / dt
			* ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] )
			+ newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) ) / al;

		vect_f[7 + i * eq_num] = /*( -( 1 / 2 ) ) * By1 * h * Jx + Pimp + ( By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] ) / ( 4 * betta * dt ) 
			- ( 1 / 2 ) * By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2 * betta * dt ) ) 
			- ( -( ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h ) / ( 2 * dx * dx * dx * dx ) ) + ( h * rho ) / ( betta * dt * dt ) 
			+ ( h * h * h * rho ) / ( 6 * betta * dt * dt * dx * dx ) + ( By1 * By1 * h * sigma_x ) / ( 8 * betta * dt ) + ( h * h * h * ( ( By1 * By1 * sigma_z ) / 4 
			+ sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) ) / ( 12 * betta * dt * dx * dx ) ) * mesh[_x].Nk[4 + i * eq_num] + h * rho * (newmark_A[4 + i * eq_num] 
			+ mesh[_x].Nk[4 + i * eq_num] / ( betta * dt * dt ) ) + ( 1 / 4 ) * By1 * By1 * h *	sigma_x * ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2 * betta * dt ) ) 
			- ( ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h ) / ( 3 * dx * dx * dx * dx ) - ( h * h * h * rho ) / ( 12 * betta * dt * dt * dx * dx ) 
			- ( h * h * h * ( ( By1 * By1 * sigma_z ) / 4 + sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) ) / ( 24 * betta * dt * dx * dx )
			- ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 48 * betta * dt * dx * dx ) ) 
			* mesh[_x].Nk[4 + r * eq_num] - ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] * (mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) 
			* ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] + ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 24 * dx * dx ) 
			- ( ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h ) / ( 3 * dx * dx * dx * dx ) - ( h * h * h * rho ) / ( 12 * betta * dt * dt * dx * dx ) 
			- ( h * h * h * ( ( By1 * By1 * sigma_z ) / 4 + sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) ) / ( 24 * betta * dt * dx * dx ) 
			+ ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 48 * betta * dt * dx * dx ) )     
			* mesh[_x].Nk[4 + j * eq_num] - ( h * h * h * rho * ( -2 * newmark_A[4 + i * eq_num] + newmark_A[4 + r * eq_num] + newmark_A[4 + j * eq_num] 
			+ ( -2 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) / ( betta * dt * dt ) ) ) 
			/ ( 12 * dx * dx ) - ( h * h * h * ( ( By1 * By1 * sigma_z ) / 4 + sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) 
			* ( -2 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] + newmark_B[4 + j * eq_num] + ( -2 * mesh[_x].Nk[4 + i * eq_num] 
			+ mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 12 * dx * dx ) 
			+ ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h * ( i == nx - 2 ? 0 : mesh[_x].Nk[4 + rr * eq_num] ) ) / ( 12 * dx * dx * dx * dx )      
			+ ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h * ( i == 1 ? 0 : mesh[_x].Nk[4 + jj * eq_num] ) ) / ( 12 * dx * dx * dx * dx )
			- ( ( -B11 + B12 * B12 / B22 + ( 2 * B12 * B66 ) / B22 ) * h * h * h * ( 6 * mesh[_x].Nk[4 + i * eq_num] 
			- 4 * mesh[_x].Nk[4 + r * eq_num] - 4 * mesh[_x].Nk[4 + j * eq_num] + ( i == nx - 2 ? 0 :  mesh[_x].Nk[4 + rr * eq_num] ) + ( i == 1 ? 0 : mesh[_x].Nk[4 + jj * eq_num] ) ) ) 
			/ ( 12 * dx * dx * dx * dx ) - ( 2 * ( B12 / B22 + ( 2 * B66 ) / B22 ) * mesh[_x].Nk[6 + i * eq_num] ) / dx * dx  + ( ( B12 / B22 + ( 2 * B66 ) / B22 ) 
			* mesh[_x].Nk[6 + r * eq_num] ) / dx * dx + ( ( B12 / B22 + ( 2 * B66 ) / B22 ) * mesh[_x].Nk[6 + j * eq_num] ) / dx / dx - ( ( B12 / B22 + ( 2 * B66 ) / B22 ) 
			* ( -2 * mesh[_x].Nk[6 + i * eq_num] + mesh[_x].Nk[6 + r * eq_num] + mesh[_x].Nk[6 + j * eq_num] ) ) / dx / dx - ( 1 / 2 ) * By1 * h * sigma_x * mesh[_x].Nk[8 + i * eq_num] 
			- eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2 * newmark_B[5 + i * eq_num] + newmark_B[5 + r * eq_num] 
			+ newmark_B[5 + j * eq_num] + ( -2 * mesh[_x].Nk[5 + i * eq_num]    
			+ mesh[_x].Nk[5 + r * eq_num] + mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 12 * dx * dx ) 
			- ( ( -( 1 / 2 ) ) * By1 * h * sigma_x - eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] 
			/ ( 2 * betta * dt ) ) - ( eps_x_0 * h * h * h * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( newmark_B[5 + r * eq_num] 
			- newmark_B[5 + j * eq_num] + ( mesh[_x].Nk[5 + r * eq_num] 
			- mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 48 * dx * dx ) - ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2 * newmark_B[5 + i * eq_num] 
			+ newmark_B[5 + r * eq_num] + newmark_B[5 + j * eq_num] + ( -2 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + r * eq_num] 
			+ mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) ) 
			/ ( 12 * dx * dx ) ) * mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[5 + i * eq_num] * ( -( ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * betta * dt ) ) + ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) 
			/ ( 12 * betta * dt * dx * dx ) ) - mesh[_x].Nk[9 + r * eq_num] * ( -( ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] 
			* ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] 
			+ ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 24 * dx * dx ) ) - ( eps_x_0 * h * h * h 
			* ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] + ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) 
			/ ( 48 * dx * dx ) ) - mesh[_x].Nk[9 + j * eq_num] * ( ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] 
			+ ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 24 * dx * dx ) + ( eps_x_0 * h * h * h 
			* ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] + (mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 48 * dx * dx ) ) + ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] 
			* ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] + ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) 
			/ ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + r * eq_num] ) / ( 48 * dx * dx ) - mesh[_x].Nk[5 + r * eq_num] * ( -( ( eps_x_0 * h * h * h 
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 24 * betta * dt * dx * dx ) ) - ( eps_x_0 * h * h * h * ( ( ( mesh[_x].Nk[9 + r * eq_num] 
			- mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * dx ) + ( mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) ) 
			/ ( 2 * dx ) ) ) / ( 48 * betta * dt * dx ) ) - mesh[_x].Nk[5 + j * eq_num] * ( -( ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) 
			/ ( 24 * betta * dt * dx * dx ) ) + ( eps_x_0 * h * h * h * ( ( ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num]) / ( 2 * dx ) 
			+ ( mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) ) / ( 2 * dx ) ) ) / ( 48 * betta * dt * dx ) ) - ( eps_x_0 
			* h * h * h * ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] + ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) 
			* ( ( ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 2 * dx ) + ( mesh[_x].Nk[9 + i * eq_num] 
			* ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) ) / ( 2 * dx ) ) ) / ( 24 * dx ) - mesh[_x].Nk[9 + i * eq_num] 
			* ( ( -( 1 / 2 ) ) * By1 * h * sigma_x * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2 * betta * dt ) ) 
			- ( h * h * h * sigma_y * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] + ( mesh[_x].Nk[4 + r * eq_num] 
			- mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) / ( 24 * dx * dx ) - ( h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( -2 * newmark_B[4 + i * eq_num] 
			+ newmark_B[4 + r * eq_num] + newmark_B[4 + j * eq_num] + ( -2 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * betta * dt ) ) ) 
			/ ( 6 * dx * dx ) - eps_x_0 * h * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- ( eps_x_0 * h * h * h * ( -2 * newmark_B[5 + i * eq_num] + newmark_B[5 + r * eq_num] + newmark_B[5 + j * eq_num] + ( -2 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + r * eq_num] 
			+ mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 12 * dx * dx ) - ( eps_x_0 * h * h * h 
			* ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] + ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) * (mesh[_x].Nk[8 + r * eq_num] 
			- mesh[_x].Nk[8 + j * eq_num] ) ) / ( 48 * dx * dx ) ) - ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] 
			+ ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] ) / ( 2 * betta * dt ) ) * mesh[_x].Nk[8 + j * eq_num] ) / ( 48 * dx * dx );*/
			( rho * h * newmark_A[4 + i * eq_num] + Pimp + sigma_x * h / 4.0 * By1 * By1 * newmark_B[4 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / betta / dt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h / 2.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) - h / 2.0 * By1 * Jx
			- rho * h * h * h / 12.0 / dx / dx * ( newmark_A[4 + r * eq_num] - 2.0 * newmark_A[4 + i * eq_num] + newmark_A[4 + j * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / 2.0 / betta / dt * ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * ( mesh[_x].Nk[4 + r * eq_num]
			- mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
			+ h * h * h / 12.0 * sigma_y / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 / 2.0 / betta / dt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
			+ mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( 1.0 / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + j * eq_num] ) )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / betta / dt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) ) / al
			- ( By1 * By1 * h * h * h * sigma_z / 48.0 * ( newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) ) / dx / dx / al;

		vect_f[8 + i * eq_num] = newmark_B[9 + i * eq_num] + mesh[_x].Nk[9 + i * eq_num] / ( 2 * betta * dt ) - ( mesh[_x].Nk[0 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) 
			/ ( 4 * betta * dt * dx ) - ( 1 / ( 2 * betta * dt ) + 2 / ( dx * dx * mu * sigma_y ) + ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] 
			+ ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] )	/ ( 2 * betta * dt ) ) / ( 2 * dx ) ) * mesh[_x].Nk[9 + i * eq_num] 
			+ ( ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) 
			/ ( 2 * betta * dt ) ) * mesh[_x].Nk[9 + i * eq_num] ) / ( 2 * dx ) + ( mesh[_x].Nk[0 + j * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) / ( 4 * betta * dt * dx ) 
			- ( -( 1 / ( dx * dx * mu * sigma_y ) ) + ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2 * betta * dt ) ) / ( 2 * dx ) ) * mesh[_x].Nk[9 + r * eq_num] 
			- ( mesh[_x].Nk[0 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 4 * betta * dt * dx ) + ( ( newmark_B[0 + i * eq_num] 
			+ mesh[_x].Nk[0 + i * eq_num] / ( 2 * betta * dt ) ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2 * dx ) - ( -( 1 / ( dx * dx * mu * sigma_y ) ) 
			- ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2 * betta * dt ) ) / ( 2 * dx ) ) * mesh[_x].Nk[9 + j * eq_num] 
			- ( -2 * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) / ( dx * dx * mu * sigma_y );

		vect_f[9 + i * eq_num] = -( ( mu * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] ) / ( 2 * betta * dt ) ) + ( By1 * mu * sigma_x 
			* mesh[_x].Nk[4 + i * eq_num] ) / ( 4 * betta * dt ) - ( 1 / 2 ) * By1 * mu * sigma_x * ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2 * betta * dt ) );
	}
}

void Solver::walkthrough( int mode )
{
	time_t tBeg, tEnd;
	time_t rgkT = 0;
	time_t orthoT = 0;

	vector<PL_NUM> baseVect;
	baseVect.resize( varNum, 0.0 );

	//integrate and orthonorm
	int _x = 0;
#pragma omp parallel firstprivate( baseVect )
	{
	for( _x; _x < Km - 1; )
	{
		if( omp_get_thread_num() == 0 )
		{
			calc_Newmark_AB( _x, mode );
			calc_system( _x );
		}

		int begIt = omp_get_thread_num() * ( varNum / 2 / NUM_OF_THREADS + 1 );
		int endIt = begIt + varNum / 2 / NUM_OF_THREADS + 1;
		if( omp_get_thread_num() == NUM_OF_THREADS - 1 )
		{
			endIt = varNum / 2;
		}

		#pragma omp for
		for( int vNum = 0; vNum < varNum / 2; ++vNum )
		{
			//rungeKutta->calc( matr_A, vect_f, dy, omp_get_thread_num(), 0, orthoBuilder->zi[_x][vNum], &( orthoBuilder->zi[_x + 1][vNum] ) );
			rungeKutta->calc( matr_A, vect_f, dy, omp_get_thread_num(), 0, orthoBuilder->zi[_x][vNum], decompVect[vNum] );
		}
		#pragma omp barrier

		if( omp_get_thread_num() == 0 )
		{
			for( int vNum = 0; vNum < varNum / 2; ++vNum )
			{
				orthoBuilder->orthonorm( vNum, _x, decompVect[vNum] );
			}

			//rungeKutta->calc( matr_A, vect_f, dy, omp_get_thread_num(), 1, orthoBuilder->z5[_x], &( orthoBuilder->z5[_x + 1] ) );
			rungeKutta->calc( matr_A, vect_f, dy, omp_get_thread_num(), 1, orthoBuilder->z5[_x], decompVect[varNum / 2] );
			orthoBuilder->orthonorm( varNum / 2, _x, decompVect[varNum / 2] );
			++_x;
		}
		#pragma omp barrier		//may be this is redundant
	}
	}

	orthoBuilder->buildSolution( &mesh );

	for( int _x = 0; _x < Km; ++_x )
	{
		orthoBuilder->flushO( _x );
	}
}

void Solver::updateDerivs()
{
	for( int i = 0; i < Km; ++i )
	{
		for( int j = 0; j < varNum; ++j )
		{
			mesh[i].d2N[j] = ( mesh[i].Nk1[j] - mesh[i].Nk0[j] ) / betta / dt / dt - mesh[i].d1N0[j] / betta / dt - ( 0.5 - betta ) / betta * mesh[i].d2N0[j];
			mesh[i].d1N[j] = mesh[i].d1N0[j] + 0.5 * dt * ( mesh[i].d2N0[j] + mesh[i].d2N[j] );
			mesh[i].Nk0[j] = mesh[i].Nk1[j];
			mesh[i].d1N0[j] = mesh[i].d1N[j];
			mesh[i].d2N0[j] = mesh[i].d2N[j];
		}
	}
}

void Solver::pre_step()
{
	for( int i = 0; i < nx; ++i )			//TODO check indexes
	{
		orthoBuilder->zi[0][i * eq_num / 2 + 0][i * eq_num + 2] = 1.0;
		orthoBuilder->zi[0][i * eq_num / 2 + 1][i * eq_num + 3] = 1.0;
		orthoBuilder->zi[0][i * eq_num / 2 + 2][i * eq_num + 5] = 1.0;
		orthoBuilder->zi[0][i * eq_num / 2 + 3][i * eq_num + 7] = 1.0;
		orthoBuilder->zi[0][( i + 1 ) * eq_num / 2 - 1][( i + 1 ) * eq_num - 1] = -1.0;
	}
	//z5s are already zeros

	calc_Newmark_AB( 0, 0 );
	calc_system( 0 );

	walkthrough( 0 );
}


void Solver::do_step()
{	
	int cont = 1;
	while( cont == 1 )
	{
		cout << " = walk\n";
		calc_Newmark_AB( 0, 1 );
		for( int i = 0; i < varNum; ++i ) 			//TODO check indexes
		{
			for( int j = 0; j < varNum / 2; ++j )
			{
				orthoBuilder->zi[0][j][i] = 0;			//TODO lags here
			}
			orthoBuilder->z5[0][i] = 0;
		}
		for( int i = 0; i < nx; ++i )			//TODO check indexes
		{
			orthoBuilder->zi[0][i * eq_num / 2 + 0][i * eq_num + 2] = 1.0;
			orthoBuilder->zi[0][i * eq_num / 2 + 1][i * eq_num + 3] = 1.0;
			orthoBuilder->zi[0][i * eq_num / 2 + 2][i * eq_num + 5] = 1.0;
			orthoBuilder->zi[0][i * eq_num / 2 + 3][i * eq_num + 7] = 1.0;
			orthoBuilder->zi[0][( i + 1 ) * eq_num / 2 - 1][( i + 1 ) * eq_num - 2] = mesh[0].Nk1[i * eq_num + 1] / betta / 2.0 / dt + newmark_B[i * eq_num + 1];
			orthoBuilder->zi[0][( i + 1 ) * eq_num / 2 - 1][( i + 1 ) * eq_num - 1] = -1.0;

			//orthoBuilder->solInfoMap[0].z5[i * eq_num + 8] = newmark_B[i * eq_num + 1] * mesh[0].Nk1[i * eq_num + 9] - newmark_B[i * eq_num + 4] * By0;
			//orthoBuilder->solInfoMap[0].z5[i * eq_num + 9] = -mesh[0].Nk1[i * eq_num + 9];
			orthoBuilder->z5[0][i * eq_num + 8] = newmark_B[i * eq_num + 1] * mesh[0].Nk1[i * eq_num + 9] - newmark_B[i * eq_num + 4] * By0;
			orthoBuilder->z5[0][i * eq_num + 9] = -mesh[0].Nk1[i * eq_num + 9];
		}

		for( int i = 0; i < Km; ++i )
		{
			for( int j = 0; j < varNum; ++j )
			{
				mesh[i].Nk[j] = mesh[i].Nk1[j];
			}
		}
		walkthrough( 0 );
		cont = checkConv();
	}
	updateDerivs();
}

int Solver::checkConv()
{
	//if( newtonIt >= maxNewtonIt )		//stopping criterion with a fixed number of iterations
	//{
	//	newtonIt = 0;
	//	return 0;
	//}
	//else
	//{
	//	int i = 15;
	//	cout << " newton iteration " << newtonIt << endl;
	//	cout << " divergence in " << Km / 2 << " " << i << " " << fabsl( ( mesh[Km / 2].Nk1[i] - mesh[Km / 2].Nk[i] ) / mesh[Km / 2].Nk[i] ) << endl;
	//	++newtonIt;
	//	return 1;
	//}

	for( int x = 0; x < Km; ++x ) //old and weird stopping criterion. I think it works only because 1-2 iterations are almost always suffice. otherwise this doed not make sense to me
	{
		for( int i = 0; i < varNum; ++i )
		{
			if( mesh[x].Nk[i] != 0.0 )
			{
				if( fabsl( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) < ALMOST_ZERO )
				{
					cout << " divergence " << x << " " << i << " " << fabsl( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) << " delta is " << ALMOST_ZERO << endl;
					return 0;
				}
			}
			else
			{
				if( fabsl( mesh[x].Nk1[i] ) < ALMOST_ZERO )
				{
					cout << " divergence " << x << " " << i << " " << fabsl( mesh[x].Nk1[i] ) << " delta is " << ALMOST_ZERO << endl;
					return 0;
				}
			}
		}
	}
	return 1;

	//for( int x = 0; x < Km; ++x ) //new stopping criterion: the max of the absolute value of the relative difference + just the abs value check
	//{
	//	for( int i = 0; i < varNum; ++i )
	//	{
	//		if( fabsl( mesh[x].Nk[i] ) > ALMOST_ZERO * 1000.0 )
	//		{
	//			if( fabsl( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) > QUASILIN_CHECK && fabsl( mesh[x].Nk1[i] - mesh[x].Nk[i] ) > QUASILIN_CHECK )
	//			{
	//				cout << " :: divergence " << x << " " << i << " " << fabsl( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) << " delta is " << QUASILIN_CHECK << endl;
	//				cout << " :: " << mesh[x].Nk1[i] << " " << mesh[x].Nk[i] << endl;
	//				return 1;
	//			}
	//		}
	//		else
	//		{
	//			if( fabsl( mesh[x].Nk1[i] ) > ALMOST_ZERO * 1000.0 )
	//			{
	//				cout << " :: divergence -- 0 -- " << x << " " << i << " " << fabsl( mesh[x].Nk1[i] ) << " delta is " << QUASILIN_CHECK << endl;
	//				return 1;
	//			}
	//		}
	//	}
	//}
	//return 0;
}

void Solver::dump_sol()
{
	ofstream dumpSol;
	stringstream ss;
	ss << "sol_" << cur_t << ".txt";
	
	dumpSol.open ( ss.str() );
	
	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < varNum; ++i )
		{
			dumpSol << mesh[x].Nk1[i] << " ";
		}
		dumpSol << endl;
	}


	dumpSol.close();
	return;
}

void Solver::dump_check_sol()	//dump numerical soln + the soln obtained analytically for 1d case
{
	PL_NUM sum = 0.0;
	PL_NUM h = hp;
	PL_NUM a = bp;
	PL_NUM t = cur_t - dt;

	int minusOne = -1;

	/*for( int i = 0; i <= 1000000; ++i )
	{
		PL_NUM omg = _MMM_PI * _MMM_PI * ( 2 * i + 1 ) * ( 2 * i + 1 ) * h / 2 / a / a * sqrtl( B22 / 3 / rho );

		minusOne = -minusOne;

		sum = sum + minusOne / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) * cosl( omg * t );
	}*/
	PL_NUM wTheor;
	wTheor = - p0 * a * a * a * a / h / h / h / B22 * ( 5.0 / 32.0 - 48.0 / _MMM_PI / _MMM_PI / _MMM_PI / _MMM_PI / _MMM_PI * sum );
	wTheor = 1.0;

	ofstream of1( "test_sol.txt", ofstream::app );
	of1 << t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] << " ; " << wTheor << " ; " << fabs( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] ) / wTheor ) << endl;
	of1.close();
}

void Solver::dump_check_sol2D()		//dump numerical soln + the soln obtained analytically for 2d case. WARNING: the plate should be square-shaped and isotropic !!!
{
	PL_NUM sum = 0.0;
	PL_NUM sum0 = 0.0;
	PL_NUM h = hp;
	PL_NUM a = bp;
	PL_NUM t = cur_t - dt;

	PL_NUM DD = E1 * h * h * h / 12 / ( 1 - nu21 * nu21 );
	PL_NUM Om = 100.0 * _MMM_PI;

	for( int m = 1; m < 50; ++m )
	{
		for( int n = 1; n < 50; ++n )
		{
			PL_NUM wmn = _MMM_PI * _MMM_PI * ( m * m / ap / ap + n * n / bp / bp ) * sqrt( DD / rho / h );
			PL_NUM Wmn = -16.0 * p0 / _MMM_PI / _MMM_PI / m / n / rho / h / ( wmn * wmn - Om * Om );
			sum += ( sin( Om * t ) - Om / wmn * sin( wmn * t ) ) * Wmn * sin( m * _MMM_PI / 2.0 ) * sin( n * _MMM_PI / 2.0 );
			sum0 = sum;
		}
	}
	PL_NUM wTheor = sum;

	ofstream of1( "test_sol.txt", ofstream::app );
	of1 << t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] << " ; " << wTheor << " ; " << fabsl( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] ) / wTheor ) << endl;
	of1.close();
}

void Solver::dump_border_vals()
{
	ofstream of1( "sol_borders.txt", ofstream::app );

	of1 << "\t" << cur_t << endl;
	of1 << " ~~ top\n";
	of1 << "u : " << mesh[Km - 1].Nk1[0 + (nx-1)/2 * eq_num] << " v : " << mesh[Km - 1].Nk1[1 + (nx-1)/2 * eq_num] << " Nxy : " << mesh[Km - 1].Nk1[2 + (nx-1)/2 * eq_num]
			<< " Nyy : " << mesh[Km - 1].Nk1[3 + (nx-1)/2 * eq_num] << " w : " << mesh[Km - 1].Nk1[4 + (nx-1)/2 * eq_num] << " W : " << mesh[Km - 1].Nk1[5 + (nx-1)/2 * eq_num]
			<< " Myy : " << mesh[Km - 1].Nk1[6 + (nx-1)/2 * eq_num] << " Nyz : " << mesh[Km - 1].Nk1[7 + (nx-1)/2 * eq_num] << " Ex : " << mesh[Km - 1].Nk1[8 + (nx-1)/2 * eq_num]
			<< " Bz : " << mesh[Km - 1].Nk1[9 + (nx-1)/2 * eq_num] << endl;
	of1 << " ~~ bottom\n";
	of1 << "u : " << mesh[0].Nk1[0 + (nx-1)/2 * eq_num] << " v : " << mesh[0].Nk1[1 + (nx-1)/2 * eq_num] << " Nxy : " << mesh[0].Nk1[2 + (nx-1)/2 * eq_num]
			<< " Nyy : " << mesh[0].Nk1[3 + (nx-1)/2 * eq_num] << " w : " << mesh[0].Nk1[4 + (nx-1)/2 * eq_num] << " W : " << mesh[0].Nk1[5 + (nx-1)/2 * eq_num]
			<< " Myy : " << mesh[0].Nk1[6 + (nx-1)/2 * eq_num] << " Nyz : " << mesh[0].Nk1[7 + (nx-1)/2 * eq_num] << " Ex : " << mesh[0].Nk1[8 + (nx-1)/2 * eq_num]
			<< " Bz : " << mesh[0].Nk1[9 + (nx-1)/2 * eq_num] << endl;
	of1 << " ~~ left\n";
	of1 << "u : " << mesh[( Km - 1 ) / 2].Nk1[0 + 0 * eq_num] << " v : " << mesh[( Km - 1 ) / 2].Nk1[1 + 0 * eq_num] << " Nxy : " << mesh[( Km - 1 ) / 2].Nk1[2 + 0 * eq_num]
			<< " Nyy : " << mesh[( Km - 1 ) / 2].Nk1[3 + 0 * eq_num] << " w : " << mesh[( Km - 1 ) / 2].Nk1[4 + 0 * eq_num] << " W : " << mesh[( Km - 1 ) / 2].Nk1[5 + 0 * eq_num]
			<< " Myy : " << mesh[( Km - 1 ) / 2].Nk1[6 + 0 * eq_num] << " Nyz : " << mesh[( Km - 1 ) / 2].Nk1[7 + 0 * eq_num] << " Ex : " << mesh[( Km - 1 ) / 2].Nk1[8 + 0 * eq_num]
			<< " Bz : " << mesh[( Km - 1 ) / 2].Nk1[9 + 0 * eq_num] << endl;
	of1 << " ~~ right\n";
	of1 << "u : " << mesh[( Km - 1 ) / 2].Nk1[0 + (nx-1) * eq_num] << " v : " << mesh[( Km - 1 ) / 2].Nk1[1 + (nx-1) * eq_num] << " Nxy : " << mesh[( Km - 1 ) / 2].Nk1[2 + (nx-1) * eq_num]
			<< " Nyy : " << mesh[( Km - 1 ) / 2].Nk1[3 + (nx-1) * eq_num] << " w : " << mesh[( Km - 1 ) / 2].Nk1[4 + (nx-1) * eq_num] << " W : " << mesh[( Km - 1 ) / 2].Nk1[5 + (nx-1) * eq_num]
			<< " Myy : " << mesh[( Km - 1 ) / 2].Nk1[6 + (nx-1) * eq_num] << " Nyz : " << mesh[( Km - 1 ) / 2].Nk1[7 + (nx-1) * eq_num] << " Ex : " << mesh[( Km - 1 ) / 2].Nk1[8 + (nx-1) * eq_num]
			<< " Bz : " << mesh[( Km - 1 ) / 2].Nk1[9 + (nx-1) * eq_num] << endl;
	of1 << "====\n";
	of1.close();
}

void Solver::dump_whole_sol( int var )
{
	stringstream ss;
	ss << "sol_whole_" << curTimeStep << ".bin";
	ofstream of1( ss.str(), ofstream::out | ofstream::binary );
	for( int y = 0; y < Km; ++y )
	{
		for( int x = 0; x < nx; ++x )
		{
			of1.write( reinterpret_cast<char*>( &(mesh[y].Nk1[var + x * eq_num]) ), sizeof(PL_NUM) );
		}
	}
	of1.close();
}

void Solver::dump_Amir_sol()
{
	stringstream ss;
	ss << "sol_Amir_" << curTimeStep << ".txt";
	ofstream of1( ss.str() );

	const int y( Km / 2 );

	for( int x = 0; x < nx; ++x )
	{
		of1 << mesh[y].Nk1[3 + x * eq_num] << endl;
	}

	of1.close();
}

void Solver::dumpMatrA( int _x )
{
	stringstream ss;
	ss << "matrA_" << curTimeStep << "_" << _x << "_" << nx << "_" << ap << ".txt";

	ofstream of( ss.str() );

	for( int i = 0; i < EQ_NUM * NUMBER_OF_LINES; ++i )
	{
		for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
		{
			of << matr_A[i][j] << " ";
		}
		of << endl;
	}
	of << endl;
	for( int i = 0; i < EQ_NUM * NUMBER_OF_LINES; ++i )
	{
		of << vect_f[i] << endl;
	}

	of << hp << " " << B66 << " " << B12 << " " << B22 << " " << B11 << " " << dx << endl;

	of.close();
}