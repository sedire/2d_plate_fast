#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include <omp.h>
#include <time.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	//omp_set_num_threads( NUM_OF_THREADS );
	
	time_t beginTime;
	time_t endTime;

	//const int size1( 4 );
	//const int size2( 3 );

	//Matrix<double, size1, size2> qrA = Matrix<double, size1, size2>::Random();
	//ColPivHouseholderQR<Matrix<double, size1, size2> > qr;
	//Matrix<double, size1, size1> qrQ;
	//Matrix<double, size1, size2> qrR;
	//Matrix<double, Dynamic, Dynamic> qrR2;

	//qr.compute( qrA );
	//qrQ = qr.householderQ();
	//qrR2 = qr.matrixQR();//.triangularView<Upper>();
	//qrR = qr.matrixQR().triangularView<Upper>();

	//cout << qrA << endl;
	//cout << " ------------------\n";
	//cout << qrQ << endl;
	//cout << " ------------------\n";
	//cout << qrR << endl;
	//cout << " ------------------\n";
	//cout << qrR2 << endl;
	//cout << " ------------------\n";
	////cout << qrQ * qrR << endl;
	//cout << " ------------------\n";
	//cout << ( qrA - qrQ * qrR ).lpNorm<Infinity>() << endl;
	//cout << " ------------------\n";

	Solver* solver = new Solver();
	//solver->setTask();

	cout << "\n doing pre_step...\n";
	beginTime = time( 0 );
	solver->pre_step();
	solver->dump_sol();
	endTime = time( 0 );

	cout << "\n pre_step done\n";
	cout << "pre_step done in " << float( endTime - beginTime ) << " ~~\n";

	while( solver->getCurTime() <= 0.05 )
	{
		for( int i = 0; i < 1; ++i )
		{
			beginTime = time( 0 );
			solver->do_step();
			endTime = time( 0 );
			cout << "step done in " << float( endTime - beginTime ) << " ~~\n";

			solver->increaseTime();
			cout << solver->getCurTime() << " -- step done\n";
		}
		//solver->dump_whole_sol( 0 );
		//solver->dump_whole_sol( 1 );
		//solver->dump_whole_sol( 2 );
		//solver->dump_whole_sol( 3 );
		//solver->dump_whole_sol( 4 );
		//solver->dump_whole_sol( 5 );
		//solver->dump_whole_sol( 6 );
		//solver->dump_whole_sol( 7 );
		//solver->dump_whole_sol( 8 );
		//solver->dump_whole_sol( 9 );
		//solver->dump_sol();
		//solver->dump_check_sol();
		solver->dump_check_sol2D();
	}
	
	delete( solver );

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}