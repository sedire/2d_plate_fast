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
	cout << " eigen lib version " << EIGEN_WORLD_VERSION << " " << EIGEN_MAJOR_VERSION << " " << EIGEN_MINOR_VERSION << endl;

	omp_set_num_threads( NUM_OF_THREADS );
	
	time_t beginTime;
	time_t endTime;

	Solver* solver = new Solver();

	/*solver->doStepThermo();
	solver->dumpThermo();
	for( int t = 0; t < 500; ++t )
	{
		cout << " ::::: " << t << endl;
		for( int i = 0; i < 1; ++i )
		{
			solver->doStepThermo();
		}
		solver->dumpThermo();
	}
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );*/

	time_t totTimeStart = time( 0 );

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
		solver->dump_whole_sol( 0 );
		solver->dump_whole_sol( 1 );
		solver->dump_whole_sol( 2 );
		//solver->dump_whole_sol( 3 );
		solver->dump_whole_sol( 4 );
		solver->dump_whole_sol( 5 );
		//solver->dump_whole_sol( 6 );
		//solver->dump_whole_sol( 7 );
		solver->dump_whole_sol( 8 );
		//solver->dump_whole_sol( 9 );
		//solver->dump_sol();
		//solver->dump_check_sol();
		//solver->dump_check_sol2D();
	}
	
	delete( solver );

	cout << ".........\n";
	cout << "... done!\n";
	cout << "... in " << time( 0 ) - totTimeStart << " sec\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}