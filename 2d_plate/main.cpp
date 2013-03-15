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

	cout << "initializing...\n";

	time_t beginTime = time( 0 );
	time_t initTime;
	time_t endTime;

	Solver* solver = new Solver();
	solver->setTask();

	initTime = time( 0 );
	cout << "initialization done in " << float( initTime - beginTime ) / 60.0 << " min\n";

	cout << "\n doing pre_step...\n";
	beginTime = time( 0 );
	solver->pre_step();
	endTime = time( 0 );
	
	cout << "hellp\n";

	solver->cur_t += solver->dt;
	++( solver->curTimeStep );
	solver->dump_check_sol2D();

	cout << "\n pre_step done\n";
	cout << "pre_step done in " << float( endTime - beginTime ) << " ~~\n";

	while( solver->cur_t <= 0.016 )
	{
		for( int i = 0; i < 1; ++i )
		{
			solver->do_step();

			solver->cur_t += solver->dt;
			++( solver->curTimeStep );
			cout << solver->cur_t << " -- step done\n";
		}
		solver->dump_check_sol();
		//solver->dump_left_border_vals();
	}
	
	free( solver );

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}