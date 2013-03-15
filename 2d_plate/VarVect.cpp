#include "VarVect.h"

VarVect::VarVect() {}
VarVect::~VarVect() {}

VarVect::VarVect( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );
	Nk1.resize( _eq_num, 0.0 );
	d1N.resize( _eq_num, 0.0 );
	d2N.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}

void VarVect::setup( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );
	Nk1.resize( _eq_num, 0.0 );
	d1N.resize( _eq_num, 0.0 );
	d2N.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}