#include "VarVect.h"

VarVect::VarVect() {}
VarVect::~VarVect() {}

VarVect::VarVect( int _varNum )
{
	Nk.resize( _varNum, 0.0 );
	Nk1.resize( _varNum, 0.0 );
	d1N.resize( _varNum, 0.0 );
	d2N.resize( _varNum, 0.0 );

	Nk0.resize( _varNum, 0.0 );
	d1N0.resize( _varNum, 0.0 );
	d2N0.resize( _varNum, 0.0 );
}

void VarVect::setup( int _varNum )
{
	Nk.resize( _varNum, 0.0 );
	Nk1.resize( _varNum, 0.0 );
	d1N.resize( _varNum, 0.0 );
	d2N.resize( _varNum, 0.0 );

	Nk0.resize( _varNum, 0.0 );
	d1N0.resize( _varNum, 0.0 );
	d2N0.resize( _varNum, 0.0 );
}