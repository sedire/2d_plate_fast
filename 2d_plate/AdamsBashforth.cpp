#include "AdamsBashforth.h"

ABPrevPhis::ABPrevPhis()
{

}

ABPrevPhis::~ABPrevPhis()
{

}

AdamsBashforth::AdamsBashforth()
{
	prevPhis.resize( 5 );
}

AdamsBashforth::AdamsBashforth( int _eq_num )
{
	eq_num = _eq_num;
	prevPhis.resize( 5 );
	for( int i = 0; i < prevPhis.size(); ++i )
	{
		prevPhis[i].F.resize( _eq_num, 0.0 );
		prevPhis[i].F_1.resize( _eq_num, 0.0 );
		prevPhis[i].F_2.resize( _eq_num, 0.0 );
		prevPhis[i].F_3.resize( _eq_num, 0.0 );
	}
}

AdamsBashforth::~AdamsBashforth()
{

}

void AdamsBashforth::calc( const vector<PL_NUM>& A, const vector<PL_NUM>& f, PL_NUM dx, int n, int hom, vector<PL_NUM>* x )
{
	if( A.size() != eq_num * eq_num ||
		f.size() != eq_num ||
		x == 0 ||
		x->size() != eq_num ||
		dx <= 0.0 || n > 5 || n < 1 || ( hom != 0 && hom != 1 ) ) {
		cout << "Error in Adams-Bashforth calc: bad input\n";
		return;
	}

	n = n - 1;
	vector<PL_NUM> newPhi;
	newPhi.resize( 8, 0.0 );

	vector<PL_NUM> Np;		//correctors
	Np.resize( 8, 0.0 );

	for( int i = 0; i < eq_num; ++i )
	{
		Np[i] = (*x)[i];
		(*x)[i] += 1.0 / 24.0 * ( 55.0 * prevPhis[n].F[i] - 59.0 * prevPhis[n].F_1[i] + 37.0 * prevPhis[n].F_2[i] - 9.0 * prevPhis[n].F_3[i] );
	}
	
	for( int i = 0; i < eq_num; ++i )
	{
		for( int j = 0; j < eq_num; ++j )
		{
			newPhi[i] += A[i * eq_num + j] * (*x)[j] * dx;
		}
		newPhi[i] += hom * f[i] * dx;
	}

	for( int i = 0; i < eq_num; ++i )
	{
		(*x)[i] = Np[i] + 1.0 / 24.0 * ( 9 * newPhi[i] + 19.0 * prevPhis[n].F[i] - 5.0 * prevPhis[n].F_1[i] + prevPhis[n].F_2[i] );
	}

	for( int i = 0; i < eq_num; ++i )
	{
		prevPhis[n].F_3[i] = prevPhis[n].F_2[i];
		prevPhis[n].F_2[i] = prevPhis[n].F_1[i];
		prevPhis[n].F_1[i] = prevPhis[n].F[i];
		prevPhis[n].F[i] = newPhi[i];
	}
}
